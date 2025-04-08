import pandas as pd
import argparse
from pathlib import Path
import logging
import os
import sys


# Set up logging configuration for informative output and debugging purposes
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


def parse_arguments():
    """
    Parse command-line arguments to get input sample paths and output filename
    """
    parser = argparse.ArgumentParser(
        description='RAR-output_handler.py: Process ResFinder and PlasmidFinder outputs into a summary report'
    )

    # User must supply either a file list or a list of sample paths, but not both
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "-f", "--filelist",
        help="Path to a text file containing a list of sample folder paths (one per line)"
    )
    group.add_argument(
        "-s", "--sampleID",
        help="One or more paths to individual sample folders",
        nargs="+"
    )

    # Output filename, default is 'report.tsv'
    parser.add_argument(
        "-o", "--outname",
        help="Name of the output TSV report file",
        default="report.tsv"
    )

    return parser.parse_args()


def load_or_na(filepath, 
               sample_id):
    """
    Load a data file (TSV) or return a DataFrame filled with "NA" 
    if file is empty or unreadable
    Parameters
    -------------
    filepath : pathlib.PosixPath object
        Object containing the path to the file
    sample_id : str object
        String containing the name of the sample ID
    """

    try:
        df = pd.read_csv(filepath, 
                         sep="\t") # read the results

        # If the file loads but is empty, create a single NA-filled row
        if df.empty:
            df = pd.DataFrame(columns = df.columns)
            df.loc[0] = ["NA"] * len(df.columns)

        # Add SampleID to keep track of origin
        df["SampleID"] = sample_id
        return df

    except Exception as e:
        # Warn User that something is wrong with files
        logger.warning(f"Failed to load file {filepath}: {e}") 
        # Return an empty Dataframe 
        df_na = pd.DataFrame([["NA"]], columns=["SampleID"]).assign(SampleID=sample_id)
        return df_na


def collect_data(sample_paths):
    """
    Collect and combine all ResFinder and PlasmidFinder data from each sample
    Parameters
    -------------
    sample_paths : list object
        list containing all the paths to the samples
    """
    

    plasmid_dfs = []
    resfinder_dfs = []
    
    # Loop the paths
    for sample_path in sample_paths:
        sample_id = os.path.basename(sample_path)  # Use folder name as SampleID

        # Construct paths to result files
        plasmid_file = Path(sample_path) / "plasmidfinder" / "results_tab.tsv"
        resfinder_file = Path(sample_path) / "resfinder" / "ResFinder_results_tab.txt"

        # If both files are missing, skip this sample with a warning
        if not plasmid_file.exists() and not resfinder_file.exists():
            logger.warning(f"Missing both result files for sample {sample_id}, skipping.")
            continue

        # Load data or substitute with "NA" row if not available
        plasmid_df = load_or_na(plasmid_file, 
                                sample_id)
        resfinder_df = load_or_na(resfinder_file, 
                                  sample_id)

        # Collect for merging later
        plasmid_dfs.append(plasmid_df)
        resfinder_dfs.append(resfinder_df)

    # Combine all dataframes across samples into one per type
    plasmid_concat = pd.concat(plasmid_dfs, 
                               ignore_index=True)
    resfinder_concat = pd.concat(resfinder_dfs,
                                 ignore_index=True)
    return plasmid_concat, resfinder_concat


def filter_plasmidfinder_data(plasmid_df):
    """
    Filter raw result data based on set tresholds
    Parameters
    -------------
    plasmid_df : pandas.DataFrame object
        Pandas dataframe containing the results from Plasmidfinder
    resfinder_df : pandas.DataFrame object
        Pandas dataframe containing the results from Resfinder
    """
    
    
    # Convert 'Identity' to numeric and check the name of the plasmid
    # if the name starts with Col then the ID threshold is set to 80 
    # otherwise set it at 95%
    
    # Ensure the 'Identity' column is numeric (non-numeric values become NaN)
    plasmid_df["Identity"] = pd.to_numeric(plasmid_df.get("Identity", 
                                                          pd.Series()), 
                                           errors="coerce")
    
    # Create two boolean masks:
    # - One for plasmids that start with 'Col' (e.g., ColRNAI, ColE1, etc.)
    # - Another for all other plasmid names
    is_col_plasmid = plasmid_df['Plasmid'].str.startswith("Col", na=False)
    is_non_col_plasmid = ~is_col_plasmid

    # Apply filtering:
    # - Keep Col plasmids if Identity > 80
    # - Keep non-Col plasmids if Identity > 95
    filtered_plasmid_df = plasmid_df[(is_col_plasmid & (plasmid_df['Identity'] >= 80)) |
                                     (is_non_col_plasmid & (plasmid_df['Identity'] >= 95))]

    
    return filtered_plasmid_df

def filter_resfinder_data(resfinder_df,
                          coverage_threshold,
                          identity_threshold):

    # Filter ResFinder hits where both Identity and Coverage are >= 95
    resfinder_df["Identity"] = pd.to_numeric(resfinder_df.get("Identity", 
                                                              pd.Series()), 
                                            errors="coerce")
    resfinder_df["Coverage"] = pd.to_numeric(resfinder_df.get("Coverage", 
                                                              pd.Series()), 
                                            errors="coerce")
    resfinder_df = resfinder_df[(resfinder_df["Identity"] >= identity_threshold) & 
                                (resfinder_df["Coverage"] >= coverage_threshold)]

    return resfinder_df
#IncI1-I(Alpha)


def build_phenotype_matrix(resfinder_df):
    """
    Transform ResFinder results for the final report:
    One row per sample, one column per phenotype, containing resistance genes
    Parameters
    -------------

    resfinder_df : pandas.DataFrame object
        Pandas dataframe containing the filtered results from Resfinder
    """    

    sample_ids = resfinder_df["SampleID"].unique()
    phenotype_columns = resfinder_df["Phenotype"].dropna().unique().tolist()

    records = []

    for sample_id in sample_ids:
        sample_data = {"SampleID": sample_id}
        sample_df = resfinder_df[resfinder_df["SampleID"] == sample_id]

        for phenotype in phenotype_columns:
            # Collect all unique resistance genes for the phenotype
            genes = sample_df[sample_df["Phenotype"] == phenotype]["Resistance gene"].dropna().unique()
            sample_data[phenotype] = ", ".join(genes) if genes.size else "NA"

        records.append(sample_data)
    records = pd.DataFrame(records)
    return records


def main():
    """
    Main pipeline logic: load, process, transform, and export the report
    """   
    args = parse_arguments()

    # Load paths from file or command-line list
    sample_paths = args.sampleID if args.sampleID else [line.strip() for line in open(args.filelist)]

    # Step 1: Load all relevant results
    plasmid_df, resfinder_df = collect_data(sample_paths)

    if plasmid_df.empty and resfinder_df.empty:
        logger.error("No valid data found. Exiting.")
        sys.exit(1)

    # Step 2: Apply quality filters to the data
    plasmid_df= filter_plasmidfinder_data(plasmid_df)
    resfinder_df = filter_resfinder_data(resfinder_df,
                                         95,
                                         95)

    # Step 3: Retain only relevant plasmid columns
    if "Plasmid" in plasmid_df.columns:
        plasmid_df = plasmid_df[["SampleID", 
                                 "Plasmid"]]
    else:
        plasmid_df = pd.DataFrame(columns=["SampleID", 
                                           "Plasmid"])

    # Step 4: Create a pivoted phenotype matrix (resistance gene per sample)
    phenotype_df = build_phenotype_matrix(resfinder_df)

    # Step 5: Merge resistance and plasmid data on SampleID
    if not plasmid_df.empty:
        
        final_df = pd.merge(phenotype_df, 
                            plasmid_df, 
                            on="SampleID",
                            how="left") 
    else:
        final_df = phenotype_df.copy()

    # Step 6: Output final summary table to TSV
    final_df.to_csv(args.outname, 
                    sep="\t", 
                    index=False)
    logger.info(f"Report saved to {args.outname}")



if __name__ == "__main__":
    main()

