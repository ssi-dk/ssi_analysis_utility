import pandas as pd
import argparse
from pathlib import Path
import logging 
import sys
import os

# pd.set_option('display.max_columns', None)  # or 1000
# pd.set_option('display.max_rows', None)  # or 1000
# pd.set_option('display.max_colwidth', None)  # or 199




if __name__ == '__main__':

    # Parse command line arguments
    parser = argparse.ArgumentParser(description='RAR-output_handler.py')
    group = parser.add_mutually_exclusive_group(required=True)

    group.add_argument("-f", "--filelist", 
                        help='Path to file list containig the path \
                              to the samples', 
                        )
    group.add_argument("-s", "--sampleID", 
                        help='Path to Samples Folders', 
                        nargs = "+")
    parser.add_argument("-o", "--outname", 
                        help='Output tsv report name', 
                        default="report.tsv")
    args = parser.parse_args()

    output_to_scan = {"plasmidfinder":"results_tab.txt",
                      "resfinder":"ResFinder_results_tab.txt"}
    
    if args.sampleID:
        path_samples = args.sampleID 
    elif args.filelist:
        path_samples = [line.rstrip() for line in open(args.filelist)]

    final_plasmid_df = pd.DataFrame()    
    final_resfinder_df = pd.DataFrame() 
    for sample in path_samples:
        # Get the sample ID from the path
        sample_ID = os.path.basename(sample)

        # Check if the plasmidfinder and resfinder output exists 
        plasmidfinder_out = Path(sample + "/plasmidfinder/results_tab.tsv")
        resfinder_out = Path(sample + "/resfinder/ResFinder_results_tab.txt")
        
        # If both files are not available then exit
        if not plasmidfinder_out.is_file() and not resfinder_out.is_file():
            sys.exit("Resfinder and Plasmidfinder output are missing check, please check your samples")
        
        try:
            plasmid_df = pd.read_csv(plasmidfinder_out,sep="\t")
            resfinder_df = pd.read_csv(resfinder_out,sep="\t")
            
            plasmid_df["SampleID"] = sample_ID
            resfinder_df["SampleID"] = sample_ID
            # If either is empty fill with na
            if plasmid_df.empty:
                plasmid_df = pd.DataFrame(columns=plasmid_df.columns)  # keeps the same columns
                plasmid_df.loc[0] = ["NA"] * len(plasmid_df.columns)
                plasmid_df["SampleID"] = sample_ID
            elif resfinder_df.empty:
                resfinder_df = pd.DataFrame(columns=resfinder_df.columns)  # keeps the same columns
                resfinder_df.loc[0] = ["NA"] * len(resfinder_df.columns)
                resfinder_df["SampleID"] = sample_ID
            
            final_plasmid_df = pd.concat([final_plasmid_df, plasmid_df ])
            final_resfinder_df = pd.concat([final_resfinder_df, resfinder_df ])

            final_plasmid_df['Identity'] = pd.to_numeric(final_plasmid_df['Identity'], errors='coerce')
            final_plasmid_df = final_plasmid_df[(final_plasmid_df['Identity'] > 80) | (final_plasmid_df['Identity'].isna())]
            final_resfinder_df = final_resfinder_df[(final_resfinder_df['Identity'] >= 95) & (final_resfinder_df['Coverage'] >= 95)]


        except:
            raise("Cannot load results")

    
    final_plasmid_df = final_plasmid_df[["SampleID","Plasmid"]]
    
    # Create a list of unique sample IDs
    sample_ids = final_resfinder_df["SampleID"].unique()

    # Get all unique phenotypes (which will be used as columns)
    phenotype_columns = final_resfinder_df["Phenotype"].unique().tolist()
    

    # Initialize dictionary with SampleID and phenotype columns
    data_dict = {col: [] for col in ["SampleID"] + phenotype_columns }
   
    # Fill the dictionary row by row
    for sample in sample_ids:
        sample_df = final_resfinder_df[final_resfinder_df["SampleID"] == sample]
        
        data_dict["SampleID"].append(sample)

        for phenotype in phenotype_columns:
            genes = sample_df[sample_df["Phenotype"] == phenotype]["Resistance gene"].unique()
            data_dict[phenotype].append(", ".join(genes) if len(genes) > 0 else "NA")

    df = pd.DataFrame.from_dict(data_dict)
    final_df = pd.merge(df, final_plasmid_df, on='SampleID')
