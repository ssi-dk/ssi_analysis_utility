import pandas as pd
import argparse
import os
from typing import List, Dict
import logging
from logging_utils import setup_logging
from thresholds import get_kma_thresholds_for_species, get_threshold

def process_res_file(res_file_path: str, thresholds: Dict[str, List[int]]) -> pd.DataFrame:
    """
    Reads and filters a KMA .res file based on predefined thresholds.

    Args:
        res_file_path (str): Path to the .res file.
        thresholds (Dict[str, List[int]]): Gene-specific thresholds.

    Returns:
        pd.DataFrame: Filtered results DataFrame.
    """
    try:
        res_df = pd.read_csv(res_file_path, sep="\t")
    except FileNotFoundError:
        raise FileNotFoundError(f"File not found: {res_file_path}")
    except pd.errors.EmptyDataError:
        raise ValueError(f"File is empty or not properly formatted: {res_file_path}")

    required_columns = {"#Template", "Template_Coverage", "Query_Identity", "Depth"}
    if not required_columns.issubset(res_df.columns):
        raise ValueError(f"Missing expected columns in {res_file_path}. Found: {', '.join(res_df.columns)}")

    res_df["threshold"] = res_df["#Template"].apply(lambda x: get_threshold(x, thresholds))
    res_df_filtered = res_df[
        (res_df["Template_Coverage"] >= res_df["threshold"].apply(lambda x: x[0])) &
        (res_df["Query_Identity"] >= res_df["threshold"].apply(lambda x: x[1]))
    ]
    return res_df_filtered

def summarize_single_sample(sample_name: str, res_path: str, verbose_flag: int = 1, thresholds: Dict[str, List[int]] = None) -> dict[str, str]:
    """
    Processes a single sample’s .res file and returns a summary dictionary.

    Args:
        sample_name (str): Sample identifier.
        res_path (str): Path to the sample's .res file.
        verbose_flag (int, optional): Include verbose info if set to 1. Default is 1.
        thresholds (Dict[str, List[int]]): Gene threshold dictionary.

    Returns:
        Dict[str, str]: Summary values extracted from the .res file.
    """
    log_dir = "examples/Log"
    setup_logging(log_dir, sample_name,"kma_ecoli")

    NA_string = "-"
    output_data = {
        "stx": NA_string,
        "OH": NA_string, "wzx": NA_string, "wzy": NA_string, "wzt": NA_string, "wzm": NA_string,
        "eae": NA_string, "ehxA": NA_string,
        "Other": NA_string
    }

    try:
        logging.info(f"Processing .res file: {res_path}")
        filtered_df = process_res_file(res_path, thresholds)
    except Exception as e:
        logging.error(f"Failed to process {res_path}: {e}")
        return output_data

    gene_map = {
        "wzx": "wzx", "wzy": "wzy", "wzt": "wzt", "wzm": "wzm",
        "eae": "eae", "ehxA": "ehxA"
    }
    toxin = "stx"
    stx_alleles = set()
    fli = NA_string
    fliC = NA_string

    for template in filtered_df["#Template"]:
        gene_list = template.split("__")
        if len(gene_list) < 3:
            continue

        gene, allele = gene_list[1], gene_list[2]

        if gene in ["eae", "ehxA"]:
            output_data[gene] = "Positive"
        elif gene in gene_map:
            output_data[gene] = allele
        elif gene == "fliC":
            fliC = allele
        elif gene == "fli":
            fli = allele
        elif gene.startswith(toxin):
            stx_alleles.add(allele)
        elif gene not in thresholds:
            output_data["Other"] = allele

    if stx_alleles:
        output_data[toxin] = ";".join(sorted(stx_alleles))

    # Determine O-type
    wzx, wzy, wzt, wzm = output_data["wzx"], output_data["wzy"], output_data["wzt"], output_data["wzm"]
    Otype = "-"
    if wzx != NA_string and wzy != NA_string and wzt == NA_string and wzm == NA_string:
        if wzx == wzy:
            Otype = wzx
            output_data["wzx"] = output_data["wzy"] = NA_string
    elif wzt != NA_string and wzm != NA_string and wzx == NA_string and wzy == NA_string:
        if wzt == wzm:
            Otype = wzt
            output_data["wzt"] = output_data["wzm"] = NA_string

    # Determine H-type
    Htype = fli if fli != NA_string else fliC
    output_data["OH"] = f"{Otype};{Htype}"

    if verbose_flag == 1:
        verbose_parts = []
        for _, row in filtered_df.iterrows():
            parts = row["#Template"].split("__")
            if len(parts) >= 3:
                gene, allele = parts[1], parts[2]
                depth = row["Depth"]
                coverage = row["Template_Coverage"]
                identity = row["Query_Identity"]
                verbose_parts.append(f"{gene}_{allele}_{depth:.2f}_{coverage:.2f}_{identity:.2f}")
        output_data["verbose"] = ";".join(verbose_parts)

    logging.info(f"Successfully processed sample: {sample_name}")
    return output_data

def main(args: argparse.Namespace) -> None:
    """
    Main function to process a samplesheet and write either merged or per-sample output files.

    Args:
        args (argparse.Namespace): Parsed command-line arguments.
    """
    samplesheet_path = os.path.abspath(args.samplesheet)

    if not os.path.exists(samplesheet_path):
        logging.info(f"Samplesheet not found: {samplesheet_path}")
        exit(1)

    # Load and normalize column names
    df = pd.read_csv(samplesheet_path, sep="\t")
    df.columns = df.columns.str.strip().str.lower()

    # Determine which columns to include based on presence/split
    if args.split == 1:
        input_columns = ["sample_name", "read1", "read2"]
        if "illumina_read_files" in df.columns:
            logging.info("Splitting illumina_read_files column into read1 and read2...")
            df[["read1", "read2"]] = df["illumina_read_files"].str.split(",", expand=True)
    elif "read1" in df.columns and "read2" in df.columns:
        input_columns = ["sample_name", "read1", "read2"]
    elif "illumina_read_files" in df.columns:
        input_columns = ["sample_name", "illumina_read_files"]
    else:
        raise ValueError("No suitable Illumina read columns found (read1/read2 or illumina_read_files).")

    input_columns += ["nanopore_read_file", "assembly_file", "organism", "variant", "notes"]
    summary_columns = ["stx", "OH", "wzx", "wzy", "wzt", "wzm", "eae", "ehxA", "Other", "verbose"]
    desired_columns = input_columns + summary_columns

    result_rows = []

    for idx, row in df.iterrows():
        row_dict = row.to_dict()
        sample = row_dict["sample_name"]
        res_path = f"examples/Results/{sample}/kmeraligner/{sample}.res"

        try:
            thresholds = get_kma_thresholds_for_species(row_dict["organism"])
        except ValueError as e:
            logging.error(f"Sample '{sample}' skipped: {e}")
            continue

        result_dict = summarize_single_sample(sample, res_path, verbose_flag=args.verbose, thresholds=thresholds)
        row_dict.update(result_dict)

        for col in desired_columns:
            row_dict.setdefault(col, "-")

        result_rows.append(row_dict)

        if args.outputfile is None:
            output_path = os.path.join("examples", "Results", sample, "kma", f"{sample}_kma.{args.suffix}")
            os.makedirs(os.path.dirname(output_path), exist_ok=True)
            df_row = pd.DataFrame([row_dict])[desired_columns]
            df_row.to_csv(output_path, sep="\t" if args.suffix == "tsv" else ",", index=False)
            logging.info(f"Per-sample summary written to: {output_path}")

    if args.outputfile:
        for row_dict in result_rows:
            for col in desired_columns:
                row_dict.setdefault(col, "-")
        out_df = pd.DataFrame(result_rows)[desired_columns]
        out_df.to_csv(args.outputfile, sep="\t" if args.suffix == "tsv" else ",", index=False)
        logging.info(f"Merged summary written to: {args.outputfile}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Append KMA .res summary info to a samplesheet.")
    parser.add_argument("-ss", "--samplesheet", required=True, help="Path to the input samplesheet .tsv file")
    parser.add_argument("-o", "--outputfile", default=None, help="Output filename. Default: per-sample output only")
    parser.add_argument("--suffix", choices=["tsv", "csv"], default="tsv", help="Output format. Default: tsv")
    parser.add_argument("--verbose", type=int, choices=[0, 1], default=1, help="Include verbose output. Default: 1")
    parser.add_argument("--split", type=int, choices=[0, 1], default=0, help="Split Illumina_read_files into Read1/Read2 (default: 0)")
    args = parser.parse_args()
    main(args)