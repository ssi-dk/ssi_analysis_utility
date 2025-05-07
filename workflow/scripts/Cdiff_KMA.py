import pandas as pd
import argparse
import os
from typing import List, Dict
import pysam
import logging
from logging_utils import setup_logging
from thresholds import get_kma_thresholds_for_species, get_threshold

def process_res_file(res_file_path: str, thresholds: Dict[str, List[int]]) -> pd.DataFrame:
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
        (res_df["Query_Identity"] >= res_df["threshold"].apply(lambda x: x[1])) &
        (res_df["Depth"] >= res_df["threshold"].apply(lambda x: x[2]))
    ]
    return res_df_filtered

def summarize_single_sample(sample_name: str, res_path: str, verbose_flag: int = 1, thresholds: Dict[str, List[int]] = None) -> dict[str, str]:
    log_dir = "examples/Log"
    setup_logging(log_dir, sample_name, "kma_cdiff")

    NA_string = "-"
    output_data = {
        "tcdA": NA_string,
        "tcdB": NA_string,
        "tcdC": NA_string,
        "cdtAB": NA_string,
        "Other": NA_string
    }

    try:
        logging.info(f"Processing .res file: {res_path}")
        filtered_df = process_res_file(res_path, thresholds)
    except Exception as e:
        logging.error(f"Failed to process {res_path}: {e}")
        return output_data

    known_genes = {"tcdA", "tcdB", "tcdC", "cdtAB", "cdtA", "cdtB"}
    other_genes = set()

    for template in filtered_df["#Template"]:
        # Attempt to split like E. coli format
        gene_list = template.split("_")
        print(gene_list)

        gene = None

        if len(gene_list) >= 2:
            # Find first match to known gene
            for part in gene_list:
                if part.lower() in ["tcdA".lower(), "tcdB".lower(), "tcdC".lower()]:
                    gene = part
                    break
                elif "cdtA" in part and "cdtB" in part:
                    gene = "cdtAB"
                    break

        if gene == "tcdA":
            output_data["tcdA"] = "Positive"
        elif gene == "tcdB":
            output_data["tcdB"] = "Positive"
        elif gene == "tcdC":
            output_data["tcdC"] = "Positive"
        elif gene == "cdtAB":
            output_data["cdtAB"] = "Positive"
        else:
            # Fallback for any other unknown gene passing thresholds
            if len(gene_list) >= 2:
                possible_gene = gene_list[1]
                if possible_gene not in known_genes:
                    other_genes.add(possible_gene)

    if other_genes:
        output_data["Other"] = ";".join(sorted(other_genes))

    if verbose_flag == 1:
        verbose_parts = []
        for _, row in filtered_df.iterrows():
            template = row["#Template"]
            if "__" in template:
                parts = template.split("__")
            else:
                parts = template.split("_")

            if len(parts) >= 2:
                gene = parts[1]
                depth = row["Depth"]
                coverage = row["Template_Coverage"]
                identity = row["Query_Identity"]
                verbose_parts.append(f"{gene}_{depth:.2f}_{coverage:.2f}_{identity:.2f}")
        output_data["verbose"] = ";".join(verbose_parts)

    logging.info(f"Successfully processed sample: {sample_name}")
    return output_data

def load_toxin_coordinates(bed6_path: str) -> Dict[str, Dict[str, str | int]]:
    """
    Load coordinates of toxin genes from a BED6 file into a dictionary.

    Args:
        bed6_path (str): Path to the BED6 file.

    Returns:
        Dict[str, Dict[str, str | int]]: Dictionary mapping gene -> {contig, start, end, length, strand}
    """
    coords = {}
    try:
        bed_df = pd.read_csv(bed6_path, sep="\t", header=None,
                             names=["contig", "start", "end", "gene", "score", "strand"])

        for _, row in bed_df.iterrows():
            gene = row["gene"]
            coords[gene] = {
                "contig": row["contig"],
                "start": int(row["start"]),
                "end": int(row["end"]),
                "length": int(row["end"]) - int(row["start"]),
                "strand": row["strand"]
            }

    except Exception as e:
        logging.error(f"Failed to load BED6 file: {e}")
        raise

    return coords

def gene_pos_to_genomic(gene_name: str, pos_in_gene: int, coord_dict: Dict[str, Dict[str, str | int]]) -> tuple[str, int]:
    """
    Convert a position within a gene to a genomic coordinate, using strand-aware logic.

    Args:
        gene_name (str): Gene symbol (e.g. 'tcdC').
        pos_in_gene (int): 1-based position within the gene.
        coord_dict (dict): Output from load_toxin_coordinates().

    Returns:
        Tuple[str, int]: (contig, genomic_position)
    """
    if gene_name not in coord_dict:
        raise ValueError(f"Gene {gene_name} not found in coordinate dictionary.")

    info = coord_dict[gene_name]
    contig = info["contig"]
    strand = info["strand"]
    start = info["start"]
    end = info["end"]

    if strand == "+":
        genomic_pos = start + (pos_in_gene - 1)
    elif strand == "-":
        genomic_pos = end - (pos_in_gene - 1)
    else:
        raise ValueError(f"Invalid strand '{strand}' for gene {gene_name}.")

    return contig, genomic_pos

def main(args: argparse.Namespace) -> None:
    samplesheet_path = os.path.abspath(args.samplesheet)

    if not os.path.exists(samplesheet_path):
        logging.info(f"Samplesheet not found: {samplesheet_path}")
        exit(1)

    df = pd.read_csv(samplesheet_path, sep="\t")
    df.columns = df.columns.str.strip().str.lower()

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
    summary_columns = ["tcdA", "tcdB", "tcdC", "cdtAB", "Other", "verbose"]
    desired_columns = input_columns + summary_columns

    result_rows = []

    for idx, row in df.iterrows():
        row_dict = row.to_dict()
        sample = row_dict["sample_name"]
        res_path = f"examples/Results/{sample}/cdiffkmeraligner/{sample}.res"

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
    parser = argparse.ArgumentParser(description="Append C. difficile KMA .res summary info to a samplesheet.")
    parser.add_argument("-ss", "--samplesheet", required=True, help="Path to the input samplesheet .tsv file")
    parser.add_argument("-o", "--outputfile", default=None, help="Output filename. Default: per-sample output only")
    parser.add_argument("--suffix", choices=["tsv", "csv"], default="tsv", help="Output format. Default: tsv")
    parser.add_argument("--verbose", type=int, choices=[0, 1], default=1, help="Include verbose output. Default: 1")
    parser.add_argument("--split", type=int, choices=[0, 1], default=0, help="Split Illumina_read_files into Read1/Read2 (default: 0)")
    args = parser.parse_args()
    main(args)
