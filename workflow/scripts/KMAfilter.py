#!/usr/bin/env python3
import pandas as pd
import argparse
from typing import List, Dict
import os
import sys
import logging
sys.path.insert(0, os.path.abspath("../scripts"))
from logging_utils import setup_logging
from thresholds import get_kma_thresholds_for_species, get_threshold

def process_kma_res(res_file: str, thresholds: Dict[str, List[int]]) -> pd.DataFrame:
    try:
        df = pd.read_csv(res_file, sep="\t")
    except FileNotFoundError:
        raise FileNotFoundError(f"File not found: {res_file}")
    except pd.errors.EmptyDataError:
        raise ValueError(f"File is empty or malformed: {res_file}")

    required_cols = {"#Template", "Template_length", "Template_Coverage", "Template_Identity", "Depth"}
    if not required_cols.issubset(df.columns):
        raise ValueError(f"Missing required columns in {res_file}: {df.columns.tolist()}")

    df["threshold"] = df["#Template"].apply(lambda x: get_threshold(x, thresholds))
    filtered_df = df[
        (df["Template_Coverage"] >= df["threshold"].apply(lambda x: x[0])) &
        (df["Template_Identity"] >= df["threshold"].apply(lambda x: x[1])) &
        (df["Depth"] >= df["threshold"].apply(lambda x: x[2]))
    ]

    logging.info("Applied thresholds:")
    for idx, row in df.iterrows():
        logging.info(f"\t {row['#Template']}: Coverage>={row['threshold'][0]}, Identity>={row['threshold'][1]}, Depth>={row['threshold'][2]}")

    return filtered_df

def summarize_filtered_hits(sample_id: str, organism: str, filtered_df: pd.DataFrame, gene_list: List[str], verbose_flag: int = 1) -> Dict[str, str]:
    result = {gene: "Negative" for gene in gene_list}
    result["Other"] = "-"
    result["verbose"] = "-"
    known_genes = set(gene.lower() for gene in gene_list)
    other_genes = set()
    verbose_parts = []

    for _, row in filtered_df.iterrows():
        template = row["#Template"]
        parts = template.split("__") if "__" in template else template.split("_")
        if len(parts) >= 2:
            gene = parts[1]
        else:
            gene = parts[0]

        if gene.lower() in known_genes:
            result[gene] = "Positive"
        else:
            other_genes.add(gene)

        if verbose_flag:
            verbose_parts.append(
                f"{gene}_{row['Template_length']}_{min(row['Template_Coverage'],100.0):.2f}_{row['Template_Identity']:.2f}_{row['Depth']:.2f}"
            )

    if other_genes:
        result["Other"] = ";".join(sorted(other_genes))
    if verbose_parts:
        result["verbose"] = ";".join(verbose_parts)

    result.update({
        "sample_id": sample_id,
        "organism": organism
    })
    return result

def main(args):
    log_dir = args.log_dir 
    setup_logging(log_dir, args.sample_id, "kma_filtering")

    try:
        logging.info(f"Loading thresholds for organism: {args.organism}")
        thresholds = get_kma_thresholds_for_species(args.organism)

        #for gene, thresh in thresholds.items():
        #    logging.info(f"Threshold for {gene}: Coverage>={thresh[0]}, Identity>={thresh[1]}, Depth>={thresh[2]}")
        logging.info(f"Processing KMA .res file: {args.KMA_res}")

        filtered_df = process_kma_res(args.KMA_res, thresholds)
        result_dict = summarize_filtered_hits(args.sample_id, args.organism, filtered_df, args.Gene_list, verbose_flag=args.verbose)
    except Exception as e:
        logging.error(f"Error processing {args.sample_id}: {e}")
        result_dict = {gene: "-" for gene in args.Gene_list}
        result_dict.update({
            "Other": "-", "verbose": "-",
            "sample_id": args.sample_id,
            "organism": args.organism
        })

    output_df = pd.DataFrame([result_dict])

    # Reorder columns: sample_id, organism, then all others
    cols = output_df.columns.tolist()
    first_cols = ["sample_id", "organism"]
    remaining_cols = [col for col in cols if col not in first_cols]
    output_df = output_df[first_cols + remaining_cols]

    sep = "\t" if args.suffix == "tsv" else "," if args.suffix == "csv" else " "
    output_df.to_csv(args.output, sep=sep, index=False)
    logging.info(f"Summary written to: {args.output}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter KMA .res file using gene thresholds.")
    parser.add_argument("--KMA_res", required=True, help="Input KMA .res file")
    parser.add_argument("--Gene_list", nargs="+", required=True, help="List of gene symbols to track")
    parser.add_argument("--organism", required=True, help="Species name for threshold loading")
    parser.add_argument("--sample_id", required=True)
    parser.add_argument("--suffix", choices=["txt", "csv", "tsv"], default="tsv")
    parser.add_argument("--output", required=True)
    parser.add_argument("--verbose", type=int, choices=[0, 1], default=1)
    parser.add_argument("--log_dir", default="examples/Log")

    args = parser.parse_args()
    main(args)

# python KMAfilter.py --KMA_res ../../examples/Results/SRR10518319/Cdiff_KMA_Toxin/SRR10518319.res --Gene_list tcdA tcdB tcdC cdtAB --organism "Clostridioides difficile" --sample_id SRR10518319 --output ../../examples/Results/SRR10518319/KMA.tsv --suffix tsv --log_dir ../../examples/Results/SRR10518319
