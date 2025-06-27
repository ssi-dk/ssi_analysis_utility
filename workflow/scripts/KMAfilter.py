#!/usr/bin/env python3
import pandas as pd
import argparse
from typing import List, Dict
import os
import sys
import logging
import yaml

# Add utility paths
sys.path.insert(0, os.path.abspath("../scripts"))
from logging_utils import setup_logging

# ------------------------- Threshold Handling ------------------------- #

def load_thresholds_from_config(organism: str, config_dir="workflow/configs_species") -> Dict[str, List[int]]:
    """
    Load thresholds from species-specific YAML config.

    Args:
        organism (str): Organism name, e.g., "Clostridioides difficile"
        config_dir (str): Path to species config directory

    Returns:
        Dict[str, List[int]]: Dictionary of gene-specific thresholds
    """
    species_map = {
        "Clostridioides difficile": "C.diff",
        "Clostridium difficile": "C.diff",
        "C. difficile": "C.diff",
        "E. coli": "E.coli",
        "Escherichia coli": "E.coli",
        "E.coli": "E.coli",
    }

    species_key = species_map.get(organism.strip(), organism.strip())
    config_path = os.path.join(config_dir, f"{species_key}.yaml")

    logging.info(f"Loading thresholds from config for: {organism} from config file {config_path}")

    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Species config not found: {config_path}")
    
    with open(config_path, "r") as f:
        config = yaml.safe_load(f)

    try:
        thresholds = config["analyses_to_run"]["KMA_filter"]["thresholds"]
        return thresholds
    except KeyError:
        raise ValueError(f"Thresholds not defined in: {config_path}")

def get_threshold(template_name: str, thresholds: Dict[str, List[int]]) -> List[int]:
    """
    Match template name to threshold in dictionary.
    """
    for key in thresholds:
        if key in template_name:
            return thresholds[key]
    return thresholds["other"]

# ------------------------- KMA Processing Logic ------------------------- #

def process_kma_res(res_file: str, thresholds: Dict[str, List[int]]) -> pd.DataFrame:
    df = pd.read_csv(res_file, sep="\t")

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
        gene = parts[1] if len(parts) >= 2 else parts[0]

        if gene.lower() in known_genes:
            result[gene] = "Positive"
        else:
            other_genes.add(gene)

        if verbose_flag:
            verbose_parts.append(
                f"{gene}_{row['Template_length']}_{min(row['Template_Coverage'], 100.0):.2f}_{row['Template_Identity']:.2f}_{row['Depth']:.2f}"
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

# ------------------------- Main CLI Entry ------------------------- #

def main(args):
    setup_logging(args.log_dir, args.sample_id, "kma_filtering")

    try:
        thresholds = load_thresholds_from_config(args.organism)

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
    cols = output_df.columns.tolist()
    first_cols = ["sample_id", "organism"]
    output_df = output_df[first_cols + [col for col in cols if col not in first_cols]]

    sep = "\t" if args.suffix == "tsv" else "," if args.suffix == "csv" else " "
    output_df.to_csv(args.output, sep=sep, index=False)
    logging.info(f"Summary written to: {args.output}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter KMA .res file using gene thresholds from YAML config.")
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
