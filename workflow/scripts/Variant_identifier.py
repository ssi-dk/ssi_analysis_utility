#!/usr/bin/env python3
import argparse
import os
import sys
import logging
import pandas as pd
import pysam
import yaml
import re
from typing import Dict, List, Tuple

# Add utility path
sys.path.insert(0, os.path.abspath("../scripts"))
from logging_utils import setup_logging

# ========================= Helper Functions ========================= #

def load_variant_config(organism: str, config_dir="workflow/configs_species") -> Tuple[dict, dict, List[str], dict]:
    species_map = {
        "Clostridioides difficile": "C.diff",
        "Clostridium difficile": "C.diff",
        "C. difficile": "C.diff",
    }
    species_key = species_map.get(organism.strip(), organism.strip())
    config_path = os.path.join(config_dir, f"{species_key}.yaml")

    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Config file not found: {config_path}")

    with open(config_path, "r") as f:
        config = yaml.safe_load(f)

    try:
        variant_cfg = config["analyses_to_run"]["Variant_detection"]
        return (
            variant_cfg["deletion_regions"],
            variant_cfg["variant_gt_thresholds"],
            variant_cfg["variant_labels"],
            variant_cfg["snp_info"]
        )
    except KeyError as e:
        raise ValueError(f"Missing Variant_detection fields in config {config_path}: {e}")

def load_bed_coordinates(bed_path: str) -> Dict[str, Dict[str, str | int]]:
    coords = {}
    df = pd.read_csv(bed_path, sep="\t", header=None, names=["contig", "start", "end", "gene", "score", "strand"])
    for _, row in df.iterrows():
        coords[row["gene"]] = {
            "contig": row["contig"],
            "start": int(row["start"]),
            "end": int(row["end"]),
            "length": int(row["end"]) - int(row["start"]),
            "strand": row["strand"]
        }
    return coords

def resolve_contig_for_gene(gene: str, available: List[str]) -> str:
    if gene in available:
        return gene

    pattern = re.compile(fr"(?:^|_){re.escape(gene)}(?:_|$)")
    for contig in available:
        if pattern.search(contig):
            logging.warning(f"Contig for gene '{gene}' not found directly; using fallback match: '{contig}'")
            return contig

    raise ValueError(f"No contig found for gene '{gene}' in available contigs.")

def extract_start_from_contig_name(contig_name: str) -> int:
    match = re.search(r"_(\d{5,})_\d{5,}$", contig_name)
    if not match:
        raise ValueError(f"Could not extract embedded start coordinate from contig: {contig_name}")
    return int(match.group(1))

def gene_pos_to_genomic(gene_name: str, pos_in_gene: int, coord_dict: Dict[str, Dict[str, str | int]], contig_list: List[str]) -> Tuple[str, int]:
    info = coord_dict[gene_name]
    base_contig = info["contig"]
    strand = info["strand"]
    gene_len = info["length"]

    embedded_contig = resolve_contig_for_gene(gene_name, contig_list)
    embedded_start = extract_start_from_contig_name(embedded_contig)

    if strand == "+":
        pos = embedded_start + (pos_in_gene - 1)
    elif strand == "-":
        pos = embedded_start + (gene_len - pos_in_gene)
    else:
        raise ValueError(f"Invalid strand: {strand}")

    return embedded_contig, pos

def check_snp_variant(bcf_path: str, contig: str, pos: int, ref: str, alt: str, window=20) -> str:
    try:
        bcf = pysam.VariantFile(bcf_path)
        for rec in bcf.fetch(contig, pos-window, pos+window):
            if rec.pos == pos:
                if rec.ref == ref and alt in rec.alts:
                    return f"{ref}>{alt}"
                elif rec.ref != ref:
                    return "other"
            elif rec.pos < pos and any(len(a) < len(rec.ref) for a in rec.alts if a):
                return "del"
        return "wt"
    except Exception as e:
        logging.warning(f"SNP check failed: {e}")
        return "-"

def evaluate_consensus_n_percent(fsa_path: str, contig: str, start: int, end: int, gene: str) -> str:
    try:
        fasta = pysam.FastaFile(fsa_path)
        contig = resolve_contig_for_gene(gene, fasta.references)
        seq = fasta.fetch(contig, start, end)
        n_percent = (seq.upper().count("N") / len(seq)) * 100 if len(seq) > 0 else 0
        return f"{seq}_{n_percent:.2f}"
    except Exception as e:
        logging.warning(f"Consensus check failed: {e}")
        return "-"

# ========================= Main Entry ========================= #

def main(args):
    setup_logging(args.log_dir, args.sample_id, "variant_checker")
    logging.info("===============================================")
    logging.info(f"Logging started for sample {args.sample_id}")
    logging.info(f"Organism: {args.organism}")
    logging.info(f"Files:")
    logging.info(f"  BED file: {args.bed}")
    logging.info(f"  Consensus FASTA: {args.consensus}")
    logging.info(f"  Mpileup BCF: {args.mpileup}")
    logging.info(f"  Variant call BCF: {args.call}")
    logging.info("===============================================")

    output_data = {
        "sample_id": args.sample_id,
        "organism": args.organism
    }

    try:
        deletion_regions, gt_thresholds, variant_labels, snp_info = load_variant_config(args.organism, args.config_dir)
        coords = load_bed_coordinates(args.bed)

        for label in variant_labels:
            output_data[label] = "-"

        contig_list = list(pysam.VariantFile(args.mpileup).header.contigs)

        for snp_id, (pos_in_gene, ref, alt, gene) in snp_info.items():
            contig, snp_pos = gene_pos_to_genomic(gene, pos_in_gene, coords, contig_list)
            embedded_start = extract_start_from_contig_name(contig)

            logging.info(f"Checking SNP: gene={gene}, position_in_gene={pos_in_gene}, ref={ref}, alt={alt}")
            logging.info(f"  → Resolved contig: {contig}")
            logging.info(f"  → Converted coordinate in contig (BCF-relative): {snp_pos}")
            logging.info(f"  → Full genomic coordinate: {snp_pos} (offset from {embedded_start})")

            result = check_snp_variant(args.mpileup, contig, snp_pos, ref, alt)
            label = f"{gene}{pos_in_gene}"
            if label in output_data:
                output_data[label] = result
            else:
                output_data[f"{label}_snp"] = result

        bcf = pysam.VariantFile(args.call)
        contig_list = list(bcf.header.contigs)

        for rec in bcf.fetch():
            if rec.qual != 0 or not rec.alts or len(rec.ref) <= len(rec.alts[0]):
                continue

            del_len = len(rec.ref) - len(rec.alts[0])
            del_start = rec.pos

            for key, val in deletion_regions.items():
                start, end, n_thresh, gene = val
                gt = gt_thresholds[key]
                if gt[-1] != gene:
                    continue
                expected_len = int(key.split("_")[-1])
                if del_len != expected_len or not (start <= del_start <= end):
                    continue

                contig = resolve_contig_for_gene(gene, contig_list)
                g_start = gene_pos_to_genomic(gene, start, coords, contig_list)[1]
                g_end = gene_pos_to_genomic(gene, end, coords, contig_list)[1]
                embedded_start = extract_start_from_contig_name(contig)

                logging.info(f"Checking deletion: gene={gene}, region=({start}-{end}), expected_len={expected_len}")
                logging.info(f"  → Resolved contig: {contig}, observed_del_start={del_start}, observed_del_len={del_len}")
                logging.info(f"  → Converted coordinates in contig (BCF-relative): start={g_start}, end={g_end}")
                logging.info(f"  → Full genomic coordinates: start={g_start}, end={g_end} (offset from {embedded_start})")

                IMF = rec.info.get("IMF", 0)
                IDV = rec.info.get("IDV", 0)
                DP = rec.info.get("DP", 0)

                logging.info(f"  IMF={IMF}, IDV={IDV}, DP={DP} | thresholds={gt[:3]}")

                if IMF >= gt[0] and IDV >= gt[1] and DP >= gt[2]:
                    label = f"{gene}del"
                    output_data[label] = key
                    seq_info = evaluate_consensus_n_percent(args.consensus, contig, min(g_start, g_end), max(g_start, g_end), gene)
                    output_data[f"{label}_consensus"] = seq_info
                    try:
                        n_pct = float(seq_info.strip().split("_")[-1])
                        if n_pct >= n_thresh:
                            output_data[label] = key.replace("ambiguous", "likely")
                    except:
                        pass
                    break

    except Exception as e:
        logging.error(f"Variant analysis failed: {e}")

    df = pd.DataFrame([output_data])
    sep = "\t" if args.suffix == "tsv" else "," if args.suffix == "csv" else " "
    df.to_csv(args.output, sep=sep, index=False)
    logging.info(f"Output written to: {args.output}")

# ========================= CLI Parser ========================= #

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generalized variant and deletion checker.")
    parser.add_argument("--sample_id", required=True)
    parser.add_argument("--organism", required=True)
    parser.add_argument("--consensus", required=True)
    parser.add_argument("--bed", required=True)
    parser.add_argument("--mpileup", required=True)
    parser.add_argument("--call", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--suffix", choices=["tsv", "csv"], default="tsv")
    parser.add_argument("--log_dir", default="examples/Log")
    parser.add_argument("--config_dir", default="workflow/configs_species")
    args = parser.parse_args()
    main(args)