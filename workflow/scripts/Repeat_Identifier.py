#!/usr/bin/env python3
import os
import sys
import argparse
import logging
import pandas as pd
from typing import Dict, List
from Bio import SeqIO
from Bio.Seq import Seq
import re

sys.path.insert(0, os.path.abspath("../scripts"))
from logging_utils import setup_logging

# ---------------------- FASTA and Type Parsers ---------------------- #

def parse_fasta(filepath: str) -> Dict[str, str]:
    return {record.id: str(record.seq) for record in SeqIO.parse(filepath, "fasta")}

def parse_repeat_types(type_file: str, fragments: Dict[str, str]) -> Dict[str, str]:
    type_map = {}
    with open(type_file) as f:
        for line in f:
            if ",\t" in line:
                type_id, pattern = line.strip().split(",\t")
                try:
                    sequence = ''.join(fragments[p] for p in pattern.split("-"))
                    type_map[type_id] = sequence
                except KeyError:
                    logging.warning(f"Fragment(s) missing for {type_id} in {type_file}")
    return type_map

def find_repeat_hits(sequence: str, type_map: Dict[str, str]) -> List[str]:
    hits = []
    for repeat_id, repeat_seq in type_map.items():
        pattern_seq = Seq(repeat_seq)
        if re.search(str(pattern_seq), sequence, re.IGNORECASE) or re.search(str(pattern_seq.reverse_complement()), sequence, re.IGNORECASE):
            hits.append(repeat_id)
    return hits

def load_combination_table(filepath: str) -> List[tuple]:
    with open(filepath) as f:
        return [tuple(line.strip().split("\t")) for line in f if len(line.strip().split("\t")) == 3]

def load_combination_table(filepath: str) -> List[tuple]:
    with open(filepath) as f:
        return [tuple(line.strip().split("\t")) for line in f if len(line.strip().split("\t")) == 3]

def infer_combined_hit(hits_a: List[str], hits_b: List[str], combo_table: List[tuple]) -> str:
    for combined, a, b in combo_table:
        if a in hits_a and b in hits_b:
            return combined
    return "Unknown"

# ------------------------- Main Typing Logic ------------------------- #

def run_repeat_typing(fasta_path: str, repeat_names: List[str], db_dir: str) -> Dict[str, str]:
    contigs = list(SeqIO.parse(fasta_path, "fasta"))
    logging.info(f"Found {len(contigs)} contigs in {fasta_path}")

    results = {}
    repeat_hits = {}

    # Step 1: Detect patterns for each repeat name
    for name in repeat_names:
        fasta_file = os.path.join(db_dir, f"{name}_repeat_sequences.fa")
        type_file = os.path.join(db_dir, f"{name}_repeat_types.txt")

        if not (os.path.exists(fasta_file) and os.path.exists(type_file)):
            logging.warning(f"Missing files for repeat type {name}. Skipping.")
            results[name] = "Unknown"
            continue

        fragments = parse_fasta(fasta_file)
        type_map = parse_repeat_types(type_file, fragments)
        hits = set()

        for record in contigs:
            hits.update(find_repeat_hits(str(record.seq), type_map))

        results[name] = ";".join(sorted(hits)) if hits else "Unknown"
        repeat_hits[name] = hits

    # Step 2: Check for any combination tables (e.g., inferred 3-column matchers)
    for file in os.listdir(db_dir):
        if file.endswith("_repeat_types.txt"):
            full_path = os.path.join(db_dir, file)
            with open(full_path) as f:
                header = f.readline().strip().split("\t")
            if len(header) == 3:
                # Detected a 3-column combination matcher
                combo_name = file.replace("_repeat_types.txt", "")
                combo_table = load_combination_table(full_path)
                _, name_a, name_b = header
                hits_a = repeat_hits.get(name_a, [])
                hits_b = repeat_hits.get(name_b, [])
                results[combo_name] = infer_combined_hit(hits_a, hits_b, combo_table)

    return results

# ---------------------------- Main CLI ---------------------------- #

def main(args):
    setup_logging(args.log_dir, args.sample_id, "repeat_typing")

    try:
        results = run_repeat_typing(
            fasta_path=args.fasta,
            repeat_names=args.repeats,
            db_dir=args.db_dir
        )
    except Exception as e:
        logging.error(f"Typing failed for {args.sample_id}: {e}")
        results = {name: "-" for name in args.repeats}

    results["sample_id"] = args.sample_id
    df = pd.DataFrame([results])
    df = df[["sample_id"] + [col for col in df.columns if col != "sample_id"]]

    sep = "\t" if args.suffix == "tsv" else "," if args.suffix == "csv" else " "
    df.to_csv(args.output, sep=sep, index=False)
    logging.info(f"Repeat typing results written to {args.output}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generalized repeat typing from assemblies using pattern matching.")
    parser.add_argument("--sample_id", required=True, help="Sample identifier")
    parser.add_argument("--fasta", required=True, help="Path to assembly FASTA file")
    parser.add_argument("--repeats", nargs="+", required=True, help="List of repeat identifiers (must match filename prefixes)")
    parser.add_argument("--db_dir", default="resources/Clostridioides_difficile_db/TRST", help="Directory containing repeat sequence and type files")
    parser.add_argument("--output", required=True, help="Output file path")
    parser.add_argument("--suffix", choices=["tsv", "csv", "txt"], default="tsv", help="Output file format")
    parser.add_argument("--log_dir", default="examples/Log", help="Logging directory")
    args = parser.parse_args()
    main(args)
