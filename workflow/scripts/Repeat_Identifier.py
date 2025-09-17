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
    """
    Parse a FASTA file and return a dictionary mapping sequence IDs to sequences.

    Args:
        filepath (str): Path to FASTA file.

    Returns:
        Dict[str, str]: Dictionary of sequence ID to sequence string.
    """
    return {record.id: str(record.seq) for record in SeqIO.parse(filepath, "fasta")}

def parse_repeat_types(ref_txt: str, fragments: Dict[str, str]) -> Dict[str, str]:
    """
    Parse a repeat type file and map type IDs to concatenated fragment sequences.

    Args:
        ref_txt (str): Path to repeat types file.
        fragments (Dict[str, str]): Dictionary of sequence fragments.

    Returns:
        Dict[str, str]: Mapping from repeat ID to full sequence pattern.
    """
    type_map = {}
    with open(ref_txt) as f:
        logging.debug(f"Loaded {len(type_map)} repeat patterns from {ref_txt}")

        for line in f:
            if ",\t" in line:
                type_id, pattern = line.strip().split(",\t")
                try:
                    sequence = ''.join(fragments[p] for p in pattern.split("-"))
                    type_map[type_id] = sequence
                except KeyError:
                    logging.warning(f"Fragment(s) missing for {type_id} in {ref_txt}")
    return type_map

def find_repeat_hits(sequence: str, id: str, type_map: Dict[str, str]) -> List[str]:
    """
    Scan a sequence for any repeat patterns defined in type_map.

    Args:
        sequence (str): Full DNA sequence to scan.
        type_map (Dict[str, str]): Dictionary mapping repeat ID to pattern sequence.

    Returns:
        List[str]: List of repeat IDs that were matched.
    """
    hits = []
    for repeat_id, repeat_seq in type_map.items():
        pattern_seq = Seq(repeat_seq)
        if re.search(str(pattern_seq), sequence, re.IGNORECASE) or re.search(str(pattern_seq.reverse_complement()), sequence, re.IGNORECASE):
            logging.info(f"found match between {id} and {str(repeat_id)}")
            hits.append(repeat_id)
    return hits

# ---------------------- Combination Matcher ---------------------- #
def match_combination_from_table(ref_combo: str, repeat_hits: Dict[str, List[str]], repeat_names: List[str]) -> str:
    """
    Match observed repeat hits against a combination lookup table to assign a combo ID.

    Args:
        ref_combo (str): Path to combo lookup file.
        repeat_hits (Dict[str, List[str]]): Observed repeat ID hits per repeat name.
        repeat_names (List[str]): List of repeat names in expected order.

    Returns:
        str: Matched combination ID or "Unknown" if no match.
    """
    with open(ref_combo) as f:
        lines = [line.strip() for line in f if line.strip()]

    if len(lines) < 2:
        logging.warning(f"Combination file {ref_combo} is empty or malformed.")
        return "Unknown"

    header = lines[0].split("\t")
    data_rows = [line.split("\t") for line in lines[1:] if len(line.split("\t")) == len(header)]

    if len(header) < 3: #need at least three columns, the combined type, and then two or more repeats to be combined.
        logging.warning(f"Combination file {ref_combo} must have at least 3 columns.")
        return "Unknown"

    # Skip first column (e.g. ST_COMB) and assume remaining columns align with repeat_names
    combo_fields = header[1:]

    if len(combo_fields) != len(repeat_names):
        logging.warning(f"Mismatch between combo fields {combo_fields} and repeat names {repeat_names}")
        return "Unknown"

    logging.info(f"Matching combo file {os.path.basename(ref_combo)} using repeat names: {repeat_names}")

    for row in data_rows:
        combo_id = row[0]              # e.g., 'tr001'
        expected_values = row[1:]      # e.g., ['A042', 'B008']

        # Check if every expected value is present in the corresponding repeat's observed hits
        matched_all = True
        for repeat_name, expected in zip(repeat_names, expected_values):
            observed_set = repeat_hits.get(repeat_name, set())

            if expected not in observed_set:
                matched_all = False
                break
            else:
                logging.info(f"Checking combo {combo_id}: repeat '{repeat_name}' needs '{expected}', observed: {observed_set}")

        if matched_all:
            logging.info(f"Matched combo {combo_id} with: {dict(zip(repeat_names, expected_values))}")
            return combo_id

    return "Unknown"

# ------------------------- Main Typing Logic ------------------------- #

def run_repeat_typing(fasta_path: str, repeat_names: List[str], combo_names: List[str], ref_seq: List[str], ref_meta: List[str]) -> Dict[str, str]:
    """
    Perform repeat typing on an assembly and optionally match combinations.

    Args:
        fasta_path (str): Path to input FASTA.
        repeat_names (List[str]): Repeat type identifiers (e.g., TR6, TR10).
        combo_names (List[str]): Combination type identifiers (e.g., TRST).
        db_dir (str): Path to database with repeat fragments and type definitions.

    Returns:
        Dict[str, str]: Typing result dictionary including repeat types and combo types.
    """
    
    contigs = list(SeqIO.parse(fasta_path, "fasta"))
    logging.info(f"Found {len(contigs)} contigs in {fasta_path}")

    results = {}
    repeat_hits = {}

    nl="\n" #fstring hack

    # Step 1: Detect patterns for each repeat name
    for name in repeat_names:
        ref_fasta = [fasta for fasta in ref_seq if name in fasta]
        if len(ref_fasta) > 1:
            logging.warning(f"Multiple TR fasta files with similar names found, will only use first of:\n - {nl.join(ref_fasta)}")
        elif len(ref_fasta) < 1:
            raise(f"Aborting. No reference fasta file found for Type Repeats: {name}")
        ref_fasta = ref_fasta[0]
        """
        ==> TR10_repeat_sequences.fa <==
        >N001
        AAATTAATTATTATATTTCTTT
        >N002
        AAATTAATTTTCTATATTTCTT

        ==> TR6_repeat_sequences.fa <==
        >R001
        CTTGCATACCACTAATAGTGC
        >R002
        CTTGCATATCACTAATAGTAC
        """

        ref_txt = [meta for meta in ref_meta if name in meta]
        if len(ref_txt) > 1:
            logging.warning(f"Multiple TR metadata files with similar names found, will only use first of:\n - {nl.join(ref_txt)}")
        elif len(ref_txt) < 1:
            raise(f"Aborting. No reference metadata file found for Type Repeats: {name}")
        ref_txt = ref_txt[0]
        """
        ==> TR10.txt <==
        B001,   N003-N002-N004-N005-N004-N006-N004-N006-N004-N007-N004-N008-N009-N010-N004-N005-N004-N011
        B002,   N003-N002-N004-N008-N009-N005-N004-N014-N015-N011

        ==> TR6.txt <==
        A001,   R044-R011-R051-R007-R023-R024-R052-R053-R004-R019-R018-R019-R018-R022-R035
        A002,   R044-R011-R051-R007-R020-R007-R015-R008-R009-R022-R035        
        """

        if not (os.path.exists(ref_fasta) and os.path.exists(ref_txt)):
            logging.warning(f"Missing files for repeat type {name}. Skipping.")
            results[name] = "Unknown"
        else:

            fragments = parse_fasta(ref_fasta)
            type_map = parse_repeat_types(ref_txt, fragments)
            hits = set()
            logging.debug(f"Loaded {len(fragments)} sequence fragments from {ref_fasta}")

            for record in contigs:
                hits.update(find_repeat_hits(str(record.seq),str(record.id),type_map))

            results[name] = ";".join(sorted(hits)) if hits else "Unknown"
            repeat_hits[name] = hits

            #Log the exact hits for this repeat name
            logging.info(f"Repeat hits for {name}: {sorted(hits) if hits else 'None'}")

    # Step 2: Check for any combination tables (≥3-column combo.txt files)
    for combo_name in combo_names:
        ref_combo = [meta for meta in ref_meta if combo_name in meta]
        if len(ref_combo) > 1:
            logging.warning(f"Multiple TR metadata files with similar names found, will only use first of:\n - {nl.join(ref_combo)}")
        elif len(ref_combo) < 1:
            raise(f"Aborting. No reference metadata file found for Type Repeats: {combo_name}")
        ref_combo = ref_combo[0]
        """
        ==> TRST.txt <==
        ST_COMB ST_A    ST_B
        tr001   A001    B001
        tr002   A002    B002
        tr003   A003    B003
        """

        if not os.path.exists(ref_combo):
            logging.warning(f"Missing combination file for {combo_name}. Skipping.")
            results[combo_name] = "Unknown"
        else:

            with open(ref_combo) as f:
                header = f.readline().strip().split("\t")

            if len(header) < 3: #need at least three columns, the combined type, and then two or more repeats to be combined.
                logging.warning(f"Combination file {ref_combo} does not have ≥3 columns. Skipping.")
                results[combo_name] = "Unknown"
                continue

            combo_id = match_combination_from_table(ref_combo, repeat_hits, repeat_names)
            results[combo_name] = combo_id

    return results
# ---------------------------- Main CLI ---------------------------- #

def main(args):
    setup_logging(log_file = args.log_file)

    try:
        results = run_repeat_typing(
            fasta_path=args.fasta,
            repeat_names=args.repeats,
            combo_names=[args.combos] if isinstance(args.combos, str) else args.combos,
            ref_seq=args.ref_seq,
            ref_meta=args.ref_meta
        )
    except Exception as e:
        logging.error(f"Typing failed for {args.sample_id}: {e}")
        raise

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
    parser.add_argument("--combos", nargs="+", required=True, help="List of combination identifiers")
    parser.add_argument("--ref_seq", nargs="+", help="Path to TR fasta files.")
    parser.add_argument("--ref_meta", nargs="+", help="Path to TR and TRST metadata files.")
    parser.add_argument("--output", required=True, help="Output file path")
    parser.add_argument("--suffix", choices=["tsv", "csv"], default="tsv", help="Output file format")
    parser.add_argument("--log_file", required = True, help = "File for logging")
    args = parser.parse_args()
    main(args)