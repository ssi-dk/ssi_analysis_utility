#!/usr/bin/env python3
import pandas as pd
import argparse
import re
import os
from typing import List, Dict
import pysam
import logging
from logging_utils import setup_logging
from thresholds import get_kma_thresholds_for_species, get_threshold, get_deletion_threshold, deletion_gt_thresholds, deletion_consensus_thresholds
from Bio import SeqIO
from Bio.Seq import Seq

# ========================= KMA FILE HANDLING =========================
def process_res_file(res_file_path: str, thresholds: Dict[str, List[int]]) -> pd.DataFrame:
    try:
        res_df = pd.read_csv(res_file_path, sep="\t")
    except FileNotFoundError:
        raise FileNotFoundError(f"File not found: {res_file_path}")
    except pd.errors.EmptyDataError:
        raise ValueError(f"File is empty or not properly formatted: {res_file_path}")

    # query_identity can be penalized more if the consensus
    required_columns = {"#Template", "Template_length","Template_Coverage", "Template_Identity", "Depth"}
    if not required_columns.issubset(res_df.columns):
        raise ValueError(f"Missing expected columns in {res_file_path}. Found: {', '.join(res_df.columns)}")

    res_df["threshold"] = res_df["#Template"].apply(lambda x: get_threshold(x, thresholds))
    res_df_filtered = res_df[
        (res_df["Template_Coverage"] >= res_df["threshold"].apply(lambda x: x[0])) &
        (res_df["Template_Identity"] >= res_df["threshold"].apply(lambda x: x[1])) &
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
        return output_data, pd.DataFrame()

    known_genes = {"tcdA", "tcdB", "tcdC", "cdtAB", "cdtA", "cdtB"}
    other_genes = set()

    for template in filtered_df["#Template"]:
        # Attempt to split like E. coli format
        gene_list = template.split("_")
        #print(f"filtered template {filtered_df['#Template']}")
        #print(f"gene list {gene_list}")

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
                length = row["Template_length"]
                coverage = row["Template_Coverage"] # breadth of coverage
                adjusted_coverage = min(coverage, 100.0)
                depth = row["Depth"] # depth of ceverage
                identity = row["Template_Identity"] 
                verbose_parts.append(f"{gene}_{length}_{adjusted_coverage:.2f}_{identity:.2f}_{depth:.2f}")
        output_data["verbose"] = ";".join(verbose_parts)

    logging.info(f"Successfully processed sample: {sample_name}")
    return output_data, filtered_df

# ========================= TOXIN VARIATION FILE HANDLING BCF =========================

def load_toxin_coordinates(bed6_path: str) -> Dict[str, Dict[str, str | int]]:
    """
    Load coordinates of toxin genes from a BED6 file into a dictionary.

    Args:
        bed6_path (str): Path to the BED6 file.

    Returns:
        Dict[str, Dict[str, str | int]]: Dictionary mapping gene -> {contig, start, end, length, strand}
    """
    #print("INSIDE TOXIN")
    coords = {}
    try:
        bed_df = pd.read_csv(bed6_path, sep="\t", header=None,
                             names=["contig", "start", "end", "gene", "score", "strand"])

        for _, row in bed_df.iterrows():
            #print(f"row in bed_df {row}")
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

    #{'tcdB': {'contig': 'AM180355.1', 'start': 787392, 'end': 794493, 'length': 7101, 'strand': '+'}, 
    # 'tcdA': {'contig': 'AM180355.1', 'start': 795842, 'end': 803975, 'length': 8133, 'strand': '+'}, 
    # 'tcdC': {'contig': 'AM180355.1', 'start': 804309, 'end': 805008, 'length': 699, 'strand': '-'}, 
    # 'cdtA': {'contig': 'AF271719.1', 'start': 0, 'end': 1392, 'length': 1392, 'strand': '+'}, 
    # 'cdtB': {'contig': 'AF271719.1', 'start': 1444, 'end': 4075, 'length': 2631, 'strand': '+'}}

    return coords

def convert_reverse_strand_regions(regions: dict[int, tuple[int, int]], gene_length: int) -> dict[int, tuple[int, int]]:
    """
    Convert reverse strand gene-relative coordinates to forward strand contig-relative coordinates.

    Args:
        regions (dict): Dictionary of deletion IDs to (start, end) coordinates (gene-relative, reverse strand).
        gene_length (int): Length of the gene.

    Returns:
        dict: Converted dictionary with forward strand contig-relative coordinates.
    """
    converted = {}
    for key, (start, end) in regions.items():
        # Reverse strand conversion
        converted_start = gene_length - (end - 1)
        converted_end = gene_length - (start - 1)
        converted[key] = (min(converted_start, converted_end), max(converted_start, converted_end))
        #print(f"original start:{start}-{end} \t converted:{converted_start}-{converted_end}")
    return converted

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
    length = end - start
    
    if strand == "+":
        genomic_pos = 1 + (pos_in_gene - 1)
    elif strand == "-":
        genomic_pos = length - (pos_in_gene - 1)
    else:
        raise ValueError(f"Invalid strand '{strand}' for gene {gene_name}.")

    return contig, genomic_pos

def check_tcdC117_variant(bcf_path: str, contig: str, pos: int, range: int,expected_ref="T", expected_alt="A") -> str:
    """
    Check for a specific variant at a position in the BCF, considering strand orientation.
    
    Args:
        bcf_path (str): Path to the BCF file.
        contig (str): Contig name.
        pos (int): 1-based genomic position.
        expected_ref (str): Expected REF base (e.g., 'T')
        expected_alt (str): Expected ALT base (e.g., 'A')
    
    Returns:
        str: '1' if variant found (T>A on minus strand), '0' otherwise.
    """

    try:
        bcf = pysam.VariantFile(bcf_path)
        #print("Opened BCF successfully")
        for rec in bcf.fetch(contig, pos-range, pos+range):
            #print(f"Record at {rec.contig}:{rec.pos}")
            #print(f"REF: {rec.ref}")
            #print(f"ALT: {rec.alts}")

            ref = rec.ref
            alt = rec.alts[0] if rec.alts else None

            #print(f"Record v2 at {rec.contig}:{rec.pos}")
            #print(f"REF v2: {rec.ref}")
            #print(f"ALT v2: {rec.alts}")

             
            if rec.pos == pos:
                # Case 1: variation fitting the A>T
                if ref == expected_ref and alt == expected_alt:
                    #print("SNP PRESENT")
                    return "A>T", ""
                # Case 2: different variation
                else:
                    #print("DIFFERENT VARIANT")
                    return "other", f"{rec.pos}_{ref}_{alt}"
            # Case 3: Deletion spanning the position
            elif rec.pos < pos:
                deletion_end = rec.pos + len(ref) - 1
                if deletion_end >= pos and any(len(a) < len(ref) for a in rec.alts if a is not None):
                    #print(f"Deletion spans position {pos}: {rec.pos}-{deletion_end}")
                    return "del", f"{rec.pos}_{ref}_{alt}"

        #print("No variant found at or around position.")
        return "wt", ""

    except Exception as e:
        #print(f"Error reading BCF at {contig}:{pos}: {e}")
        return "-", ""

def filter_deletion_details(detail_str: str) -> str:
    """
    Filters a deletion detail string by removing entries with 0-length deletions.

    Args:
        detail_str (str): A string of the format "pos1_REF_ALT_len+pos2_REF_ALT_len+..."

    Returns:
        str: Filtered string with only actual deletions (len > 0), or '-' if none remain.
    """
    if not detail_str:
        return "-"

    parts = detail_str.split("+")
    valid_parts = []

    for part in parts:
        fields = part.split("_")
        if len(fields) == 4:
            try:
                if int(fields[3]) > 0:
                    valid_parts.append(part)
            except ValueError:
                continue  # skip malformed entries
        else:
            continue  # skip incomplete entries

    return "+".join(valid_parts) if valid_parts else "-"

def merge_overlapping_deletions_sliced(deletions: List[tuple[int, int, str]]) -> tuple[int, str, int]:
    """
    Merge overlapping deletions by slicing sequences based on genomic overlap.

    Args:
        deletions (List[Tuple[int, int, str]]): Each entry is (start, end, ref_sequence).

    Returns:
        Tuple[int, str, int]: Merged start position, merged sequence, merged length.
    """
    if not deletions:
        return -1, "", 0

    deletions.sort(key=lambda x: x[0])
    merged_start, merged_end, merged_seq = deletions[0]
    for start, end, seq in deletions[1:]:
        overlap = merged_end - start + 1
        if overlap > 0:
            merged_seq += seq[overlap:]
        else:
            merged_seq += seq
        merged_end = max(merged_end, end)

    merged_len = merged_end - merged_start + 1
    return merged_start, merged_seq, merged_len

def check_deletions_in_region(
    bcf_path: str,
    contig: str,
    gene_name: str,
    orig_regions: dict[int, tuple[int, int]],
    target_regions: dict[int, tuple[int, int]],
    region_buffer: int = 5,
    length_tolerance: int = 1,
    min_overlap_fraction: float = 0.4,
    use_indels_thresholds: bool = False 
) -> tuple[str, str, int]:
    """
    Scan for deletions in a BCF, merge overlapping ones, and match merged spans to expected regions.
    """
    print(f"\n[INFO] Scanning deletions in {bcf_path} for {gene_name}")
    try:
        bcf = pysam.VariantFile(bcf_path)
        raw_deletions = []
        filtered_records = []

        # Step 0: Filter on variant entries.
        print("[STEP 0] Filtering deletions based on thresholds for QUAL == 0 and Depth")
        for rec in bcf.fetch(contig):
            ref = rec.ref
            alt = rec.alts[0] if rec.alts else None
            pos = rec.pos
            qual = rec.qual

            IMF = rec.info.get("IMF", 0)
            IDV = rec.info.get("IDV", 0)
            DP = rec.info.get("DP", 0)
                       
            if float(qual) > 0:
                if not use_indels_thresholds:
                    filtered_records.append((pos, ref, alt))
                continue
            else:
                print(f"variant: \t pos:{pos} : ref:{ref} \t alt:{alt} \t qual:{qual} \t IMF:{IMF} \t IDV:{IDV} \t DP:{DP}")
                # Skip non-deletions
                if not alt or len(ref) <= len(alt):
                    continue
                
                # Now apply deletion-specific filtering
                del_len = len(ref) - len(alt)
                matched = False

                for expected_len, (region_start, region_end) in target_regions.items():
                    if del_len != expected_len:
                        continue

                    del_len = len(ref) - len(alt)
                    del_start = pos
                    del_end = pos + del_len - 1

                    orig_del_reg = orig_regions[del_len]
                    
                    deletion_key = f"del{orig_del_reg[0]}_{orig_del_reg[1]}_{del_len}"
                    print(f"deletion key: {deletion_key}")

                    if use_indels_thresholds:
                        try:
                            thresholds = get_deletion_threshold(deletion_key, deletion_gt_thresholds)
                            print(f"  → Loaded thresholds for {deletion_key}: IMF ≥ {thresholds[0]}, IDV ≥ {thresholds[1]}, DP ≥ {thresholds[2]}")
                        except ValueError:
                            continue  # no thresholds for this deletion

                        if IMF >= thresholds[0] and IDV >= thresholds[1] and DP >= thresholds[2]:
                            filtered_records.append((pos, ref, alt))
                            print(f"  [PASS] variant {pos}:{ref}-{alt} passed {deletion_key} threshold check")
                        else:
                            print(f"  [FAIL] variant {pos}:{ref}-{alt} failed {deletion_key} failed threshold check")                
                        matched = True
                        break  # use only first matching region
                    else:
                        filtered_records.append((pos, ref, alt))
                        print(f"  [NO FILTER] accepted variant {pos}:{ref}-{alt} without thresholds")
                        matched = True
                        break  # use only first matching region

                if not matched:
                    print("  [SKIP] No matching deletion region found for filtering.")

        # Step 1: Extract all deletions.
        print("[STEP 1] Extracting deletions")
        raw_deletions = []
        for pos, ref, alt in filtered_records:
            del_len = len(ref) - len(alt)
            del_start = pos
            del_end = pos + del_len - 1
            detail_str = f"{pos}_{ref}_{alt}_{del_len}"
            raw_deletions.append((del_start, del_end, del_len, detail_str))
            if del_len >0:
                print(f"  → Deletion: {del_len}bp at {del_start}-{del_end}")

        if not raw_deletions:
            print("[INFO] No deletions found.")
            return "-", "-",0

        # Step 2: Merge overlapping deletions.
        print("[STEP 2] Merging overlapping deletions")
        merged_deletions = []
        raw_deletions.sort()
        current = list(raw_deletions[0])
        for d in raw_deletions[1:]:
            if d[0] <= current[1] + 1:
                current[1] = max(current[1], d[1])
                current[2] += d[2]
                current[3] += f"+{d[3]}"
            else:
                merged_deletions.append(tuple(current))
                current = list(d)
        merged_deletions.append(tuple(current))

        # Step 3: Exact match.
        print("[STEP 3] Checking for exact matches")
        padded_regions = {
            expected_len: (
                max(1, region_start - region_buffer),
                region_end + region_buffer
            )
            for expected_len, (region_start, region_end) in target_regions.items()
        }

        for del_start, del_end, _, detail_str in merged_deletions:
            actual_len = del_end - del_start + 1

            if actual_len > 0:
                print(f"  → Merged deletion: {actual_len}bp at {del_start}-{del_end}")
    
            for expected_len, (region_start, region_end) in padded_regions.items():
                if actual_len > 0:
                    print(f"    Comparing to expected {expected_len}bp region (± 5 nt): {region_start}-{region_end}")
                
                if actual_len == expected_len and region_start <= del_start and del_end <= region_end:
                    print(f"[MATCH] Exact match: {expected_len}bp deletion")
                    return str(expected_len), detail_str, int(expected_len)

        # Step 4: Fallback partial match.
        print("[STEP 4] Checking partial overlaps")
        for del_start, del_end, _, detail_str in merged_deletions:
            actual_len = del_end - del_start + 1
    
            if actual_len > 0:
                print(f"  → Checking merged deletion {actual_len}bp at {del_start}-{del_end}")
    
            for expected_len, (region_start, region_end) in padded_regions.items():
                if abs(actual_len - expected_len) > length_tolerance:
                    if actual_len > 0:
                        print(f"    Skipped {actual_len}bp vs {expected_len}bp: outside length tolerance")
                    continue
    
                overlap_start = max(del_start, region_start)
                overlap_end = min(del_end, region_end)
                overlap = max(0, overlap_end - overlap_start + 1)
                target_len = region_end - region_start + 1
                overlap_fraction = overlap / target_len
                print(f"    → Fallback for {expected_len}bp: overlap = {overlap}/{target_len} = {overlap_fraction:.2f}")
                if overlap_fraction >= min_overlap_fraction:
                    print(f"[MATCH] Fallback partial match for {expected_len}bp (≥ {min_overlap_fraction})")
                    return f"partial_{expected_len}", detail_str, int(expected_len)

        # Step 5: Recovered early match (start-proximal rescue).
        print("[STEP 5] Checking for rescued match based on early alignment")
        for expected_len, (region_start, region_end) in padded_regions.items():
            for del_start, del_end, _, detail_str in merged_deletions:
                actual_len = del_end - del_start + 1
                dist = abs(del_start - region_start)

                if actual_len > 0:
                    print(f"  → Checking {actual_len}bp deletion at {del_start}-{del_end} vs region {region_start}-{region_end} with distance between start coordinate {dist} nt")

                if actual_len == expected_len and dist <= 5:
                    print(f"[MATCH] Rescued early match for {expected_len}bp deletion at {del_start}-{del_end}")
                    return f"rescued_{expected_len}", detail_str, int(expected_len)

        # Step 6: Ambiguous overlap (last resort).
        print("[STEP 6] Checking for ambiguous overlaps")
        #print(f"  → merged deletions {merged_deletions}")
        for del_start, del_end, _, detail_str in merged_deletions:
            if del_len <= 0:
                continue
            
            for expected_len, (region_start, region_end) in padded_regions.items():
                if del_end >= region_start and del_start <= region_end:
                    print(f"  → [AMBIGUOUS] Deletion at {min(del_start,del_end)}-{max(del_start,del_end)} overlaps expected {expected_len}bp region")
                    # Extract actual ref seqs from detail_str like: "321_REF1_ALT1_LEN+330_REF2_ALT2_LEN"
                    deletion_details = detail_str.split("+")
                    parsed = []
                    for d in deletion_details:
                        parts = d.split("_")
                        #print(f"PARTS {parts}")
                        if len(parts) >= 4:
                            s, ref, alt, del_len = parts
                            parsed.append((int(s), int(s)+len(ref)-1, ref))

                    merged_start, merged_seq, merged_len = merge_overlapping_deletions_sliced(parsed)
                    merged_info = f"{merged_start}_{merged_seq}_{merged_seq[0]}_{merged_len-1}"
                    return f"ambiguous_{expected_len}", merged_info, int(expected_len)
        
        # Step 7: Nothing matched.
        print("[STEP 7] No deletions matched any target region.")
        all_details = ";".join([d[3] for d in merged_deletions])
        # for the merged deletions think about what you want to insert here
        return "-", all_details, int(merged_deletions[2][2])

    except Exception as e:
        print(f"[ERROR] While scanning deletions: {e}")
        return "-", "-", 0

def verify_bcf(bcf_path: str, indels_bcf_path: str, contig: str, 
               orig_regions: dict[int, tuple[int, int]], target_regions: dict[int, tuple[int, int]], 
               region_buffer: int = 5, length_tolerance: int = 1) -> tuple[str, str,int]:
    """
    Verify deletions in two BCF files: the main BCF and the indels BCF.
    
    Args:
        bcf_path (str): Path to the main BCF file.
        indels_bcf_path (str): Path to the indels BCF file.
        contig (str): Contig name.
        target_regions (dict): Target deletion regions (already converted to forward strand).
        region_buffer (int): Buffer size (in nt) applied symmetrically to region boundaries.
        length_tolerance (int): Allowed deviation for matching deletion lengths.

    Returns:
        tuple[str, str]: Deletion status and details (e.g., "18;36" for matched deletions and "pos_REF_ALT_len" for deletion details).
    """

    #print(f"Verifying deletions in BCF files: {bcf_path} and {indels_bcf_path}")
    #print("\n[INFO] First attempt: scanning primary BCF")

    del_status, del_details, expected_del_len = check_deletions_in_region(
        bcf_path=bcf_path,
        contig=contig,
        gene_name="tcdC",
        orig_regions=orig_regions,
        target_regions=target_regions,
        region_buffer=region_buffer,
        length_tolerance=length_tolerance,
        use_indels_thresholds=False
    )

    # If result is ambiguous or no confident match, try the indels file
    if del_status.startswith("ambiguous") or del_status in ("-", "", None):
        print(f"[INFO] No confident deletion in main BCF {bcf_path}. Now checking indels BCF: {indels_bcf_path}")
        del_status_indels, del_details_indels, expected_del_len = check_deletions_in_region(
            bcf_path=indels_bcf_path,
            contig=contig,
            gene_name="tcdC",
            orig_regions=orig_regions,
            target_regions=target_regions,
            region_buffer=region_buffer,
            length_tolerance=length_tolerance,
            use_indels_thresholds=True
        )
        #print(f"orig region {orig_regions} target region {target_regions}")
        if del_status_indels not in ("-", "", None):
            if del_status_indels.startswith("ambiguous"):
                print(f"[INFO] Found ambiguous deletion in indels BCF: {del_status_indels}")
            else:
                print(f"[INFO] Found confident deletion in indels BCF with length {del_status_indels} and expected length {expected_del_len}")
            return del_status_indels, del_details_indels, expected_del_len
        else:
            print("[INFO] Still no confident match even in indels BCF.")

    return del_status, del_details, expected_del_len

def find_matching_contig(filtered_df: pd.DataFrame, gene_name: str) -> str:
    """
    Searches the filtered .res file dataframe for a contig name containing the given gene.
    
    Args:
        filtered_df (pd.DataFrame): Filtered .res file.
        gene_name (str): Gene of interest (e.g. 'tcdC').

    Returns:
        str: Matching contig name (e.g. 'AM180355.1_tcdC_804309_805008')
    """
    matches = filtered_df[filtered_df["#Template"].str.contains(gene_name, case=False)]
    if not matches.empty:
        return matches.iloc[0]["#Template"]
    raise ValueError(f"No contig found for gene: {gene_name}")

# ========================= CONSENSUS SEQUENCE HANDLING - KMA alignment =========================
def extract_consensus_indel_seq(
    fasta_path: str,
    contig: str,
    del_expected_len: int,
    converted_regions: dict[int, tuple[int, int]],
) -> str:
    """
    Extracts the gene region from FASTA using deletion length key, and annotates with %N.

    Args:
        fasta_path (str): Path to the FASTA file.
        contig (str): Contig name containing the gene (e.g., AM180355.1_tcdC_...).
        del_expected_len (int): Expected deletion length key (e.g., 54).
        converted_regions (dict): Dictionary mapping deletion length → (start, end).

    Returns:
        str: Annotated sequence string like ATGCNNNN_33.33 or "-" on failure.
    """
    try:
        if del_expected_len not in converted_regions:
            raise ValueError(f"Deletion length {del_expected_len} not found in regions")

        print("[STEP 8] Checking the consensus sequence for deleted regions")

        start, end = converted_regions[del_expected_len]
        fasta = pysam.FastaFile(fasta_path)
        seq = fasta.fetch(contig, start, end)

        print(f"  → Consensus sequence region: {start}-{end}")

        n_count = seq.upper().count("N")
        n_percent = (n_count / len(seq)) * 100 if len(seq) > 0 else 0
        return f"{seq}_{n_percent:.2f}"

    except Exception as e:
        logging.warning(f"Failed to extract/annotate from {contig}:{del_expected_len} in {fasta_path}: {e}")
        return "-"

def evaluate_ambiguous_consensus(deletion_label: str, consensus_str: str, thresholds: dict) -> str:
    """
    Upgrade 'ambiguous_XX' to 'likely_XX' if %N ≥ threshold in deletion_consensus_thresholds.

    Args:
        deletion_label (str): e.g. 'ambiguous_39'
        consensus_str (str): e.g. 'ATGCNNNN_88.88'
        thresholds (dict): mapping of deletion keys to [N% threshold]

    Returns:
        str: updated deletion label ('likely_XX' or original ambiguous)
    """
    if not deletion_label.startswith("ambiguous") or "_" not in consensus_str:
        return deletion_label

    try:
        n_percent = float(consensus_str.strip().split("_")[-1])
    except Exception:
        return deletion_label

    try:
        print("[STEP 9] Checking support for deletions within consensus sequence using thresholds")

        del_len = deletion_label.split("_")[1]
        for key, value in thresholds.items():
            if key.endswith(f"_{del_len}"):
                required_n = value[0]
                print(f"  → Checking consensus sequence with current label {deletion_label}")
                if n_percent >= required_n:  # upgraded if enough Ns
                    deletion_label = deletion_label.replace("ambiguous", "likely")
                    print(f"  → Consensus sequence with {n_percent}% N's passed thresholds of {required_n}% N's - supporting deletion, changing label: {deletion_label}")
                    return deletion_label
                else:
                    print(f"  → Consensus sequence with {n_percent}% N's failed thresholds of {required_n}% N's - keeping label: {deletion_label}")
                break
    except Exception:
        pass

    return deletion_label

# ========================= TRST FILE HANDLING - assembly =========================

def parse_fasta(filepath):
    return {record.id: str(record.seq) for record in SeqIO.parse(filepath, "fasta")}

def parse_types(filepath, fragments):
    types = {}
    with open(filepath) as f:
        for line in f:
            if ",\t" in line:
                key, pattern = line.strip().split(",\t")
                try:
                    types[key] = ''.join(fragments[p] for p in pattern.split("-"))
                except KeyError:
                    continue
    return types

def find_matches(sequence, patterns):
    hits = []
    for name, pattern in patterns.items():
        pattern_seq = Seq(pattern)
        if re.search(str(pattern_seq), sequence, re.IGNORECASE) or re.search(str(pattern_seq.reverse_complement()), sequence, re.IGNORECASE):
            hits.append(name)
    return hits

def load_trst_types(filepath):
    trst_table = []
    with open(filepath) as f:
        for line in f:
            trst, tr6, tr10 = line.strip().split("\t")
            trst_table.append((trst, tr6, tr10))
    return trst_table

def match_trst(rTR6, rTR10, trst_table):
    for trst, tr6, tr10 in trst_table:
        if tr6 in rTR6 and tr10 in rTR10:
            return trst
    return "Unknown"

def run_trst_typing_on_fasta(fasta_path, db_dir):
    TR6_frags = parse_fasta(os.path.join(db_dir, "TR6_repeat_sequences.fa"))
    TR10_frags = parse_fasta(os.path.join(db_dir, "TR10_repeat_sequences.fa"))

    TR6_types = parse_types(os.path.join(db_dir, "TR6_repeat_types.txt"), TR6_frags)
    TR10_types = parse_types(os.path.join(db_dir, "TR10_repeat_types.txt"), TR10_frags)

    trst_table = load_trst_types(os.path.join(db_dir, "TRST_repeat_types.txt"))

    all_TR6_hits = set()
    all_TR10_hits = set()
    contig_count = 0

    for record in SeqIO.parse(fasta_path, "fasta"):
        contig_count += 1
        seq = str(record.seq)

        rTR6 = find_matches(seq, TR6_types)
        rTR10 = find_matches(seq, TR10_types)

        all_TR6_hits.update(rTR6)
        all_TR10_hits.update(rTR10)

    # Try to find a matching TRST
    for trst, tr6, tr10 in trst_table:
        if tr6 in all_TR6_hits and tr10 in all_TR10_hits:
            return {"TRST": trst, "TR6": tr6, "TR10": tr10}

    # No match found → preserve raw TR6/TR10 hits if available
    TR6_val = ";".join(sorted(all_TR6_hits)) if all_TR6_hits else "Unknown"
    TR10_val = ";".join(sorted(all_TR10_hits)) if all_TR10_hits else "Unknown"

    return {"TRST": "Unknown", "TR6": TR6_val, "TR10": TR10_val}

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

    # All expected columns
    input_columns += ["nanopore_read_file", "assembly_file", "organism", "variant", "notes"]
    summary_columns = ["tcdA", "tcdB", "tcdC", "cdtAB", "Other", "verbose"]
    variation_columns = ["tcdC117", "tcdCdel", "deletion_details", "deletion_consensus"]
    TRST_columns = ["TRST", "TR6", "TR10"]
    desired_columns = input_columns + summary_columns + variation_columns + TRST_columns

    # Load coordinates from BED
    coord_dict = load_toxin_coordinates("resources/Clostridioides_difficile_db/Toxin/Cdiff_Toxin.bed6")

    # Reverse-strand gene-relative deletion regions
    """
    18 + 39 https://pmc.ncbi.nlm.nih.gov/articles/PMC1828959/ 
    18 https://journals.asm.org/doi/pdf/10.1128/jcm.02340-15
    
    ALL
    https://journals.asm.org/doi/10.1128/jcm.00866-08
    """
    tcdC_deletion_regions = {
        18: (330, 347),
        36: (301, 336),
        39: (341, 379),
        54: (313, 366)
    }

    result_rows = []

    for idx, row in df.iterrows():
        row_dict = row.to_dict()
        sample = row_dict["sample_name"]
        print(f"\n================= Processing sample: {sample} =================")

        # -------------------------- KMA --------------------------
        res_path = f"examples/Results/{sample}/Cdiff_KMA_Toxin/{sample}.res"

        try:
            thresholds = get_kma_thresholds_for_species(row_dict["organism"])
        except ValueError as e:
            logging.error(f"Sample '{sample}' skipped: {e}")
            continue

        result_dict, filtered_df = summarize_single_sample(sample, res_path, verbose_flag=args.verbose, thresholds=thresholds)
        row_dict.update(result_dict)

        # ----------------------- BCFtools variations -----------------------
        bcf_path = f"examples/Results/{sample}/GenotypeCalls/{sample}.Cdiff_KMA_Toxin.calls.bcf"
        bcf_indel_path = f"examples/Results/{sample}/GenotypeCalls/{sample}.Cdiff_KMA_Toxin.indels.bcf"

        if filtered_df.empty:
            row_dict["tcdC117"] = "-"
            row_dict["tcdCdel"] = "-"
            row_dict["deletion_consensus"] = "-"
        else:
            try:
                contig = find_matching_contig(filtered_df, "tcdC")
                _, genomic_pos = gene_pos_to_genomic("tcdC", 117, coord_dict)
                #print(f"Position 117 in tcdC maps to contig pos: {genomic_pos}")

                variant_status, extra_verbose = check_tcdC117_variant(
                    bcf_path, contig, genomic_pos, range=20, expected_ref="T", expected_alt="A"
                )
                row_dict["tcdC117"] = variant_status

                if extra_verbose:
                    row_dict["verbose"] = (
                        extra_verbose if row_dict["verbose"] in ("-", "") else f"{row_dict['verbose']};{extra_verbose}"
                    )
            except Exception as e:
                logging.warning(f"Failed tcdC117 SNP check for {sample}: {e}")
                row_dict["tcdC117"] = "-"

            try:
                gene_length = coord_dict["tcdC"]["length"]
                converted_regions = convert_reverse_strand_regions(tcdC_deletion_regions, gene_length)

                del_status, del_details, del_expected_len = verify_bcf(
                    bcf_path=bcf_path,
                    indels_bcf_path=bcf_indel_path,
                    contig=contig,
                    orig_regions=tcdC_deletion_regions,
                    target_regions=converted_regions,
                    region_buffer=5,
                    length_tolerance=1
                )
                row_dict["tcdCdel"] = del_status
                
                row_dict["deletion_details"] = filter_deletion_details(del_details)
                #print(f"deletion region keys {tcdC_deletion_regions.keys()} and values {tcdC_deletion_regions.values()}")

                if del_expected_len in tcdC_deletion_regions:
                    fasta_path = f"examples/Results/{sample}/Cdiff_KMA_Toxin/{sample}.fsa"
                    tcdC_seq_with_n =extract_consensus_indel_seq(fasta_path,contig,del_expected_len,converted_regions)
                    row_dict["deletion_consensus"] = tcdC_seq_with_n

                    # Evaluate if ambiguous should be upgraded to likely
                    if row_dict["tcdCdel"].startswith("ambiguous") and row_dict["deletion_consensus"] != "-":
                        row_dict["tcdCdel"] = evaluate_ambiguous_consensus(
                            row_dict["tcdCdel"],
                            row_dict["deletion_consensus"],
                            deletion_consensus_thresholds
                        )
                    
                else:
                    row_dict["deletion_consensus"] = "-"
            except Exception as e:
                logging.warning(f"Failed tcdC deletion check or FASTA extraction for {sample}: {e}")
                row_dict["tcdCdel"] = "-"
                row_dict["deletion_consensus"] = "-"

        # ----------------------- TRST Typing -----------------------
        assembly_path = f"examples/Results/{sample}/skesa/{sample}.contigs.fasta"
        trst_db_path = os.path.join("resources", "Clostridioides_difficile_db", "TRST")

        try:
            trst_result = run_trst_typing_on_fasta(assembly_path, trst_db_path)
            row_dict["TRST"] = trst_result["TRST"]
            row_dict["TR6"] = trst_result["TR6"]
            row_dict["TR10"] = trst_result["TR10"]
        except Exception as e:
            logging.warning(f"TRST typing failed for {sample}: {e}")
            row_dict["TRST"] = "-"
            row_dict["TR6"] = "-"
            row_dict["TR10"] = "-"

        for col in desired_columns:
            row_dict.setdefault(col, "-")
        result_rows.append(row_dict)

        if args.outputfile is None:
            output_path = os.path.join("examples", "Results", sample, "kma", f"{sample}_kma.{args.suffix}")
            os.makedirs(os.path.dirname(output_path), exist_ok=True)
            pd.DataFrame([row_dict])[desired_columns].to_csv(output_path, sep="\t" if args.suffix == "tsv" else ",", index=False)
            logging.info(f"Wrote per-sample summary: {output_path}")

    if args.outputfile:
        for row_dict in result_rows:
            for col in desired_columns:
                row_dict.setdefault(col, "-")
        pd.DataFrame(result_rows)[desired_columns].to_csv(args.outputfile, sep="\t" if args.suffix == "tsv" else ",", index=False)
        logging.info(f"Wrote merged summary: {args.outputfile}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Append C. difficile KMA .res summary info to a samplesheet.")
    parser.add_argument("-ss", "--samplesheet", required=True, help="Path to the input samplesheet .tsv file")
    parser.add_argument("-o", "--outputfile", default=None, help="Output filename. Default: per-sample output only")
    parser.add_argument("--suffix", choices=["tsv", "csv"], default="tsv", help="Output format. Default: tsv")
    parser.add_argument("--verbose", type=int, choices=[0, 1], default=1, help="Include verbose output. Default: 1")
    parser.add_argument("--split", type=int, choices=[0, 1], default=0, help="Split Illumina_read_files into Read1/Read2 (default: 0)")
    args = parser.parse_args()
    main(args)
