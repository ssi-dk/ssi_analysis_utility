#!/usr/bin/env python3
import pandas as pd
import argparse
import os
import sys
from typing import List, Dict
import pysam
import logging
from logging_utils import setup_logging
import yaml
#sys.path.insert(0, os.path.abspath("../scripts"))

# ========================= Helper files =========================

def get_deletion_threshold(deletion_key: str, thresholds: Dict[str, List[float]]) -> List[float]:
    """
    Returns the [IMF, IDV, DP] threshold list for a given deletion ID.
    
    Args:
        deletion_key (str): Deletion ID, e.g. 'del330_347_18'

    Returns:
        List[float]: List containing [IMF, IDV, DP] thresholds

    Raises:
        ValueError if the key is not found
    """
    for key in thresholds:
        if key in deletion_key:
            return thresholds[key]
    raise ValueError(f"No deletion thresholds found for: {deletion_key}")

def load_variant_detection_config(organism: str, config_dir="workflow/configs_species") -> dict:
    """
    Load variant detection thresholds and region info from a species-specific YAML config file.

    Args:
        organism (str): Full organism name, e.g. "Clostridioides difficile"
        config_dir (str): Directory containing species YAML files.

    Returns:
        dict: Dictionary with keys:
            - genomic_coord (str)
            - snp_info (dict)
            - deletion_regions (dict[int] -> (start, end))
            - variant_gt_thresholds (dict[str] -> list[float|int])
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

    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Species config not found: {config_path}")
    
    with open(config_path, "r") as f:
        config = yaml.safe_load(f)

    try:
        variant_block = config["analyses_to_run"]["Variant_detection"]
        bed_path = variant_block["genomic_coord"]
        snp_info = variant_block.get("snp_info", {})
        raw_deletions = variant_block.get("deletion_regions", {})
        raw_thresholds = variant_block.get("variant_gt_thresholds", {})
        consensus = variant_block.get("deletion_consensus_thresholds",{})
    except KeyError as e:
        raise ValueError(f"Missing required Variant_detection fields in config: {e}")

    deletion_regions = {}
    for key, val in raw_deletions.items():
        # del_330_347_18: [330, 347, 70, "tcdC"] → 18: (330, 347)
        length = int(key.split("_")[-1])
        start, end = val[0], val[1]
        deletion_regions[length] = (start, end)

    return {
        "genomic_coord": bed_path,
        "snp_info": snp_info,
        "deletion_regions": deletion_regions,
        "variant_gt_thresholds": raw_thresholds,
        "deletion_consensus_thresholds": consensus
    }

# ========================= KMA FILE HANDLING =========================
def process_res_file(res_file_path: str) -> pd.DataFrame:
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

    return res_df  # no filtering

# ========================= TOXIN VARIATION FILE HANDLING BCF =========================

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
        logging.info(f"Converted regions - original start:{start}-{end} \t converted:{converted_start}-{converted_end}")
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
        for rec in bcf.fetch(contig, pos-range, pos+range):

            ref = rec.ref
            alt = rec.alts[0] if rec.alts else None
             
            if rec.pos == pos:
                # Case 1: variation fitting the A>T
                if ref == expected_ref and alt == expected_alt:
                    logging.info("Case 1: A>T SNP present at 117")
                    return "A>T", ""
                # Case 2: different variation
                else:
                    logging.info("Case 2: Different SNP present at 117")
                    return "other", f"{rec.pos}_{ref}_{alt}"
            # Case 3: Deletion spanning the position
            elif rec.pos < pos:
                deletion_end = rec.pos + len(ref) - 1
                if deletion_end >= pos and any(len(a) < len(ref) for a in rec.alts if a is not None):
                    logging.info("Case 3: A deletion spans the position at 117")
                    return "del", f"{rec.pos}_{ref}_{alt}"

        logging.info("No variant found at or around position.")
        return "wt", ""

    except Exception as e:
        logging.error(f"Error reading BCF at {contig}:{pos}: {e}")
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
    deletion_gt_thresholds: dict[str, list[float,int, int]],
    region_buffer: int = 5,
    length_tolerance: int = 1,
    min_overlap_fraction: float = 0.4,
    use_indels_thresholds: bool = False,
) -> tuple[str, str, int]:
    """
    Scan for deletions in a BCF, merge overlapping ones, and match merged spans to expected regions.
    """
    logging.info(f"\n[INFO] Scanning deletions in {bcf_path} for {gene_name}")
    try:
        bcf = pysam.VariantFile(bcf_path)
        raw_deletions = []
        filtered_records = []

        # Step 0: Filter on variant entries.
        logging.info("[STEP 0] Filtering deletions based on thresholds for QUAL == 0 and Depth")
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
                logging.info(f"variant: \t pos:{pos} : ref:{ref} \t alt:{alt} \t qual:{qual} \t IMF:{IMF} \t IDV:{IDV} \t DP:{DP}")
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

                    if use_indels_thresholds:
                        try:
                            thresholds = get_deletion_threshold(deletion_key, deletion_gt_thresholds)
                            print(thresholds)
                            print("\n\nAFSFSAFFSAFSAFSFSAFSA\n\n")
                            logging.info(f"  → Loaded thresholds for {deletion_key}: IMF ≥ {thresholds[0]}, IDV ≥ {thresholds[1]}, DP ≥ {thresholds[2]}")
                        except ValueError:
                            continue  # no thresholds for this deletion

                        if IMF >= thresholds[0] and IDV >= thresholds[1] and DP >= thresholds[2]:
                            filtered_records.append((pos, ref, alt))
                            logging.info(f"  [PASS] variant {pos}:{ref}-{alt} passed {deletion_key} threshold check")
                        else:
                            logging.info(f"  [FAIL] variant {pos}:{ref}-{alt} failed {deletion_key} failed threshold check")                
                        matched = True
                        break  # use only first matching region
                    else:
                        filtered_records.append((pos, ref, alt))
                        logging.info(f"  [NO FILTER] accepted variant {pos}:{ref}-{alt} without thresholds")
                        matched = True
                        break  # use only first matching region

                if not matched:
                    logging.info("  [SKIP] No matching deletion region found for filtering.")

        # Step 1: Extract all deletions.
        logging.info("[STEP 1] Extracting deletions")
        raw_deletions = []
        for pos, ref, alt in filtered_records:
            del_len = len(ref) - len(alt)
            del_start = pos
            del_end = pos + del_len - 1
            detail_str = f"{pos}_{ref}_{alt}_{del_len}"
            raw_deletions.append((del_start, del_end, del_len, detail_str))
            if del_len >0:
                logging.info(f"  → Deletion: {del_len}bp at {del_start}-{del_end}")

        if not raw_deletions:
            logging.info("[INFO] No deletions found.")
            return "-", "-",0

        # Step 2: Merge overlapping deletions.
        logging.info("[STEP 2] Merging overlapping deletions")
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
        logging.info("[STEP 3] Checking for exact matches")
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
                logging.info(f"  → Merged deletion: {actual_len}bp at {del_start}-{del_end}")
    
            for expected_len, (region_start, region_end) in padded_regions.items():
                if actual_len > 0:
                    logging.info(f"    Comparing to expected {expected_len}bp region (± 5 nt): {region_start}-{region_end}")
                
                if actual_len == expected_len and region_start <= del_start and del_end <= region_end:
                    logging.info(f"[MATCH] Exact match: {expected_len}bp deletion")
                    return str(expected_len), detail_str, int(expected_len)

        # Step 4: Fallback partial match.
        logging.info("[STEP 4] Checking partial overlaps")
        for del_start, del_end, _, detail_str in merged_deletions:
            actual_len = del_end - del_start + 1
    
            if actual_len > 0:
                logging.info(f"  → Checking merged deletion {actual_len}bp at {del_start}-{del_end}")
    
            for expected_len, (region_start, region_end) in padded_regions.items():
                if abs(actual_len - expected_len) > length_tolerance:
                    if actual_len > 0:
                        logging.info(f"    Skipped {actual_len}bp vs {expected_len}bp: outside length tolerance")
                    continue
    
                overlap_start = max(del_start, region_start)
                overlap_end = min(del_end, region_end)
                overlap = max(0, overlap_end - overlap_start + 1)
                target_len = region_end - region_start + 1
                overlap_fraction = overlap / target_len
                logging.info(f"    → Fallback for {expected_len}bp: overlap = {overlap}/{target_len} = {overlap_fraction:.2f}")
                if overlap_fraction >= min_overlap_fraction:
                    logging.info(f"[MATCH] Fallback partial match for {expected_len}bp (≥ {min_overlap_fraction})")
                    return f"partial_{expected_len}", detail_str, int(expected_len)

        # Step 5: Recovered early match (start-proximal rescue).
        logging.info("[STEP 5] Checking for rescued match based on early alignment")
        for expected_len, (region_start, region_end) in padded_regions.items():
            for del_start, del_end, _, detail_str in merged_deletions:
                actual_len = del_end - del_start + 1
                dist = abs(del_start - region_start)

                if actual_len > 0:
                    logging.info(f"  → Checking {actual_len}bp deletion at {del_start}-{del_end} vs region {region_start}-{region_end} with distance between start coordinate {dist} nt")

                if actual_len == expected_len and dist <= 5:
                    logging.info(f"[MATCH] Rescued early match for {expected_len}bp deletion at {del_start}-{del_end}")
                    return f"rescued_{expected_len}", detail_str, int(expected_len)

        # Step 6: Ambiguous overlap (last resort).
        logging.info("[STEP 6] Checking for ambiguous overlaps")
        for del_start, del_end, _, detail_str in merged_deletions:
            if del_len <= 0:
                continue
            
            for expected_len, (region_start, region_end) in padded_regions.items():
                if del_end >= region_start and del_start <= region_end:
                    logging.info(f"  → [AMBIGUOUS] Deletion at {min(del_start,del_end)}-{max(del_start,del_end)} overlaps expected {expected_len}bp region")
                    # Extract actual ref seqs from detail_str like: "321_REF1_ALT1_LEN+330_REF2_ALT2_LEN"
                    deletion_details = detail_str.split("+")
                    parsed = []
                    for d in deletion_details:
                        parts = d.split("_")

                        if len(parts) >= 4:
                            s, ref, alt, del_len = parts
                            parsed.append((int(s), int(s)+len(ref)-1, ref))

                    merged_start, merged_seq, merged_len = merge_overlapping_deletions_sliced(parsed)
                    merged_info = f"{merged_start}_{merged_seq}_{merged_seq[0]}_{merged_len-1}"
                    return f"ambiguous_{expected_len}", merged_info, int(expected_len)
        
        # Step 7: Nothing matched.
        logging.info("[STEP 7] No deletions matched any target region.")
        all_details = ";".join([d[3] for d in merged_deletions])
        # for the merged deletions think about what you want to insert here
        if merged_deletions:
            return "-", all_details, int(merged_deletions[0][2])
        else:
            return "-", "-", 0

    except Exception as e:
        logging.error(f"[ERROR] While scanning deletions: {e}")
        return "-", "-", 0

def verify_bcf(bcf_path: str, indels_bcf_path: str, contig: str, 
               orig_regions: dict[int, tuple[int, int]], target_regions: dict[int, tuple[int, int]],deletion_thresholds: dict[str, list[float,int, int]],
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

    del_status, del_details, expected_del_len = check_deletions_in_region(
        bcf_path=bcf_path,
        contig=contig,
        gene_name="tcdC",
        orig_regions=orig_regions,
        target_regions=target_regions,
        deletion_gt_thresholds=deletion_thresholds,
        region_buffer=region_buffer,
        length_tolerance=length_tolerance,
        use_indels_thresholds=True
    )

    # If result is ambiguous or no confident match, try the indels file
    if del_status.startswith("ambiguous") or del_status in ("-", "", None):
        logging.info(f"[INFO] No confident deletion in main BCF {bcf_path}. Now checking indels BCF: {indels_bcf_path}")
        del_status_indels, del_details_indels, expected_del_len = check_deletions_in_region(
            bcf_path=indels_bcf_path,
            contig=contig,
            gene_name="tcdC",
            orig_regions=orig_regions,
            target_regions=target_regions,
            deletion_gt_thresholds=deletion_thresholds,
            region_buffer=region_buffer,
            length_tolerance=length_tolerance,
            use_indels_thresholds=True
        )

        if del_status_indels not in ("-", "", None):
            if del_status_indels.startswith("ambiguous"):
                logging.info(f"[INFO] Found ambiguous deletion in indels BCF: {del_status_indels}")
            else:
                logging.info(f"[INFO] Found confident deletion in indels BCF with length {del_status_indels} and expected length {expected_del_len}")
            return del_status_indels, del_details_indels, expected_del_len
        else:
            logging.info("[INFO] Still no confident match even in indels BCF.")

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

        logging.info("[STEP 8] Checking the consensus sequence for deleted regions")

        start, end = converted_regions[del_expected_len]
        fasta = pysam.FastaFile(fasta_path)
        seq = fasta.fetch(contig, start, end)

        logging.info(f"  → Consensus sequence region: {start}-{end}")

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
        logging.info("[STEP 9] Checking support for deletions within consensus sequence using thresholds")

        del_len = deletion_label.split("_")[1]
        for key, value in thresholds.items():
            if key.endswith(f"_{del_len}"):
                required_n = value[0]
                logging.info(f"  → Checking consensus sequence with current label {deletion_label}")
                if n_percent >= required_n:  # upgraded if enough Ns
                    deletion_label = deletion_label.replace("ambiguous", "likely")
                    logging.info(f"  → Consensus sequence with {n_percent}% N's passed thresholds of {required_n}% N's - supporting deletion, changing label: {deletion_label}")
                    return deletion_label
                else:
                    logging.info(f"  → Consensus sequence with {n_percent}% N's failed thresholds of {required_n}% N's - keeping label: {deletion_label}")
                break
    except Exception:
        pass

    return deletion_label

def main(args: argparse.Namespace) -> None:
    setup_logging(args.log_dir, args.sample_id, "variant_handling")

    sample = args.sample_id
    organism = args.organism
    logging.info(f"\n================= Processing sample: {sample} =================")

    row_dict = {
        "sample_id": sample,
        "organism": organism,
    }

    # Output columns
    first_cols = ["sample_id", "organism"]
    variation_columns = ["tcdC117", "tcdCdel", "deletion_details", "deletion_consensus"]
    desired_columns = first_cols + variation_columns

    # Load BED file based on organism
    # Load variant detection config
    variant_config = load_variant_detection_config(organism)

    print(variant_config)
    print(variant_config["variant_gt_thresholds"]["del330_347_18"][0:3])
    print("aSF ")
    #print(deletion_consensus_thresholds)
    #exit(1)
    bed_path = variant_config["genomic_coord"]
    coord_dict = load_toxin_coordinates(bed_path)

    tcdC_deletion_regions = variant_config["deletion_regions"]
    deletion_gt_thresholds2 = variant_config["variant_gt_thresholds"]

    try:
        res_path = f"examples/Results/{sample}/Cdiff_KMA_Toxin/{sample}.res"
        res_df = process_res_file(res_path)

        bcf_path = f"examples/Results/{sample}/GenotypeCalls/{sample}.Cdiff_KMA_Toxin.calls.bcf"
        bcf_indel_path = f"examples/Results/{sample}/GenotypeCalls/{sample}.Cdiff_KMA_Toxin.indels.bcf"

        if res_df.empty:
            row_dict["tcdC117"] = "-"
            row_dict["tcdCdel"] = "-"
            row_dict["deletion_details"] = "-"
            row_dict["deletion_consensus"] = "-"
        else:
            contig = find_matching_contig(res_df, "tcdC")
            _, genomic_pos = gene_pos_to_genomic("tcdC", 117, coord_dict)

            variant_status, extra_verbose = check_tcdC117_variant(
                bcf_path, contig, genomic_pos, range=20, expected_ref="T", expected_alt="A"
            )
            row_dict["tcdC117"] = variant_status
            if variant_status != "wt":
                extra_info = f"tcdC117:{variant_status}" + (f"_{extra_verbose}" if extra_verbose else "")
                row_dict["deletion_details"] = extra_info
            else:
                row_dict["deletion_details"] = "-"

            gene_length = coord_dict["tcdC"]["length"]
            converted_regions = convert_reverse_strand_regions(tcdC_deletion_regions, gene_length)

            del_status, del_details, del_expected_len = verify_bcf(
                bcf_path=bcf_path,
                indels_bcf_path=bcf_indel_path,
                contig=contig,
                orig_regions=tcdC_deletion_regions,
                target_regions=converted_regions,
                deletion_thresholds = deletion_gt_thresholds2,
                region_buffer=5,
                length_tolerance=1
            )
            row_dict["tcdCdel"] = del_status
            row_dict["deletion_details"] = filter_deletion_details(del_details)

            if del_expected_len in tcdC_deletion_regions:
                fasta_path = f"examples/Results/{sample}/Cdiff_KMA_Toxin/{sample}.fsa"
                consensus_str = extract_consensus_indel_seq(fasta_path, contig, del_expected_len, converted_regions)
                row_dict["deletion_consensus"] = consensus_str

                if del_status.startswith("ambiguous") and consensus_str != "-":
                    row_dict["tcdCdel"] = evaluate_ambiguous_consensus(
                        del_status,
                        consensus_str,
                        variant_config["deletion_consensus_thresholds"]
                    )
            else:
                row_dict["deletion_consensus"] = "-"

    except Exception as e:
        logging.warning(f"[ERROR] Sample {sample} processing failed: {e}")
        row_dict["tcdC117"] = "-"
        row_dict["tcdCdel"] = "-"
        row_dict["deletion_details"] = "-"
        row_dict["deletion_consensus"] = "-"

    for col in desired_columns:
        row_dict.setdefault(col, "-")

    pd.DataFrame([row_dict])[desired_columns].to_csv(
        args.outputfile,
        sep="\t" if args.suffix == "tsv" else ",",
        index=False
    )
    logging.info(f"Wrote summary to: {args.outputfile}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Append C. difficile KMA .res summary info to a samplesheet.")
    parser.add_argument("--sample_id", required=True, help="Sample name to process")
    parser.add_argument("--organism", required=True, help="Organism name (used to load species-specific resources)")
    parser.add_argument("-o", "--outputfile", default=None, help="Output filename. Default: per-sample output only")
    parser.add_argument("--suffix", choices=["tsv", "csv"], default="tsv", help="Output format. Default: tsv")
    parser.add_argument("--verbose", type=int, choices=[0, 1], default=1, help="Include verbose output. Default: 1")
    parser.add_argument("--log_dir", default="examples/Log")

    args = parser.parse_args()
    main(args)