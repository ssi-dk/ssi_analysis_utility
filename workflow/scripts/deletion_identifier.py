#!/usr/bin/env python3
from __future__ import annotations

import argparse
from typing import Dict, List, Tuple
import pandas as pd
import pysam

# -------------------------- Helpers --------------------------

def reverse_complement(seq: str) -> str:
    comp = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
    return "".join(comp.get(b.upper(), b) for b in reversed(seq))

def load_bed6(bed6_path: str) -> Dict[str, Dict[str, object]]:
    """Load BED6 into {gene: {contig, start, end, length, strand}}."""
    bedfile_dict: Dict[str, Dict[str, object]] = {}

    try:
        df = pd.read_csv(
            bed6_path,
            sep="\t",
            header=0,
            names=["contig", "start", "end", "gene", "score", "strand"],
            dtype={"contig": "string", "start": int, "end": int, "gene": "string", "score": "float", "strand": "string"},
            )
        for idx, row in df.iterrows():
            bedfile_dict[str(row["gene"])]= {
                "contig": str(row["contig"]),
                "start": int(row["start"]),
                "end": int(row["end"]),
                "length": (int(row["end"]) - int(row["start"])),
                "strand": str(row["strand"])
            }
        print(bedfile_dict)
    except Exception as e:
        print(f"Failed to load BED6 file: {e}")
        raise
    
    # 'tcdC': {'contig': 'AM180355.1', 'start': 804309, 'end': 805008, 'length': 699, 'strand': '-'}, 
    return bedfile_dict

def process_res_file(res_file_path: str) -> pd.DataFrame:
    # i need these excepts if the .res file if later changed to temp() or somehow altered during the snakemake pipeline
    try:
        res_df = pd.read_csv(res_file_path, sep="\t")
    except FileNotFoundError:
        raise FileNotFoundError(f"File not found: {res_file_path}")
    except pd.errors.EmptyDataError:
        raise ValueError(f"File is empty or not properly formatted: {res_file_path}")

    required_columns = {"#Template", "Template_length","Template_Coverage", "Template_Identity", "Depth"}
    if not required_columns.issubset(res_df.columns):
        raise ValueError(f"Missing expected columns in {res_file_path}. Found: {', '.join(res_df.columns)}")

    return res_df

def find_contig_for_gene(res_df: pd.DataFrame, gene: str) -> str:
    hits = res_df[res_df["#Template"].str.contains(gene, case=False, na=False)]
    if hits.empty:
        raise ValueError(f"No contig with gene '{gene}' in .res")
    return str(hits.iloc[0]["#Template"])

def read_meta(meta_path: str, organism: str) -> pd.DataFrame:
    """Read deletion meta TSV and filter by species."""
    df = pd.read_csv(meta_path, sep="\t")
    
    required_col = {"species", "gene", "del_start", "del_end", "gt_IMF", "gt_IDV", "gt_DP", "consensus_N"}
    
    missing_col = required_col - set(df.columns)
    if missing_col:
        raise ValueError(f"provided snp information file to (--meta) is missing required columns: {', '.join(sorted(missing_col))}")
    
    # normalize dtypes
    df["species"] = df["species"].astype(str)
    df = df[df["species"] == organism].copy()
    if df.empty:
        return df

    df["gene"] = df["gene"].astype(str)
    for pos in ("del_start", "del_end"):
        df[pos] = df[pos].astype(int)
    for threshold in ("gt_IMF", "gt_IDV", "gt_DP", "consensus_N"):
        df[threshold] = pd.to_numeric(df[threshold])

    # add length and a key like del330_347_17
    df["del_len"] = (df["del_end"] - df["del_start"] + 1).astype(int)
    print(f"the length {df['del_len']}")
    df["del_key"] = df.apply(
        lambda r: f"del{r.del_start}_{r.del_end}_{r.del_len}", axis=1
    )
    print(f"key is {df['del_key']}")
    return df

def convert_reverse_strand_regions(
    regions: Dict[int, Tuple[int, int]], gene_length: int
) -> Dict[int, Tuple[int, int]]:
    """Convert reverse-strand gene-relative coords to forward contig coords."""
    converted: Dict[int, Tuple[int, int]] = {}
    for length_key, (start, end) in regions.items():
        converted_start = gene_length - (end - 1)
        converted_end = gene_length - (start - 1)
        converted[length_key] = (min(converted_start, converted_end), max(converted_start, converted_end))
    return converted

def is_in_any_region(pos: int, spans: List[Tuple[int, int]]) -> bool:
    return any(a <= pos <= b for a, b in spans)

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

def get_consensus_threshold_for_length(meta_subset: pd.DataFrame, del_len: int) -> float:
    """
    Return consensus_N threshold for a given deletion length.
    - If exactly one row for this length: use it.
    - If multiple rows: return the MAX consensus_N (conservative).
    - If none: return 100.0 (effectively disables auto-upgrade).
    """
    rows = meta_subset[meta_subset["del_len"] == del_len]
    if rows.empty:
        return 100.0
    return float(rows["consensus_N"].max())

# -------------------- Core BCF scanning ----------------------
def get_thresholds_for_deletion(
    meta_subset: pd.DataFrame,
    pos: int,
    del_len: int,
) -> Tuple[float, float, float, float, str]:
    """
    Returns (IMF_thr, IDV_thr, DP_thr, consensus_N, deletion_key).

    Tries exact window first: start=pos, end=pos+del_len-1.
    Falls back to per-length iff there is exactly one row for that length.
    Raises KeyError if ambiguous or not found.
    """
    start = pos
    end = pos + del_len - 1

    row = meta_subset[(meta_subset["del_start"] == start) & (meta_subset["del_end"] == end)]
    if not row.empty:
        r = row.iloc[0]
        key = f"del{int(r.del_start)}_{int(r.del_end)}_{del_len}"
        return float(r.gt_IMF), float(r.gt_IDV), float(r.gt_DP), float(r.consensus_N), key

    by_len = meta_subset[meta_subset["del_len"] == del_len]
    if len(by_len) == 1:
        r = by_len.iloc[0]
        key = f"del{int(r.del_start)}_{int(r.del_end)}_{del_len}"
        return float(r.gt_IMF), float(r.gt_IDV), float(r.gt_DP), float(r.consensus_N), key

    raise KeyError(f"No unambiguous thresholds for pos={pos}, len={del_len}")

def check_deletions_in_region(
    bcf_path: str,
    contig: str,
    target_regions: Dict[int, Tuple[int, int]],
    meta_subset: pd.DataFrame,
    region_buffer: int = 5,
    length_tolerance: int = 1,
    min_overlap_fraction: float = 0.3,
    use_indels_thresholds: bool = False,
) -> Tuple[str, str, int]:
    """
    Scan BCF and try to match deletions to expected regions (by length + position).
    Returns: (label(s) ';'-joined or '-', detail(s) ';'-joined or '-', best_matched_length)
    """
    try:
        bcf = pysam.VariantFile(bcf_path)
    except Exception as e:
        print(f"[ERROR] Cannot open BCF {bcf_path}: {e}")
        return "-", "-", 0

    # Step 0: Compute buffered scan regions for quick filtering of indels
    print("[STEP 0] Computed desired regions to examine")
    padded_regions = {
        expected_len: (
            max(1, region_start - region_buffer),
            region_end + region_buffer
        )
        for expected_len, (region_start, region_end) in target_regions.items()
    }
    region_spans = list(padded_regions.values())
    print(f"{padded_regions}")
    print(f"[INFO] Will scan BCF only within {len(region_spans)} buffered target regions")

    # Step 1: Filter on variant entries.
    print("\n[STEP 1] Filtering deletions based on thresholds for QUAL == 0, thresholds and the spanning region")

    # Collect candidate records inside any region span, with length + thresholds
    candidates: List[Tuple[int, str, str, int]] = []  # (pos, ref, alt, del_len)

    for rec in bcf.fetch(contig):
        if not is_in_any_region(rec.pos, region_spans):
            continue

        ref = rec.ref
        alt = rec.alts[0] if rec.alts else None
        pos = rec.pos
        qual = float(rec.qual) if rec.qual is not None else 0.0

        # deletions only
        if alt is None or len(ref) <= len(alt):
            continue

        del_len = len(ref) - len(alt)

        # only consider expected lengths
        if del_len not in target_regions:
            continue

        if qual > 0.0:
            if not use_indels_thresholds:
                # High-quality deletion passes without IMF/IDV/DP checks
                candidates.append((rec.pos, ref, alt, del_len))
                print(f"variant at pos:{rec.pos} ref:{ref} alt:{alt} passed QUAL>{qual}")
            # When use_indels_thresholds=True we SKIP high-QUAL records (parity with original)
            continue
        
        # QUAL == 0 → evaluate IMF/IDV/DP if using thresholds
        if use_indels_thresholds:

            try:
                IMF_thr, IDV_thr, DP_thr, _, deletion_key = get_thresholds_for_deletion(meta_subset, pos, del_len)
            except KeyError:
                # ambiguous/missing thresholds → skip (same spirit as old code)
                continue

            IMF = float(rec.info.get("IMF", 0.0))
            IDV = float(rec.info.get("IDV", 0.0))
            DP  = float(rec.info.get("DP", 0.0))
            print(f"variant: pos:{pos} ref:{ref} alt:{alt} qual:{qual} IMF:{IMF} IDV:{IDV} DP:{DP}")

            if IMF >= IMF_thr and IDV >= IDV_thr and DP >= DP_thr:
                candidates.append((pos, ref, alt, del_len))
                print(f"  [PASS] {deletion_key} IMF/IDV/DP >= {IMF_thr}/{IDV_thr}/{DP_thr}")
            else:
                print(f"  [FAIL] {deletion_key} IMF/IDV/DP <  {IMF_thr}/{IDV_thr}/{DP_thr}")
        else:
            # Not using thresholds & QUAL==0 → accept (matches original else-branch)
            candidates.append((pos, ref, alt, del_len))
            print(f"  [NO FILTER] accepted QUAL==0 {deletion_key} at pos:{pos}")

    if not candidates:
        print("[INFO] No deletions passed filters.")
        return "-", "-", 0

    # Step 2: sort/merge adjacent-or-overlapping candidates
    print("[STEP 2] Merging overlapping deletions")
    candidates.sort()
    merged: List[Tuple[int, int, int, str]] = []  # (start, end, total_len, details_str)
    current = None
    for pos, ref, alt, dlen in candidates:
        start = pos
        end = pos + dlen - 1
        detail_str = f"{pos}_{ref}_{alt}_{dlen}"
        if current is None:
            current = [start, end, dlen, detail_str]
        else:
            if start <= current[1] + 1:
                # merge
                current[1] = max(current[1], end)
                current[2] += dlen
                current[3] += f"+{detail_str}"
            else:
                merged.append(tuple(current))
                current = [start, end, dlen, detail_str]
    if current is not None:
        merged.append(tuple(current))

    # Step 3: match (exact → partial → ambiguous)
    labels: List[str] = []
    details: List[str] = []
    best_len = 0

    for match_start, match_end, match_len, det in merged:
        actual_len = match_end - match_start + 1
        matched = False

        # (1) exact: same length and fully within padded window
        print("[STEP 3] Checking for exact matches")
        for L, (window_start, window_end) in padded_regions.items():
            if actual_len == L and window_start <= match_start and match_end <= window_end:
                labels.append(str(L))
                details.append(det)
                best_len = max(best_len, L)
                matched = True
                break
        if matched:
            continue

        # (2) partial: near length and overlap fraction satisfied
        print(f"[STEP 4] Checking partial overlaps with overlap percentage {min_overlap_fraction*100}%")
        for expected_len, (window_start, window_end) in padded_regions.items():
            if abs(actual_len - expected_len) > length_tolerance:
                if actual_len > 0:
                    print(f"Skipped {actual_len}bp vs {expected_len}bp: outside length tolerance")
                continue
        
            overlap_start = max(match_start, window_start)
            overlap_end = min(match_end, window_end)
            overlap = max(0, overlap_end - overlap_start + 1)
            target_len = window_end - window_start + 1
            overlap_fraction = overlap / target_len

            if target_len > 0 and overlap_fraction >= min_overlap_fraction:
                labels.append(f"partial_{expected_len}")
                details.append(det)
                best_len = max(best_len, expected_len)
                matched = True
                break
        if matched:
            continue

        print("[STEP 5] Checking for rescued match based on early alignment")
        for expected_len, (window_start, window_end) in padded_regions.items():

            if actual_len == expected_len and abs(match_start - window_start) <= region_buffer:
                print(f"[MATCH] Rescued early match for {expected_len}bp deletion at {match_start}-{window_start}")
                labels.append(f"rescued_{expected_len}")
                details.append(det)
                best_len = max(best_len, expected_len)
                matched = True
                break
        if matched:
            continue

        # (3) ambiguous: any overlap with any window
        print("[STEP 6] Checking for ambiguous overlaps")
        overlapped = [L for L, (window_start, window_end) in padded_regions.items()
                      if not (match_end < window_start or match_start > window_end)]
        if overlapped:
            # synthesize compact detail using merged ref sequence
            parts = []
            for part in det.split("+"):
                p = part.split("_")
                if len(p) >= 3:
                    s = int(p[0])
                    ref = p[1]
                    parts.append((s, s + len(ref) - 1, ref))
            ms, mseq, mlen = merge_overlapping_deletions_sliced(parts)
            # choose window whose expected length is closest
            L = min(overlapped, key=lambda x: abs(x - actual_len))
            labels.append(f"ambiguous_{L}")
            details.append(f"{ms}_{mseq}_{mseq[:1] if mseq else ''}_{max(0, mlen-1)}")
            best_len = max(best_len, L)
        else:
            # no overlap with expectations; still carry raw details
            details.append(det)

    if not labels:
        return "-", ";".join(d[3] for d in merged), best_len

    return ";".join(labels), ";".join(details), best_len

def extract_consensus_window(
    fasta_path: str,
    contig: str,
    del_len_key: int,
    converted_regions: Dict[int, Tuple[int, int]],
    consensus_threshold_N: float
) -> str:
    """Extract region by deletion length and annotate %N. Returns '-' if below threshold."""
    try:
        if del_len_key not in converted_regions:
            return "-"
        s, e = converted_regions[del_len_key]
        fa = pysam.FastaFile(fasta_path)
        seq = fa.fetch(contig, s, e)
        n_pct = (seq.upper().count("N") / max(1, len(seq))) * 100.0
        if n_pct >= consensus_threshold_N:
            return f"{s}-{e}_{seq}_{n_pct:.2f}"
        return "-"
    except Exception as e:
        print(f"[WARN] consensus extraction failed: {e}")
        return "-"

def extract_best_ambiguous_consensus(
    fasta_path: str,
    contig: str,
    converted_regions: Dict[int, Tuple[int, int]],
    meta_subset: pd.DataFrame
) -> str:
    """
    Scan all expected regions and return ONLY a threshold-passing region:
      - return "start-end_SEQUENCE_%N" for the region above its length-specific threshold
        with the highest %N (ties naturally broken by iteration order)
      - return "-" if NO region meets its threshold
    """
    try:
        fasta = pysam.FastaFile(fasta_path)
    except Exception as e:
        print(f"[WARN] cannot open FASTA {fasta_path}: {e}")
        return "-"

    best = ("-", -1.0, None, "-")  # (region_str, n_pct, L, seq)

    for L, (s, e) in converted_regions.items():
        try:
            seq = fasta.fetch(contig, s, e)
        except Exception:
            continue

        """
        # scoring the region with the highest %N (only checking it passes the threshold), ignoring length.
        n_pct = (seq.upper().count("N") / max(1, len(seq))) * 100.0
        consN = get_consensus_threshold_for_length(meta_subset, L)

        # Only consider regions that meet their threshold
        if n_pct >= consN and n_pct > best[1]:
            best = (f"{s}-{e}", n_pct, L, seq)
        """

        """
        scored candidates by score = %N * region_len and picked the highest score among those above threshold considering the length.

        39 bp at ~71.79% (score ≈ 2800) beats 18 bp at 100% (score = 1800).
        """
        region_len = max(1, e - s + 1)
        n_pct = (seq.upper().count("N") / region_len) * 100.0
        consN = get_consensus_threshold_for_length(meta_subset, L)

        if n_pct >= consN:
            score = n_pct * region_len  # length-aware like the old code
            if score > best[1]:
                best = (f"{s}-{e}", n_pct, L, seq, score)

    if best[2] is not None:
        return f"{best[0]}_{best[3]}_{best[1]:.2f}"
    return "-"

def evaluate_ambiguous_consensus(label: str, consensus_str: str, meta_subset: pd.DataFrame) -> str:
    """
    If label is 'ambiguous_L' and the consensus %N ≥ threshold_N for length L,
    upgrade to 'likely_L'.
    """
    if not label.startswith("ambiguous_"):
        return label
    try:
        n_pct = float(consensus_str.strip().split("_")[-1])
    except Exception:
        return label
    try:
        L = int(label.split("_")[1])
        consN = get_consensus_threshold_for_length(meta_subset, L)
        if n_pct >= consN:
            return f"likely_{L}"
    except Exception:
        pass
    return label

# --------------------------- Identifying deletions ---------------------------

def run(
    organism: str,
    res_path: str,
    call_bcf: str,
    indels_bcf: str,
    fasta_path: str,
    bed_path: str,
    meta_path: str,
    output_path: str,
    region_buffer: int,
    length_tolerance: int,
    partial_overlap: float,
) -> None:

    res_df = process_res_file(res_path)
    bed = load_bed6(bed_path)
    meta = read_meta(meta_path, organism)

    # If no rows for this organism → emit dashes
    if meta.empty:
        pd.DataFrame([{
            "organism": organism,
            # emit no gene-specific columns; downstream joiners can fill if needed
        }]).to_csv(output_path, sep="\t", index=False)
        return

    # Group meta by gene
    genes = list(dict.fromkeys(meta["gene"]))  # stable unique order
    row: Dict[str, str] = {"organism": organism}

    for gene in genes:
        meta_g = meta[meta["gene"] == gene].copy()

        # find contig and strand/length from BED
        try:
            contig = find_contig_for_gene(res_df, gene)
        except Exception as e:
            print(f"[WARN] {gene}: contig not found in .res → '-' ({e})")
            row[f"{gene}del"] = "-"
            row[f"{gene}_deletion_details"] = "-"
            row[f"{gene}_deletion_consensus"] = "-"
            continue

        if gene not in bed:
            print(f"[WARN] {gene}: not in BED6 → '-'")
            row[f"{gene}del"] = "-"
            row[f"{gene}_deletion_details"] = "-"
            row[f"{gene}_deletion_consensus"] = "-"
            continue

        strand = str(bed[gene]["strand"])
        gene_len = int(bed[gene]["length"])

        # Build expected regions dicts keyed by length
        # gene-relative (meta file is gene coordinates)
        orig_regions: Dict[int, Tuple[int, int]] = {}
        for _, r in meta_g.iterrows():
            L = int(r["del_len"])
            s, e = int(r["del_start"]), int(r["del_end"])
            # Keep max span per length if duplicates
            if L in orig_regions:
                s0, e0 = orig_regions[L]
                orig_regions[L] = (min(s0, s), max(e0, e))
            else:
                orig_regions[L] = (s, e)

        # If minus strand, convert to forward (contig) coords
        if strand == "-":
            converted = convert_reverse_strand_regions(orig_regions, gene_len)
        elif strand == "+":
            # plus-strand contigs are 1..gene_len; meta is already gene-relative → shift to 1..gene_len
            converted = orig_regions.copy()
        else:
            print(f"[WARN] {gene}: invalid strand '{strand}' → '-'")
            row[f"{gene}del"] = "-"
            row[f"{gene}_deletion_details"] = "-"
            row[f"{gene}_deletion_consensus"] = "-"
            continue

        # 1) Main BCF
        label, detail, best_len = check_deletions_in_region(
            bcf_path=call_bcf,
            contig=contig,
            target_regions=converted,
            meta_subset=meta_g,
            region_buffer=region_buffer,
            length_tolerance=length_tolerance,
            min_overlap_fraction=partial_overlap,
            use_indels_thresholds=True
        )

        # 2) If ambiguous or none, try indels BCF
        if label.startswith("ambiguous") or label in ("-", "", None):
            label2, detail2, best_len2 = check_deletions_in_region(
                bcf_path=indels_bcf,
                contig=contig,
                target_regions=converted,
                meta_subset=meta_g,
                region_buffer=region_buffer,
                length_tolerance=length_tolerance,
                min_overlap_fraction=partial_overlap,
                use_indels_thresholds=True
            )
            if label2 not in ("-", "", None):
                label, detail, best_len = label2, detail2, best_len2

        # 3) Consensus extraction
        #    - if we have a confident (or at least length-bearing) label, use that length
        #    - else scan all and pick best ambiguous region
        if isinstance(best_len, int) and best_len in converted:
            # consensus threshold from meta
            try:
                consN = get_consensus_threshold_for_length(meta_g, best_len)
            except KeyError:
                consN = 100.0
            consensus = extract_consensus_window(
                fasta_path=fasta_path,
                contig=contig,
                del_len_key=best_len,
                converted_regions=converted,
                consensus_threshold_N=consN,
            )
            # optional upgrade ambiguous→likely if consensus supports
            if label.startswith("ambiguous") and consensus != "-":
                label = evaluate_ambiguous_consensus(label, consensus, meta_g)
        else:
            consensus = extract_best_ambiguous_consensus(
                fasta_path=fasta_path,
                contig=contig,
                converted_regions=converted,
                meta_subset=meta_g,
            )
            # If we had no call but consensus suggests a specific region, emit a "likely_<len>"
            if label in ("-", "", None) and consensus not in ("-", "", None):
                try:
                    reg = consensus.split("_")[0]
                    s, e = map(int, reg.split("-"))
                    L = e - s + 1
                    label = f"likely_{L}"
                except Exception:
                    pass

        row[f"{gene}del"] = label if label else "-"
        row[f"{gene}_deletion_details"] = filter_deletion_details(detail) if detail else "-"
        row[f"{gene}_deletion_consensus"] = consensus if consensus else "-"

    pd.DataFrame([row]).to_csv(output_path, sep="\t", index=False)

# ----------------------------- CLI -----------------------------

def main():
    ap = argparse.ArgumentParser(description="Identify configured deletions using a meta TSV (no YAML configs).")
    ap.add_argument("--organism", required=True)
    ap.add_argument("--res", required=True, help="KMA .res file (for contig lookup)")
    ap.add_argument("--call", required=True, help="Main genotype BCF")
    ap.add_argument("--indels", required=True, help="Indels BCF")
    ap.add_argument("--fsa", required=True, help="Consensus FASTA (kmer consensus)")
    ap.add_argument("--bed", required=True, help="BED6 with gene coordinates/strand")
    ap.add_argument("--metafile", required=True, help="TSV with expected deletion windows and thresholds")
    ap.add_argument("-o", "--output", required=True, help="Output TSV path")
    # matching/scan parameters
    ap.add_argument("--partial_overlap", type=float, default=0.3, help="Overlap fraction for partial matches (default 0.3)")
    ap.add_argument("--partial_match_length_tolerance", type=int, default=1, help="±bp tolerance for partial length matches (default 1)")
    ap.add_argument("--deletion_region_buffer", type=int, default=5, help="bp to pad each expected region (default 5)")
    args = ap.parse_args()

    run(
        organism=args.organism,
        res_path=args.res,
        call_bcf=args.call,
        indels_bcf=args.indels,
        fasta_path=args.fsa,
        bed_path=args.bed,
        meta_path=args.metafile,
        output_path=args.output,
        region_buffer=args.deletion_region_buffer,
        length_tolerance=args.partial_match_length_tolerance,
        partial_overlap=args.partial_overlap,
    )

if __name__ == "__main__":
    main()
