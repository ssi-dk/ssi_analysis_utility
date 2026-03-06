#!/usr/bin/env python3

from __future__ import annotations

import argparse
import re
from typing import Dict, List, Optional, Tuple

import pandas as pd
import pysam


# ========================= FASTA helpers =========================


def load_fasta_sequences(fasta_path: str) -> Dict[str, str]:
    """
    Load a (possibly multi-)FASTA file into a dict: {header: sequence}.

    The header is taken as the first whitespace-separated token after '>'.
    Sequences are uppercased.
    """
    seqs: Dict[str, str] = {}
    header: Optional[str] = None
    chunks: List[str] = []

    print(f"[INFO] Loading FASTA sequences from: {fasta_path}")
    with open(fasta_path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                # flush previous
                if header is not None:
                    seqs[header] = "".join(chunks).upper()
                header = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        # flush last
        if header is not None:
            seqs[header] = "".join(chunks).upper()

    print(f"[INFO] Loaded {len(seqs)} sequences from FASTA.")
    return seqs


# ========================= Step 1: load_metafile =========================


def load_deletion_metafile(meta_path: str) -> pd.DataFrame:
    """
    Load the deletion metafile.

    Expected columns:
      species, gene, del_start, del_end, del_type, gt_IMF, gt_IDV, gt_DP, consensus_N
    """
    df = pd.read_csv(meta_path, sep="\t")
    required_cols = {
        "species",
        "gene",
        "del_start",
        "del_end",
        "del_type",
        "gt_IMF",
        "gt_IDV",
        "gt_DP",
        "consensus_N",
    }

    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(
            f"Metafile is missing required columns: {', '.join(sorted(missing))}"
        )

    # Normalize types
    df["species"] = df["species"].astype(str)
    df["gene"] = df["gene"].astype(str)
    df["del_start"] = df["del_start"].astype(int)
    df["del_end"] = df["del_end"].astype(int)
    df["del_type"] = df["del_type"].astype(int)
    df["gt_IMF"] = pd.to_numeric(df["gt_IMF"])
    df["gt_IDV"] = pd.to_numeric(df["gt_IDV"])
    df["gt_DP"] = pd.to_numeric(df["gt_DP"])
    df["consensus_N"] = pd.to_numeric(df["consensus_N"])

    return df


# ========================= Step 2: Check_organism =========================


def check_organism(meta_df: pd.DataFrame, organism: str) -> pd.DataFrame:
    """
    Filter deletion metafile rows to the requested organism.

    Returns a new DataFrame containing only rows where species == organism.
    """
    organism = str(organism).strip()
    filtered_metadf = meta_df[meta_df["species"].str.strip() == organism].copy()
    print(
        f"[INFO] Filtered metafile for organism '{organism}', "
        f"found {len(filtered_metadf)} rows."
    )
    return filtered_metadf


# ========================= Step 3: Match_gene_bcf =========================


def match_gene_bcf(
    bcf: pysam.VariantFile,
    genes: List[str],
) -> Dict[str, str]:
    """
    For each gene, find a contig in the BCF header whose name contains the gene
    as a substring (regex, case-insensitive).

    Returns a dict: {gene -> contig_name}.

    If multiple contigs match a gene, the first match is used.
    Genes with no contig match are omitted from the dictionary.
    """

    contig_names = list(bcf.header.contigs)
    gene_to_contig: Dict[str, str] = {}

    print(
        f"[INFO] Matching {len(genes)} genes to {len(contig_names)} contigs in BCF header."
    )
    for gene in genes:
        pattern = re.compile(re.escape(gene), flags=re.IGNORECASE)
        matched_contig: Optional[str] = None
        for contig in contig_names:
            if pattern.search(contig):
                matched_contig = contig
                print(
                    f"[DEBUG] gene from metafile '{gene}' matched contig "
                    f"'{contig}' in BCF header."
                )
                break

        if matched_contig is None:
            print(f"[WARN] No contig in BCF matches gene '{gene}' (substring regex).")
            continue

        gene_to_contig[gene] = matched_contig

    return gene_to_contig


# ========================= Step 4: Extract deletion variants in region =========================


def extract_deletion_variants_in_region(
    bcf: pysam.VariantFile,
    contig: str,
    start: int,
    end: int,
    region_buffer: int,
) -> List[pysam.VariantRecord]:
    """
    Extract all deletion variant records from `bcf` on `contig` that overlap
    [start - buffer, end + buffer] (1-based, inclusive).

    A deletion is defined as having at least one ALT shorter than REF.
    """

    # Build padded window
    padded_start_1based = max(1, start - region_buffer)
    padded_end_1based = end + region_buffer

    # Convert to 0-based half-open for pysam.fetch
    fetch_start0 = padded_start_1based - 1
    fetch_end0 = padded_end_1based  # pysam end is exclusive

    print(
        f"[INFO] Extracting deletion variants from contig '{contig}' "
        f"in padded region {padded_start_1based}-{padded_end_1based} "
        f"(original region {start}-{end}, buffer {region_buffer})."
    )

    records: List[pysam.VariantRecord] = []
    try:
        for rec in bcf.fetch(contig, fetch_start0, fetch_end0):
            ref = rec.ref
            alts = rec.alts
            print(
                f"[DEBUG] record on {contig} in padded region: "
                f"pos={rec.pos}, search_window={fetch_start0}-{fetch_end0}, "
                f"ref={ref}, alts={alts}"
            )
            if not alts:
                continue

            # keep only deletion records
            if any(alt is not None and len(alt) < len(ref) for alt in alts):
                # compute deletion span for the primary alt (first shorter alt)
                pos = rec.pos  # 1-based
                for alt in alts:
                    if alt is not None and len(alt) < len(ref):
                        raw_del_len = len(ref) - len(alt)
                        del_start = pos
                        del_end = pos + raw_del_len - 1
                        # require overlap with padded region
                        if (
                            del_end >= padded_start_1based
                            and del_start <= padded_end_1based
                        ):
                            print(
                                f"[DEBUG] Found deletion at {contig}:{del_start}-{del_end} "
                                f"(len={raw_del_len}) overlapping {start}-{end} ± {region_buffer}"
                            )
                            records.append(rec)
                        else:
                            print(
                                f"[DEBUG] Deletion at {contig}:{del_start}-{del_end} "
                                f"(len={raw_del_len}) did NOT overlap region {start}-{end} "
                                f"± {region_buffer}"
                            )
                        break  # only consider first deletion alt per record
    except ValueError as e:
        print(
            f"[ERROR] Could not fetch {contig}:{padded_start_1based}-"
            f"{padded_end_1based} from BCF: {e}"
        )

    print(f"[INFO] Total deletion-like records found in region: {len(records)}")
    return records


# ========================= Step 5: INFO helpers =========================


def _get_info_float(rec: pysam.VariantRecord, key: str) -> float:
    """
    Safely extract a float-like INFO value from a VariantRecord.

    Handles single values and 1-element tuples/lists. Returns 0.0 on failure.
    """
    val = rec.info.get(key)
    if isinstance(val, (int, float)):
        return float(val)
    if isinstance(val, (list, tuple)) and val:
        inner = val[0]
        if isinstance(inner, (int, float)):
            return float(inner)
    return 0.0


def _get_info_int(rec: pysam.VariantRecord, key: str) -> int:
    """
    Safely extract an int-like INFO value from a VariantRecord.

    Handles single values and 1-element tuples/lists. Returns 0 on failure.
    """
    val = rec.info.get(key)
    if isinstance(val, int):
        return val
    if isinstance(val, (list, tuple)) and val:
        inner = val[0]
        if isinstance(inner, int):
            return inner
    return 0


# ========================= Deletion span helper =========================


def get_deletion_span(rec: pysam.VariantRecord) -> Optional[Tuple[int, int, int]]:
    """
    From a deletion-like VariantRecord, return (start, end, length) for
    the first ALT that is shorter than REF. Returns None if no deletion ALT.
    """
    ref = rec.ref
    alts = rec.alts or []
    for alt in alts:
        if alt is not None and len(alt) < len(ref):
            del_len = len(ref) - len(alt)
            del_start = rec.pos
            del_end = rec.pos + del_len - 1
            return del_start, del_end, del_len
    return None


# ========================= Call-based canonical assignment (categories 1–3) =========================


def assign_best_canonical_for_call(
    variants: List[pysam.VariantRecord],
    meta_gene: pd.DataFrame,
    min_frac: float = 0.6,
) -> Tuple[
    Optional[int],
    Optional[int],
    Optional[int],
    Optional[int],
    Optional[int],
    Optional[float],
    Optional[float],
    Optional[bool],
]:
    """
    Given 'call' deletion variants for a single gene and that gene's
    metafile subset (with del_start, del_end, del_type, gt_IMF, gt_IDV, gt_DP),
    pick the best canonical deletion type and annotate whether it passes
    IMF/IDV/DP thresholds.

    Scoring rules (canonical matching):

      1) We first filter on the expected-side fraction:
         frac_expected = overlap / expected_len
         and require frac_expected >= min_frac.

      2) For all remaining candidates, we build a score tuple:

           score_tuple = (
               frac_expected,
               frac_observed,
               -len_diff,
               -exp_type,
           )

         Meaning:

           - Prefer higher frac_expected
               → more of the canonical region is deleted.

           - If tie, prefer higher frac_observed
               → more of the observed deletion lies inside the canonical region.

           - If tie, prefer smaller absolute length difference
               (len_diff = |obs_len - exp_len|).

           - If still tie, prefer smaller canonical length (exp_type)
               → e.g. 39 over 54.

    We maintain two "best" candidates:
      - best_pass: passes IMF/IDV/DP thresholds
      - best_fail: fails at least one of IMF/IDV/DP

    Categories (decided in `run`) are:

      Category 1 (high confidence):
        - best_pass exists AND len_diff <= 1

      Category 2 (good confidence):
        - best_pass exists AND len_diff > 1

      Category 3 (lower confidence):
        - no best_pass, but best_fail exists (any len_diff)

    This function itself does not assign category numbers; it just returns:

      (best_del_type,
       best_len_diff,
       best_obs_start,
       best_obs_end,
       best_obs_len,
       best_frac_expected,
       best_frac_observed,
       best_passes_thresholds)
    """

    if not variants:
        print("[INFO] No call variants provided for canonical assignment.")
        return None, None, None, None, None, None, None, None

    # Best candidates among threshold-passing and threshold-failing
    best_pass_score: Optional[Tuple[float, float, int, int]] = None
    best_fail_score: Optional[Tuple[float, float, int, int]] = None

    best_pass_data = {
        "type": None,
        "len_diff": None,
        "obs_start": None,
        "obs_end": None,
        "obs_len": None,
        "frac_expected": None,
        "frac_observed": None,
    }
    best_fail_data = {
        "type": None,
        "len_diff": None,
        "obs_start": None,
        "obs_end": None,
        "obs_len": None,
        "frac_expected": None,
        "frac_observed": None,
    }

    print(
        f"[INFO] Assigning best canonical deletion from {len(variants)} call variants "
        f"and {len(meta_gene)} canonical windows with min_frac={min_frac}."
    )

    for rec in variants:
        span = get_deletion_span(rec)
        if span is None:
            continue
        obs_start, obs_end, obs_len = span

        IMF = _get_info_float(rec, "IMF")
        IDV = _get_info_float(rec, "IDV")
        DP = _get_info_int(rec, "DP")

        print(
            f"[DEBUG] Considering call deletion {obs_start}-{obs_end} (len={obs_len}) "
            f"with IMF={IMF}, IDV={IDV}, DP={DP}"
        )

        for _, r in meta_gene.iterrows():
            exp_start = int(r["del_start"])
            exp_end = int(r["del_end"])
            exp_type = int(r["del_type"])
            exp_len = exp_end - exp_start + 1
            gt_IMF = float(r["gt_IMF"])
            gt_IDV = float(r["gt_IDV"])
            gt_DP = int(r["gt_DP"])

            overlap_start = max(obs_start, exp_start)
            overlap_end = min(obs_end, exp_end)
            overlap = max(0, overlap_end - overlap_start + 1)
            if overlap <= 0:
                continue

            frac_expected = overlap / exp_len
            if frac_expected < min_frac:
                print(
                    f"[DEBUG] Overlap between obs {obs_start}-{obs_end} and "
                    f"exp {exp_start}-{exp_end} (type={exp_type}) has frac_expected="
                    f"{frac_expected:.2f} < {min_frac}, skipping."
                )
                continue

            frac_observed = overlap / obs_len
            len_diff = abs(obs_len - exp_len)

            passes_thresholds = IMF >= gt_IMF and IDV >= gt_IDV and DP >= gt_DP

            score_tuple = (frac_expected, frac_observed, -len_diff, -exp_type)

            print(
                f"[DEBUG] Call candidate: obs {obs_start}-{obs_end} len={obs_len} "
                f"-> canonical {exp_start}-{exp_end} len={exp_len} (type={exp_type}), "
                f"overlap={overlap}, frac_expected={frac_expected:.3f}, "
                f"frac_observed={frac_observed:.3f}, len_diff={len_diff}, "
                f"passes_thresholds={passes_thresholds}, score_tuple={score_tuple}"
            )

            if passes_thresholds:
                # Compete for best_pass
                if best_pass_score is None or score_tuple > best_pass_score:
                    best_pass_score = score_tuple
                    best_pass_data = {
                        "type": exp_type,
                        "len_diff": len_diff,
                        "obs_start": obs_start,
                        "obs_end": obs_end,
                        "obs_len": obs_len,
                        "frac_expected": frac_expected,
                        "frac_observed": frac_observed,
                    }
                    print(
                        f"[INFO] New best PASS candidate (call): type={exp_type}, "
                        f"len_diff={len_diff}, score_tuple={score_tuple}"
                    )
            else:
                # Compete for best_fail
                if best_fail_score is None or score_tuple > best_fail_score:
                    best_fail_score = score_tuple
                    best_fail_data = {
                        "type": exp_type,
                        "len_diff": len_diff,
                        "obs_start": obs_start,
                        "obs_end": obs_end,
                        "obs_len": obs_len,
                        "frac_expected": frac_expected,
                        "frac_observed": frac_observed,
                    }
                    print(
                        f"[INFO] New best FAIL candidate (call): type={exp_type}, "
                        f"len_diff={len_diff}, score_tuple={score_tuple}"
                    )

    # Decide which candidate wins
    if best_pass_score is not None:
        t = best_pass_data["type"]
        ld = best_pass_data["len_diff"]
        os = best_pass_data["obs_start"]
        oe = best_pass_data["obs_end"]
        ol = best_pass_data["obs_len"]
        fe = best_pass_data["frac_expected"]
        fo = best_pass_data["frac_observed"]

        fe_str = f"{fe:.3f}" if fe is not None else "NA"
        fo_str = f"{fo:.3f}" if fo is not None else "NA"

        print(
            f"[INFO] Final CALL canonical assignment (threshold-passing): "
            f"type={t}, len_diff={ld}, obs={os}-{oe}, "
            f"frac_expected={fe_str}, frac_observed={fo_str}"
        )
        return t, ld, os, oe, ol, fe, fo, True

    if best_fail_score is not None:
        t = best_fail_data["type"]
        ld = best_fail_data["len_diff"]
        os = best_fail_data["obs_start"]
        oe = best_fail_data["obs_end"]
        ol = best_fail_data["obs_len"]
        fe = best_fail_data["frac_expected"]
        fo = best_fail_data["frac_observed"]

        fe_str = f"{fe:.3f}" if fe is not None else "NA"
        fo_str = f"{fo:.3f}" if fo is not None else "NA"

        print(
            f"[INFO] Final CALL canonical assignment (threshold-failing): "
            f"type={t}, len_diff={ld}, obs={os}-{oe}, "
            f"frac_expected={fe_str}, frac_observed={fo_str}"
        )
        return t, ld, os, oe, ol, fe, fo, False

    print("[INFO] No call canonical deletion passed the overlap filter.")
    return None, None, None, None, None, None, None, None


# ========================= Mpileup canonical assignment (categories 4–5) =========================


def assign_best_canonical_for_mpileup(
    variants: List[pysam.VariantRecord],
    meta_gene: pd.DataFrame,
    min_frac: float = 0.6,
) -> Tuple[
    Optional[int],
    Optional[int],
    Optional[int],
    Optional[int],
    Optional[int],
    Optional[float],
    Optional[float],
]:
    """
    Given mpileup deletion variants for a single gene and that gene's
    metafile subset (with del_start, del_end, del_type), pick the best
    canonical deletion type using the same scoring as call
    (without IMF/IDV/DP thresholds).
    """

    if not variants:
        print("[INFO] No mpileup variants provided for canonical assignment.")
        return None, None, None, None, None, None, None

    best_type: Optional[int] = None
    best_len_diff: Optional[int] = None
    best_obs_start: Optional[int] = None
    best_obs_end: Optional[int] = None
    best_obs_len: Optional[int] = None
    best_frac_expected: Optional[float] = None
    best_frac_observed: Optional[float] = None

    best_score_tuple: Optional[Tuple[float, float, int, int]] = None

    print(
        f"[INFO] Assigning best canonical deletion from {len(variants)} mpileup variants "
        f"and {len(meta_gene)} canonical windows with min_frac={min_frac}."
    )

    for rec in variants:
        span = get_deletion_span(rec)
        if span is None:
            continue
        obs_start, obs_end, obs_len = span
        print(
            f"[DEBUG] Considering mpileup deletion {obs_start}-{obs_end} (len={obs_len})"
        )

        for _, r in meta_gene.iterrows():
            exp_start = int(r["del_start"])
            exp_end = int(r["del_end"])
            exp_type = int(r["del_type"])
            exp_len = exp_end - exp_start + 1

            overlap_start = max(obs_start, exp_start)
            overlap_end = min(obs_end, exp_end)
            overlap = max(0, overlap_end - overlap_start + 1)
            if overlap <= 0:
                continue

            frac_expected = overlap / exp_len
            if frac_expected < min_frac:
                print(
                    f"[DEBUG] Overlap between obs {obs_start}-{obs_end} and "
                    f"exp {exp_start}-{exp_end} (type={exp_type}) has frac_expected="
                    f"{frac_expected:.2f} < {min_frac}, skipping."
                )
                continue

            # Reciprocal: fraction of observed deletion inside canonical window
            frac_observed = overlap / obs_len

            len_diff = abs(obs_len - exp_len)

            score_tuple = (frac_expected, frac_observed, -len_diff, -exp_type)

            print(
                f"[DEBUG] Candidate mapping: obs {obs_start}-{obs_end} len={obs_len} "
                f"-> canonical {exp_start}-{exp_end} len={exp_len} (type={exp_type}), "
                f"overlap={overlap}, frac_expected={frac_expected:.3f}, "
                f"frac_observed={frac_observed:.3f}, len_diff={len_diff}, "
                f"score_tuple={score_tuple}"
            )

            if best_score_tuple is None or score_tuple > best_score_tuple:
                best_score_tuple = score_tuple
                best_type = exp_type
                best_len_diff = len_diff
                best_obs_start = obs_start
                best_obs_end = obs_end
                best_obs_len = obs_len
                best_frac_expected = frac_expected
                best_frac_observed = frac_observed
                print(
                    f"[INFO] New best canonical candidate (mpileup): type={best_type}, "
                    f"len_diff={best_len_diff}, score_tuple={best_score_tuple}"
                )

    if best_type is None:
        print("[INFO] No mpileup canonical deletion passed the overlap filter.")
        return None, None, None, None, None, None, None

    fe_str = f"{best_frac_expected:.3f}" if best_frac_expected is not None else "NA"
    fo_str = f"{best_frac_observed:.3f}" if best_frac_observed is not None else "NA"

    print(
        f"[INFO] Final mpileup canonical assignment: type={best_type}, "
        f"len_diff={best_len_diff}, obs={best_obs_start}-{best_obs_end}, "
        f"frac_expected={fe_str}, frac_observed={fo_str}"
    )
    return (
        best_type,
        best_len_diff,
        best_obs_start,
        best_obs_end,
        best_obs_len,
        best_frac_expected,
        best_frac_observed,
    )


# ========================= SAM / assembly helpers =========================


def parse_cigar_deletions(cigar: str, pos: int) -> List[Tuple[int, int, int]]:
    """
    Parse a CIGAR string and return a list of deletion events on the reference.

    Each deletion is returned as (del_start, del_end, del_len), all 1-based
    on the reference coordinate system.

    POS is the 1-based leftmost mapping position of the read on the reference.
    """
    deletions: List[Tuple[int, int, int]] = []
    ref_pos = pos  # 1-based reference position

    for length_str, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        length = int(length_str)

        if op == "D":
            del_start = ref_pos
            del_end = ref_pos + length - 1
            deletions.append((del_start, del_end, length))
            ref_pos += length  # deletion consumes reference
        elif op in ("M", "=", "X", "N"):
            # Aligned match/mismatch or ref skip consumes reference
            ref_pos += length
        else:
            # I, S, H, P do not consume reference
            pass

    return deletions


def assign_best_canonical_for_assembly(
    del_spans: List[Tuple[int, int, int]],
    meta_gene: pd.DataFrame,
    min_frac: float = 0.6,
) -> Tuple[
    Optional[int],
    Optional[int],
    Optional[int],
    Optional[int],
    Optional[int],
    Optional[float],
    Optional[float],
]:
    """
    Canonical mapping for deletions derived from SAM CIGAR strings ("assembly").

    del_spans is a list of (obs_start, obs_end, obs_len) tuples on the reference.

    Scoring and overlap logic is identical to the mpileup-based assignment:

      1) Require frac_expected = overlap/expected_len >= min_frac.
      2) Among passing candidates, maximise:

           (frac_expected, frac_observed, -len_diff, -exp_type)
    """
    if not del_spans:
        print("[INFO] No SAM/CIGAR deletions provided for canonical assignment.")
        return None, None, None, None, None, None, None

    best_type: Optional[int] = None
    best_len_diff: Optional[int] = None
    best_obs_start: Optional[int] = None
    best_obs_end: Optional[int] = None
    best_obs_len: Optional[int] = None
    best_frac_expected: Optional[float] = None
    best_frac_observed: Optional[float] = None
    best_score_tuple: Optional[Tuple[float, float, int, int]] = None

    print(
        f"[INFO] Assigning best canonical deletion from {len(del_spans)} SAM/CIGAR deletions "
        f"and {len(meta_gene)} canonical windows with min_frac={min_frac}."
    )

    for obs_start, obs_end, obs_len in del_spans:
        print(
            f"[DEBUG] Considering SAM/CIGAR deletion {obs_start}-{obs_end} (len={obs_len})"
        )

        for _, r in meta_gene.iterrows():
            exp_start = int(r["del_start"])
            exp_end = int(r["del_end"])
            exp_type = int(r["del_type"])
            exp_len = exp_end - exp_start + 1

            overlap_start = max(obs_start, exp_start)
            overlap_end = min(obs_end, exp_end)
            overlap = max(0, overlap_end - overlap_start + 1)
            if overlap <= 0:
                continue

            frac_expected = overlap / exp_len
            if frac_expected < min_frac:
                print(
                    f"[DEBUG] Overlap between SAM obs {obs_start}-{obs_end} and "
                    f"exp {exp_start}-{exp_end} (type={exp_type}) has frac_expected="
                    f"{frac_expected:.2f} < {min_frac}, skipping."
                )
                continue

            frac_observed = overlap / obs_len
            len_diff = abs(obs_len - exp_len)

            score_tuple = (frac_expected, frac_observed, -len_diff, -exp_type)

            print(
                f"[DEBUG] SAM candidate: obs {obs_start}-{obs_end} len={obs_len} "
                f"-> canonical {exp_start}-{exp_end} len={exp_len} (type={exp_type}), "
                f"overlap={overlap}, frac_expected={frac_expected:.3f}, "
                f"frac_observed={frac_observed:.3f}, len_diff={len_diff}, "
                f"score_tuple={score_tuple}"
            )

            if best_score_tuple is None or score_tuple > best_score_tuple:
                best_score_tuple = score_tuple
                best_type = exp_type
                best_len_diff = len_diff
                best_obs_start = obs_start
                best_obs_end = obs_end
                best_obs_len = obs_len
                best_frac_expected = frac_expected
                best_frac_observed = frac_observed
                print(
                    f"[INFO] New best canonical candidate (SAM): type={best_type}, "
                    f"len_diff={best_len_diff}, score_tuple={best_score_tuple}"
                )

    if best_type is None:
        print("[INFO] No SAM/CIGAR canonical deletion passed the overlap filter.")
        return None, None, None, None, None, None, None

    fe_str = f"{best_frac_expected:.3f}" if best_frac_expected is not None else "NA"
    fo_str = f"{best_frac_observed:.3f}" if best_frac_observed is not None else "NA"

    print(
        f"[INFO] Final SAM/CIGAR canonical assignment: type={best_type}, "
        f"len_diff={best_len_diff}, obs={best_obs_start}-{best_obs_end}, "
        f"frac_expected={fe_str}, frac_observed={fo_str}"
    )
    return (
        best_type,
        best_len_diff,
        best_obs_start,
        best_obs_end,
        best_obs_len,
        best_frac_expected,
        best_frac_observed,
    )


def get_support_source(
    has_call: bool,
    has_mpileup: bool,
    has_consensus: bool,
    has_assembly: bool,
    used_consensus_for_upgrade: bool,
    used_assembly_for_upgrade: bool,
    category: Optional[str],
) -> str:
    """
    Derive a human-readable support string describing which sources
    contributed to the final category.
    """
    if category is None or category == "0":
        return "-"

    parts: List[str] = []

    # Base sources
    if has_call:
        parts.append("call")
    if has_mpileup:
        parts.append("mpileup")

    # Consensus contribution
    if has_consensus:
        if not has_call and not has_mpileup and not has_assembly and category == "6":
            parts.append("consensus_only")
        else:
            parts.append("consensus")

    # Assembly (SAM/CIGAR) contribution
    if has_assembly:
        if not has_call and not has_mpileup and not has_consensus and category == "6":
            parts.append("assembly_only")
        else:
            parts.append("assembly")

    support = "+".join(parts) if parts else "unknown"

    # Conflict flag for Category 7
    if category == "7":
        if support != "unknown":
            support = support + "+conflict"
        else:
            support = "conflict"

    return support


# ========================= Main routine =========================


def run(
    organism: str,
    bcf_path: str,
    mpileup_bcf_path: str,
    meta_path: str,
    fasta_path: str,
    output_path: str,
    deletion_region_buffer: int = 5,
    overlap_fraction: float = 0.6,
    sam_path: Optional[str] = None,
) -> None:
    # Step 1: load_deletion_metafile
    print(
        "\n# ========================= Step 1: load_deletion_metafile =========================\n"
    )
    meta_df = load_deletion_metafile(meta_path)

    # Step 2: Check_organism
    print(
        "\n# ========================= Step 2: Check_organism =========================\n"
    )
    org_meta = check_organism(meta_df, organism)

    if org_meta.empty:
        print(
            f"[WARN] No rows for organism '{organism}' in deletion metafile. Writing empty output."
        )
        out_df = pd.DataFrame(
            columns=[
                "species",
                "contig_name",
                "gene",
                "expected_start",
                "expected_end",
                "expected_variant",
                "deletion_start",
                "deletion_end",
                "deletion_length",
                "assembly_deletion_start",
                "assembly_deletion_end",
                "assembly_deletion_length",
                "expected_overlap_pct",
                "observed_overlap_pct",
                "assembly_expected_overlap_pct",
                "assembly_observed_overlap_pct",
                "consensus_N_pct",
                "consensus_N_length",
                "classified_deletion_variant",
                "classified_consensus_variant",
                "classified_assembly_variant",
                "category",
                "support_source",
            ]
        )
        out_df.to_csv(output_path, sep="\t", index=False)
        return

    print(
        f"[INFO] Using overlap_fraction={overlap_fraction} for canonical assignments."
    )

    # Load FASTA
    print("\n# ========================= Load FASTA =========================\n")
    try:
        fasta_seqs = load_fasta_sequences(fasta_path)
    except Exception as e:
        raise RuntimeError(f"Failed to open consensus FASTA '{fasta_path}': {e}")

    print("\n# ========================= Check BCF files =========================\n")

    # Open call BCF once
    try:
        bcf = pysam.VariantFile(bcf_path)
        print(f"[INFO] Successfully opened call BCF: {bcf_path}")
    except Exception as e:
        raise RuntimeError(f"Failed to open call BCF file '{bcf_path}': {e}")

    # Open mpileup BCF once
    try:
        mpileup_bcf = pysam.VariantFile(mpileup_bcf_path)
        print(f"[INFO] Successfully opened mpileup BCF: {mpileup_bcf_path}")
    except Exception as e:
        raise RuntimeError(f"Failed to open mpileup BCF file '{mpileup_bcf_path}': {e}")

    print(
        "\n# ========================= Step 3: Match_gene_bcf using the call BCF header =========================\n"
    )

    # Step 3: Match_gene_bcf using the call BCF header
    genes = sorted(org_meta["gene"].astype(str).unique())
    gene_to_contig = match_gene_bcf(bcf, genes)

    results: List[dict] = []

    # Process per gene
    for gene in genes:
        meta_g = org_meta[org_meta["gene"] == gene].copy()
        print(
            f"[INFO] ----- Processing gene '{gene}' with {len(meta_g)} deletion windows -----"
        )

        if gene not in gene_to_contig:
            print(
                f"[WARN] No contig mapping found in call BCF for gene '{gene}', "
                f"all windows will be wt/category=0 (and may be fully filtered later)."
            )
            contig_name = "-"
        else:
            contig_name = gene_to_contig[gene]
            print(f"[INFO] Using contig '{contig_name}' for gene '{gene}'.")

        # Per-gene state
        chosen_type: Optional[int] = None
        chosen_len_diff: Optional[int] = None
        chosen_obs_start: Optional[int] = None
        chosen_obs_end: Optional[int] = None
        chosen_obs_len: Optional[int] = None
        chosen_frac_expected: Optional[float] = None
        chosen_frac_observed: Optional[float] = None
        chosen_category: Optional[str] = None  # "1".."7"
        chosen_source: Optional[str] = None  # "call", "mpileup", "consensus", "assembly"

        classified_deletion_variant: Optional[int] = None   # from call/mpileup
        classified_consensus_variant: Optional[int] = None  # from consensus N%
        classified_assembly_variant: Optional[int] = None   # from SAM/CIGAR

        consensus_chosen_N_pct: Optional[float] = None
        consensus_chosen_N_length: Optional[int] = None

        # Assembly-specific canonical mapping (for reporting)
        assembly_deletion_start: Optional[int] = None
        assembly_deletion_end: Optional[int] = None
        assembly_deletion_length: Optional[int] = None
        assembly_frac_expected: Optional[float] = None
        assembly_frac_observed: Optional[float] = None

        # Support flags
        has_call = False
        has_mpileup = False
        has_consensus = False
        has_assembly = False
        used_consensus_for_upgrade = False
        used_assembly_for_upgrade = False
        has_conflict = False  # not used directly, but kept for clarity

        # ---------------------- CALL BCF (categories 1–3) ----------------------
        print(
            "\n# ========================= Step 4: CALL BCF – canonical classification (categories 1–3) =========================\n"
        )

        if contig_name != "-":
            gene_min_start = int(meta_g["del_start"].min())
            gene_max_end = int(meta_g["del_end"].max())

            print(
                f"[INFO] Call BCF: extracting candidate deletions spanning all windows for gene '{gene}' "
                f"({gene_min_start}-{gene_max_end})"
            )

            call_deletion_variants = extract_deletion_variants_in_region(
                bcf=bcf,
                contig=contig_name,
                start=gene_min_start,
                end=gene_max_end,
                region_buffer=deletion_region_buffer,
            )

            (
                call_type,
                call_len_diff,
                call_obs_start,
                call_obs_end,
                call_obs_len,
                call_frac_expected,
                call_frac_observed,
                call_passes_thresholds,
            ) = assign_best_canonical_for_call(
                variants=call_deletion_variants,
                meta_gene=meta_g,
                min_frac=overlap_fraction,
            )

            if call_type is not None:
                has_call = True
                if call_passes_thresholds:
                    if call_len_diff is not None and call_len_diff <= 1:
                        chosen_category = "1"
                        print(
                            f"[INFO] Call-based classification for gene '{gene}': "
                            f"Category 1 (high confidence; len_diff={call_len_diff})"
                        )
                    else:
                        chosen_category = "2"
                        print(
                            f"[INFO] Call-based classification for gene '{gene}': "
                            f"Category 2 (good confidence; len_diff={call_len_diff})"
                        )
                else:
                    chosen_category = "3"
                    print(
                        f"[INFO] Call-based classification for gene '{gene}': "
                        f"Category 3 (thresholds not met; len_diff={call_len_diff})"
                    )

                chosen_type = call_type
                chosen_len_diff = call_len_diff
                chosen_obs_start = call_obs_start
                chosen_obs_end = call_obs_end
                chosen_obs_len = call_obs_len
                chosen_frac_expected = call_frac_expected
                chosen_frac_observed = call_frac_observed
                chosen_source = "call"
                classified_deletion_variant = call_type
            else:
                print(
                    f"[INFO] Call BCF: no canonical deletion found for gene '{gene}' "
                    f"(no Category 1–3)."
                )
        else:
            print(
                f"[INFO] Gene '{gene}' has no contig in call BCF, skipping call-based classification."
            )

        # ---------------------- MPILEUP (categories 4–5, fallback) ----------------------
        print(
            "\n# ========================= Step 5: MPILEUP BCF – canonical classification (categories 4–5) =========================\n"
        )

        if chosen_category is None and contig_name != "-":
            print(
                f"[INFO] No call-based Category 1–3 for gene '{gene}', "
                f"checking mpileup BCF for canonical assignment (Category 4–5)."
            )

            gene_min_start = int(meta_g["del_start"].min())
            gene_max_end = int(meta_g["del_end"].max())

            print(
                f"[INFO] Mpileup BCF: extracting candidate deletions spanning all windows for gene '{gene}' "
                f"({gene_min_start}-{gene_max_end})"
            )

            mpileup_deletion_variants = extract_deletion_variants_in_region(
                bcf=mpileup_bcf,
                contig=contig_name,
                start=gene_min_start,
                end=gene_max_end,
                region_buffer=deletion_region_buffer,
            )

            (
                mp_type,
                mp_len_diff,
                mp_obs_start,
                mp_obs_end,
                mp_obs_len,
                mp_frac_expected,
                mp_frac_observed,
            ) = assign_best_canonical_for_mpileup(
                variants=mpileup_deletion_variants,
                meta_gene=meta_g,
                min_frac=overlap_fraction,
            )

            if mp_type is not None:
                has_mpileup = True
                if mp_len_diff is not None and mp_len_diff == 0:
                    chosen_category = "4"
                    print(
                        f"[INFO] Mpileup-based classification for gene '{gene}': "
                        f"Category 4 (exact canonical length; len_diff={mp_len_diff})"
                    )
                else:
                    chosen_category = "5"
                    print(
                        f"[INFO] Mpileup-based classification for gene '{gene}': "
                        f"Category 5 (inexact canonical length; len_diff={mp_len_diff})"
                    )

                chosen_type = mp_type
                chosen_len_diff = mp_len_diff
                chosen_obs_start = mp_obs_start
                chosen_obs_end = mp_obs_end
                chosen_obs_len = mp_obs_len
                chosen_frac_expected = mp_frac_expected
                chosen_frac_observed = mp_frac_observed
                chosen_source = "mpileup"
                classified_deletion_variant = mp_type
            else:
                print(
                    f"[INFO] Mpileup BCF: no canonical deletion assigned for gene '{gene}'."
                )
        elif chosen_category is not None:
            print(
                f"[INFO] Gene '{gene}' already has call-based category {chosen_category}; "
                f"mpileup canonical mapping skipped."
            )
        else:
            print(
                f"[INFO] Gene '{gene}' has no contig mapping; mpileup classification skipped."
            )

        # Decide whether to skip consensus + assembly entirely
        skip_secondary = (
            chosen_category == "1" and chosen_source == "call"
        )

        # ---------------------- Step 6: Build per-window consensus stats (if not skipped) ----------------------
        window_consensus_info: Dict[Tuple[int, int, int], dict] = {}

        if not skip_secondary:
            print(
                "\n# ========================= Step 6: Build per-window consensus stats for this gene (N% and threshold) =========================\n"
            )

            for _, row in meta_g.iterrows():
                expected_start = int(row["del_start"])
                expected_end = int(row["del_end"])
                expected_variant = int(row["del_type"])
                consensus_threshold = float(row["consensus_N"])

                consensus_N_pct = None
                consensus_N_length = None

                if contig_name != "-" and contig_name in fasta_seqs:
                    seq = fasta_seqs[contig_name]
                    start_idx = expected_start - 1
                    end_idx = expected_end
                    if start_idx < 0 or end_idx > len(seq):
                        print(
                            f"[WARN] Expected window {gene}:{expected_start}-{expected_end} "
                            f"out of FASTA bounds for contig '{contig_name}' (len={len(seq)}). "
                            f"Skipping consensus N% calculation."
                        )
                    else:
                        region = seq[start_idx:end_idx]
                        if region:
                            n_count = sum(1 for base in region if base in ("N", "n"))
                            consensus_N_length = n_count
                            consensus_N_pct = (n_count / len(region)) * 100.0
                            print(
                                f"[INFO] Consensus N% for {gene}:{expected_start}-{expected_end} on "
                                f"{contig_name} = {consensus_N_pct:.2f}% ({n_count}/{len(region)} N), "
                                f"threshold={consensus_threshold:.2f}%"
                            )
                elif contig_name == "-":
                    print(
                        f"[INFO] No contig for gene '{gene}' in FASTA; consensus_N_pct left as None."
                    )
                else:
                    print(
                        f"[WARN] Contig '{contig_name}' not found in FASTA; "
                        f"consensus_N_pct left as None."
                    )

                key = (expected_start, expected_end, expected_variant)
                window_consensus_info[key] = {
                    "consensus_N_pct": consensus_N_pct,
                    "consensus_N_length": consensus_N_length,
                    "consensus_threshold": consensus_threshold,
                }
        else:
            print(
                "\n# ========================= Skipping consensus (Step 7) and assembly (Step 8) because gene is Category 1 from CALL =========================\n"
            )

        # ---------------------- Step 7: Consensus-based classification & adjustment (if not skipped) ----------------------
        if not skip_secondary:
            print(
                "\n# ========================= Step 7: Consensus-based classification and category adjustment =========================\n"
            )

            consensus_candidates: List[Tuple[int, int, int, float, int]] = []

            for _, row in meta_g.iterrows():
                expected_start = int(row["del_start"])
                expected_end = int(row["del_end"])
                expected_variant = int(row["del_type"])
                key = (expected_start, expected_end, expected_variant)
                info = window_consensus_info.get(key, {})
                c_pct = info.get("consensus_N_pct")
                c_len = info.get("consensus_N_length")
                c_thr = info.get("consensus_threshold", 0.0)

                if c_pct is not None and c_pct >= c_thr:
                    consensus_candidates.append(
                        (expected_variant, expected_start, expected_end, c_pct, c_len or 0)
                    )

            if consensus_candidates:
                consensus_candidates.sort(key=lambda x: x[0], reverse=True)
                (
                    c_exp_variant,
                    c_exp_start,
                    c_exp_end,
                    c_pct,
                    c_n_len,
                ) = consensus_candidates[0]
                classified_consensus_variant = c_exp_variant
                consensus_chosen_N_pct = c_pct
                consensus_chosen_N_length = c_n_len
                has_consensus = True

                print(
                    f"[INFO] Consensus classification for gene '{gene}': "
                    f"selected expected_variant={c_exp_variant}, "
                    f"window={c_exp_start}-{c_exp_end}, "
                    f"consensus_N_pct={c_pct:.2f}%, "
                    f"consensus_N_length={c_n_len}"
                )
            else:
                print(
                    f"[INFO] No consensus window for gene '{gene}' passed its consensus_N threshold; "
                    f"no consensus-based deletion classification."
                )

            # Consensus rules
            if classified_deletion_variant is not None and chosen_category is not None:
                if chosen_category == "1":
                    # This case now won't happen because we skip consensus if Category 1 from CALL
                    print(
                        f"[INFO] Gene '{gene}' has Category 1 from CALL; "
                        f"consensus will NOT change this classification."
                    )
                else:
                    if classified_consensus_variant is not None:
                        if classified_consensus_variant == classified_deletion_variant:
                            expected_len = classified_deletion_variant
                            if (
                                consensus_chosen_N_length is not None
                                and consensus_chosen_N_length == expected_len
                            ):
                                boost = 2
                                support_type = "exact"
                            else:
                                boost = 1
                                support_type = "inexact"

                            original_cat_int = int(chosen_category)
                            new_cat_int = max(1, original_cat_int - boost)
                            if new_cat_int != original_cat_int:
                                used_consensus_for_upgrade = True
                            print(
                                f"[INFO] Consensus {support_type} support for gene '{gene}' "
                                f"(deletion_variant={classified_deletion_variant}): "
                                f"original_category={chosen_category} → new_category={new_cat_int}"
                            )
                            chosen_category = str(new_cat_int)
                        else:
                            print(
                                f"[INFO] Conflict for gene '{gene}': "
                                f"deletion_variant={classified_deletion_variant} from {chosen_source}, "
                                f"but consensus_variant={classified_consensus_variant}. "
                                f"Setting category=7 (conflict)."
                            )
                            chosen_category = "7"
                            has_conflict = True
                    else:
                        print(
                            f"[INFO] Gene '{gene}' has Category {chosen_category} from {chosen_source} "
                            f"but consensus did not classify a deletion; category unchanged."
                        )
            else:
                if classified_consensus_variant is not None:
                    chosen_category = "6"
                    chosen_type = classified_consensus_variant
                    chosen_source = "consensus"
                    print(
                        f"[INFO] Gene '{gene}' has no call/mpileup deletion but consensus "
                        f"classifies expected_variant={classified_consensus_variant}; "
                        f"setting category=6 (consensus-only)."
                    )
                else:
                    print(
                        f"[INFO] Gene '{gene}' has no deletion classification from call/mpileup "
                        f"and no consensus classification (all remain category 0)."
                    )
        elif skip_secondary:
            print(
                "[INFO] Gene has Category 1 from CALL; consensus classification skipped."
            )
        

        # ---------------------- Step 8: SAM / assembly-based adjustment (if not skipped) ----------------------
        sam_deletion_spans: List[Tuple[int, int, int]] = []

        if not skip_secondary and sam_path is not None and contig_name != "-":
            print(
                "\n# ========================= Step 8: SAM gene-substring matches and assembly-based category adjustment =========================\n"
            )

            print(
                f"[INFO] Scanning SAM file for gene-substring matches (genes=1): {sam_path}"
            )
            try:
                with open(sam_path, "r") as sam_fh:
                    for line in sam_fh:
                        line = line.rstrip("\n")
                        if not line or line.startswith("@"):
                            continue
                        fields = line.split("\t")
                        if len(fields) < 11:
                            continue
                        qname = fields[0]
                        flag = int(fields[1])
                        rname = fields[2]
                        pos = int(fields[3])
                        mapq = int(fields[4])
                        cigar = fields[5]

                        if rname != contig_name:
                            continue
                        if re.search(
                            re.escape(gene), rname, flags=re.IGNORECASE
                        ) is None:
                            continue

                        cs_tag = ""
                        for f in fields[11:]:
                            if f.startswith("cs:Z:"):
                                cs_tag = f
                                break

                        print(
                            f"[SAM_MATCH] gene='{gene}' RNAME='{rname}' FLAG={flag} POS={pos} "
                            f"MAPQ={mapq} CIGAR='{cigar}' CS='{cs_tag}'"
                        )

                        del_events = parse_cigar_deletions(cigar, pos)
                        for del_start, del_end, del_len in del_events:
                            print(
                                f"[SAM_DEL] gene='{gene}' RNAME='{rname}' "
                                f"DEL_START={del_start} DEL_LEN={del_len}"
                            )
                            sam_deletion_spans.append((del_start, del_end, del_len))

                print(
                    f"[INFO] Finished scanning SAM. Total matching alignments: "
                    f"{len(sam_deletion_spans)} deletions (events, across all hits)."
                )
            except Exception as e:
                print(f"[ERROR] Failed to read SAM file '{sam_path}': {e}")
        elif sam_path is None:
            print("[INFO] No SAM file provided; skipping assembly-based classification.")
        elif skip_secondary:
            print(
                "[INFO] Gene has Category 1 from CALL; assembly-based classification skipped."
            )
        else:
            print(
                f"[INFO] Gene '{gene}' has no contig mapping; skipping assembly-based classification."
            )

        if not skip_secondary and sam_deletion_spans:
            # Two modes:
            #  1) Support mode (we already have a canonical type from call/mpileup/consensus)
            #  2) Fallback mode (no canonical type at all -> assembly-only Category 6)
            expected_lengths = sorted(set(meta_g["del_type"].astype(int)))

            if chosen_type is not None and chosen_category is not None:
                # -------- SUPPORT MODE: pure length check vs chosen_type --------
                print(
                    f"[INFO] Assembly support mode for gene '{gene}': "
                    f"checking for CIGAR deletions of length {chosen_type}."
                )
                support_events = [
                    (s, e, l)
                    for (s, e, l) in sam_deletion_spans
                    if l == chosen_type
                ]

                if support_events:
                    # Choose leftmost event for reporting
                    support_events.sort(key=lambda x: x[0])
                    asm_obs_start, asm_obs_end, asm_obs_len = support_events[0]

                    has_assembly = True
                    classified_assembly_variant = chosen_type
                    assembly_deletion_start = asm_obs_start
                    assembly_deletion_end = asm_obs_end
                    assembly_deletion_length = asm_obs_len

                    # We treat assembly as full support; we do not use overlap_fraction here
                    assembly_frac_expected = None
                    assembly_frac_observed = None

                    if chosen_category != "1":
                        original_cat_int = int(chosen_category)
                        # Exact-length assembly support: boost 2 categories closer to 1
                        new_cat_int = max(1, original_cat_int - 2)
                        if new_cat_int != original_cat_int:
                            used_assembly_for_upgrade = True
                        print(
                            f"[INFO] Assembly exact-length support for gene '{gene}' "
                            f"(deletion_variant={chosen_type}): "
                            f"original_category={chosen_category} → new_category={new_cat_int}"
                        )
                        chosen_category = str(new_cat_int)
                    else:
                        print(
                            f"[INFO] Gene '{gene}' is Category 1; assembly support will not change this."
                        )
                else:
                    print(
                        f"[INFO] Assembly found no deletions of length {chosen_type} for gene '{gene}'; "
                        f"ignoring assembly for this gene (no conflict)."
                    )

            else:
                # -------- FALLBACK MODE: use mpileup-like canonical scoring --------
                print(
                    f"[INFO] No canonical deletion for gene '{gene}' from call/mpileup/consensus; "
                    f"using assembly-only canonical assignment (Category 6 fallback)."
                )
                (
                    asm_type,
                    asm_len_diff,
                    asm_obs_start,
                    asm_obs_end,
                    asm_obs_len,
                    asm_frac_expected,
                    asm_frac_observed,
                ) = assign_best_canonical_for_assembly(
                    del_spans=sam_deletion_spans,
                    meta_gene=meta_g,
                    min_frac=overlap_fraction,
                )

                if asm_type is not None:
                    has_assembly = True
                    classified_assembly_variant = asm_type
                    assembly_deletion_start = asm_obs_start
                    assembly_deletion_end = asm_obs_end
                    assembly_deletion_length = asm_obs_len
                    assembly_frac_expected = asm_frac_expected
                    assembly_frac_observed = asm_frac_observed

                    chosen_category = "6"
                    chosen_type = asm_type
                    chosen_source = "assembly"
                    chosen_len_diff = asm_len_diff
                    chosen_obs_start = asm_obs_start
                    chosen_obs_end = asm_obs_end
                    chosen_obs_len = asm_obs_len
                    chosen_frac_expected = asm_frac_expected
                    chosen_frac_observed = asm_frac_observed

                    print(
                        f"[INFO] Gene '{gene}' has assembly-only classification: "
                        f"variant={asm_type}, category=6 (assembly-only)."
                    )
                else:
                    print(
                        f"[INFO] No assembly canonical deletion passed overlap_fraction for gene '{gene}'."
                    )
        else:
            if not skip_secondary:
                print(
                    f"[INFO] No SAM/CIGAR deletions for gene '{gene}' after filtering; "
                    f"no assembly-based adjustment."
                )

        # ---------------------- support_source ----------------------
        support_source_str = get_support_source(
            has_call=has_call,
            has_mpileup=has_mpileup,
            has_consensus=has_consensus,
            has_assembly=has_assembly,
            used_consensus_for_upgrade=used_consensus_for_upgrade,
            used_assembly_for_upgrade=used_assembly_for_upgrade,
            category=chosen_category,
        )

        # Compute assembly overlap percentages for reporting (if we had canonical assembly mapping)
        assembly_expected_overlap_pct = None
        assembly_observed_overlap_pct = None
        if assembly_frac_expected is not None:
            assembly_expected_overlap_pct = assembly_frac_expected * 100.0
        if assembly_frac_observed is not None:
            assembly_observed_overlap_pct = assembly_frac_observed * 100.0

        # ---------------------- Step 9: Build per-window rows for this gene ----------------------
        print(
            "\n# ========================= Step 9: Build per-window rows for this gene =========================\n"
        )

        assigned_for_gene = False
        for _, row in meta_g.iterrows():
            species = str(row["species"])
            expected_start = int(row["del_start"])
            expected_end = int(row["del_end"])
            expected_variant = int(row["del_type"])

            deletion_start = None
            deletion_end = None
            deletion_length = None
            expected_overlap_pct = None
            observed_overlap_pct = None
            consensus_N_pct = None
            consensus_N_length = None
            category = "0"

            key = (expected_start, expected_end, expected_variant)
            c_info = window_consensus_info.get(key, {})
            if c_info:
                consensus_N_pct = c_info.get("consensus_N_pct")
                consensus_N_length = c_info.get("consensus_N_length")

            if (
                not assigned_for_gene
                and chosen_category is not None
                and chosen_type is not None
                and expected_variant == chosen_type
            ):
                category = chosen_category
                if (
                    chosen_obs_start is not None
                    and chosen_obs_end is not None
                    and chosen_obs_len is not None
                    and chosen_source in ("call", "mpileup", "assembly")
                ):
                    deletion_start = int(chosen_obs_start)
                    deletion_end = int(chosen_obs_end)
                    deletion_length = int(chosen_obs_len)

                if (
                    chosen_frac_expected is not None
                    and chosen_frac_observed is not None
                ):
                    expected_overlap_pct = chosen_frac_expected * 100.0
                    observed_overlap_pct = chosen_frac_observed * 100.0

                assigned_for_gene = True

                print(
                    f"[INFO] For gene '{gene}', assigning "
                    f"{(chosen_source or 'UNKNOWN').upper()} "
                    f"canonical type={chosen_type} to window {expected_start}-{expected_end}, "
                    f"category={category}, expected_overlap_pct={expected_overlap_pct}, "
                    f"observed_overlap_pct={observed_overlap_pct}"
                )

            results.append(
                {
                    "species": species,
                    "contig_name": contig_name,
                    "gene": gene,
                    "expected_start": expected_start,
                    "expected_end": expected_end,
                    "expected_variant": expected_variant,
                    "deletion_start": deletion_start,
                    "deletion_end": deletion_end,
                    "deletion_length": deletion_length,
                    "assembly_deletion_start": assembly_deletion_start,
                    "assembly_deletion_end": assembly_deletion_end,
                    "assembly_deletion_length": assembly_deletion_length,
                    "expected_overlap_pct": expected_overlap_pct,
                    "observed_overlap_pct": observed_overlap_pct,
                    "assembly_expected_overlap_pct": assembly_expected_overlap_pct,
                    "assembly_observed_overlap_pct": assembly_observed_overlap_pct,
                    "consensus_N_pct": consensus_N_pct,
                    "consensus_N_length": consensus_N_length,
                    "classified_deletion_variant": classified_deletion_variant,
                    "classified_consensus_variant": classified_consensus_variant,
                    "classified_assembly_variant": classified_assembly_variant,
                    "category": category,
                    "support_source": support_source_str,
                }
            )

    # Build final DataFrame
    out_df = pd.DataFrame(results)

    # Only keep rows where category != "0"
    nonzero_df = out_df[out_df["category"] != "0"].copy()

    if nonzero_df.empty:
        print("[INFO] No deletions detected (all category 0). Writing single dash row.")
        first_row = org_meta.iloc[0]
        species = str(first_row["species"])
        gene = str(first_row["gene"])
        contig_name = "-"
        if gene in gene_to_contig:
            contig_name = gene_to_contig[gene]

        fallback = pd.DataFrame(
            [
                {
                    "species": species,
                    "contig_name": contig_name,
                    "gene": gene,
                    "expected_start": "-",
                    "expected_end": "-",
                    "expected_variant": "-",
                    "deletion_start": "-",
                    "deletion_end": "-",
                    "deletion_length": "-",
                    "assembly_deletion_start": "-",
                    "assembly_deletion_end": "-",
                    "assembly_deletion_length": "-",
                    "expected_overlap_pct": "-",
                    "observed_overlap_pct": "-",
                    "assembly_expected_overlap_pct": "-",
                    "assembly_observed_overlap_pct": "-",
                    "consensus_N_pct": "-",
                    "consensus_N_length": "-",
                    "classified_deletion_variant": "-",
                    "classified_consensus_variant": "-",
                    "classified_assembly_variant": "-",
                    "category": "-",
                    "support_source": "-",
                }
            ]
        )
        fallback.to_csv(output_path, sep="\t", index=False)
    else:
        int_cols = [
            "expected_start",
            "expected_end",
            "expected_variant",
            "deletion_start",
            "deletion_end",
            "deletion_length",
            "assembly_deletion_start",
            "assembly_deletion_end",
            "assembly_deletion_length",
            "consensus_N_length",
            "classified_deletion_variant",
            "classified_consensus_variant",
            "classified_assembly_variant",
        ]
        for col in int_cols:
            if col in nonzero_df.columns:
                nonzero_df[col] = nonzero_df[col].astype("Int64")

        pct_cols = [
            "expected_overlap_pct",
            "observed_overlap_pct",
            "assembly_expected_overlap_pct",
            "assembly_observed_overlap_pct",
            "consensus_N_pct",
        ]
        for col in pct_cols:
            if col in nonzero_df.columns:
                nonzero_df[col] = nonzero_df[col].map(
                    lambda x: f"{x:.2f}" if isinstance(x, (float, int)) else x
                )

        nonzero_df = nonzero_df[
            [
                "species",
                "contig_name",
                "gene",
                "expected_start",
                "expected_end",
                "expected_variant",
                "deletion_start",
                "deletion_end",
                "deletion_length",
                "assembly_deletion_start",
                "assembly_deletion_end",
                "assembly_deletion_length",
                "expected_overlap_pct",
                "observed_overlap_pct",
                "assembly_expected_overlap_pct",
                "assembly_observed_overlap_pct",
                "consensus_N_pct",
                "consensus_N_length",
                "classified_deletion_variant",
                "classified_consensus_variant",
                "classified_assembly_variant",
                "category",
                "support_source",
            ]
        ]
        nonzero_df.to_csv(output_path, sep="\t", index=False)

    print(f"[INFO] Deletion summary written to: {output_path}")


def main() -> None:
    arg = argparse.ArgumentParser(
        description=(
            "Identify configured deletions from a call BCF and a mpileup BCF "
            "using a deletion metafile (no KMA/BED), and report %N in the "
            "consensus FASTA for each expected deletion window, with "
            "consensus- and assembly-supported category adjustments."
        )
    )
    arg.add_argument(
        "--organism", required=True, help="Organism name to filter metafile by."
    )
    arg.add_argument(
        "--call", required=True, help="BCF file with genotype calls (deletions)."
    )
    arg.add_argument(
        "--mpileup",
        required=True,
        help="BCF file with mpileup indel calls (fallback).",
    )
    arg.add_argument(
        "--metafile",
        required=True,
        help=(
            "TSV deletion metafile with columns: species, gene, del_start, "
            "del_end, del_type, gt_IMF, gt_IDV, gt_DP, consensus_N."
        ),
    )
    arg.add_argument(
        "--fsa",
        required=True,
        help=(
            "Consensus FASTA; used to compute %N for each expected deletion "
            "window defined in the metafile."
        ),
    )
    arg.add_argument(
        "--sam",
        required=False,
        help=(
            "Optional SAM alignment file; used to extract CIGAR-based deletions "
            "for additional assembly support and classification."
        ),
    )
    arg.add_argument(
        "-o",
        "--output",
        required=True,
        help=(
            "Output TSV path "
            "(species, contig_name, gene, expected_start, expected_end, "
            "expected_variant, deletion_start, deletion_end, deletion_length, "
            "assembly_deletion_start, assembly_deletion_end, assembly_deletion_length, "
            "expected_overlap_pct, observed_overlap_pct, "
            "assembly_expected_overlap_pct, assembly_observed_overlap_pct, "
            "consensus_N_pct, consensus_N_length, classified_deletion_variant, "
            "classified_consensus_variant, classified_assembly_variant, "
            "category, support_source)."
        ),
    )
    arg.add_argument(
        "--deletion_region_buffer",
        type=int,
        default=5,
        help="bp to pad each expected deletion region on each side (default 5).",
    )
    arg.add_argument(
        "--overlap_fraction",
        type=float,
        default=0.6,
        help=(
            "Minimum fraction of expected deletion covered to consider mapping "
            "from call/mpileup/SAM (default 0.6)."
        ),
    )

    args = arg.parse_args()
    run(
        organism=args.organism,
        bcf_path=args.call,
        mpileup_bcf_path=args.mpileup,
        meta_path=args.metafile,
        fasta_path=args.fsa,
        output_path=args.output,
        deletion_region_buffer=args.deletion_region_buffer,
        overlap_fraction=args.overlap_fraction,
        sam_path=args.sam,
    )


if __name__ == "__main__":
    main()
