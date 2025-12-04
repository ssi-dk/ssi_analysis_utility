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


# ========================= Category → support_source helper =========================


def get_support_source(category: str) -> str:
    """
    Map a category code to a human-readable support source label.

    Semantics (fallback when no per-gene label is set):

      1,2,3 → "call"
      4,5   → "mpileup"
      6     → "consensus_only"   (no call/mpileup or consensus-only classification)
      7     → "conflict"         (call/mpileup vs consensus classification differ)

    '-' or '0' or anything else → "none" / "unknown"
    """
    if pd.isna(category):
        return "none"
    cat_str = str(category)

    if cat_str in {"1", "2", "3"}:
        return "call"
    if cat_str in {"4", "5"}:
        return "mpileup"
    if cat_str == "6":
        return "consensus_only"
    if cat_str == "7":
        return "conflict"
    if cat_str in {"0", "-"}:
        return "none"
    return "unknown"


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

         This is expected-centric: it is optimised for
         "how well does this observed deletion fulfil the canonical window
         we defined?", then refined by how well the observed event is
         contained by that canonical window.

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

    If no candidate passes min_frac, all values are None.
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

            # Scoring scheme:
            #
            # 1) Prefer higher frac_expected → more of the canonical region is deleted.
            # 2) If tie, prefer higher frac_observed → more of the observed deletion
            #    lies inside the canonical region.
            # 3) If tie, prefer smaller abs(len_diff) = |obs_len - exp_len|.
            # 4) If still tie, prefer smaller canonical length (exp_type),
            #    e.g. 39 over 54.
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

    # Decide which candidate wins overall according to the category rules:
    #
    # Category 1 / 2 prefer "best_pass" if it exists,
    # otherwise Category 3 may use "best_fail".
    if best_pass_score is not None:
        # Use best_pass_data
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
        # Use best_fail_data
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
    canonical deletion type.

    Scoring rules (same as for call, but without IMF/IDV/DP):

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

         This is expected-centric: it is optimised for
         "how well does this observed deletion fulfil the canonical window
         we defined?", then refined by how well the observed event is
         contained by that canonical window.

    Returns:
      (best_del_type,
       best_len_diff,
       best_obs_start,
       best_obs_end,
       best_obs_len,
       best_frac_expected,
       best_frac_observed)

      If no candidate passes min_frac, all values are None.
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

    # score_tuple = (frac_expected, frac_observed, -len_diff, -exp_type)
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

            # Reciprocal: what fraction of the observed deletion lies inside
            # the canonical (expected) window?
            frac_observed = overlap / obs_len

            len_diff = abs(obs_len - exp_len)

            # Scoring scheme:
            #
            # 1) Prefer higher frac_expected
            # 2) If tie, prefer higher frac_observed
            # 3) If tie, prefer smaller abs(len_diff)
            # 4) If still tie, prefer smaller canonical length (exp_type)
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

    # Prepare nicely formatted strings for logging
    if best_frac_expected is not None:
        fe_str = f"{best_frac_expected:.3f}"
    else:
        fe_str = "NA"

    if best_frac_observed is not None:
        fo_str = f"{best_frac_observed:.3f}"
    else:
        fo_str = "NA"

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


# ========================= Main routine =========================


def run(
    organism: str,
    bcf_path: str,
    mpileup_bcf_path: str,
    meta_path: str,
    fasta_path: str,  # now used to compute N%
    output_path: str,
    deletion_region_buffer: int = 5,
    overlap_fraction: float = 0.6,
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
                "expected_overlap_pct",
                "observed_overlap_pct",
                "consensus_N_pct",
                "consensus_N_length",
                "classified_deletion_variant",
                "classified_consensus_variant",
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

    # Process per gene so we produce at most one canonical call per gene
    # from call BCF (categories 1–3) and, if needed, from mpileup (4–5),
    # and possibly override / upgrade with consensus classification (including category 6/7).
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

        # Prepare storage for the final chosen assignment for this gene
        chosen_type: Optional[int] = None
        chosen_len_diff: Optional[int] = None
        chosen_obs_start: Optional[int] = None
        chosen_obs_end: Optional[int] = None
        chosen_obs_len: Optional[int] = None
        chosen_frac_expected: Optional[float] = None
        chosen_frac_observed: Optional[float] = None
        chosen_category: Optional[str] = None  # "1".."7"
        chosen_source: Optional[str] = None  # "call" or "mpileup" or "consensus"

        classified_deletion_variant: Optional[int] = None
        classified_consensus_variant: Optional[int] = None
        consensus_chosen_N_pct: Optional[float] = None
        consensus_chosen_N_length: Optional[int] = None

        # Track support source for this gene (to write into all rows for this gene)
        support_source_for_gene: Optional[str] = None
        # Track whether consensus actually adjusted the category (2–5 → better)
        consensus_adjusted: bool = False

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
                # Decide category from call:
                #   Category 1 (high confidence, exact-ish):
                #     - passes thresholds AND len_diff <= 1
                #   Category 2 (good confidence, inexact length):
                #     - passes thresholds AND len_diff > 1
                #   Category 3 (low support):
                #     - fails thresholds (any len_diff)
                if call_passes_thresholds:
                    if call_len_diff is not None and call_len_diff <= 1:
                        chosen_category = "1"  # Category 1: high confidence call
                        print(
                            f"[INFO] Call-based classification for gene '{gene}': "
                            f"Category 1 (high confidence; len_diff={call_len_diff})"
                        )
                    else:
                        chosen_category = "2"  # Category 2: good confidence call
                        print(
                            f"[INFO] Call-based classification for gene '{gene}': "
                            f"Category 2 (good confidence; len_diff={call_len_diff})"
                        )
                else:
                    chosen_category = "3"  # Category 3: weak call, below IMF/IDV/DP
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
                support_source_for_gene = "call"
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
                # Category 4: exact canonical length (len_diff == 0)
                # Category 5: inexact canonical length (len_diff > 0)
                if mp_len_diff is not None and mp_len_diff == 0:
                    chosen_category = "4"  # Category 4: mpileup exact canonical
                    print(
                        f"[INFO] Mpileup-based classification for gene '{gene}': "
                        f"Category 4 (exact canonical length; len_diff={mp_len_diff})"
                    )
                else:
                    chosen_category = "5"  # Category 5: mpileup inexact canonical
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
                support_source_for_gene = "mpileup"
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
            # contig_name == "-"
            print(
                f"[INFO] Gene '{gene}' has no contig mapping; mpileup classification skipped."
            )

        # ---------------------- Step 6: Compute consensus N% per window (for this gene) ----------------------
        print(
            "\n# ========================= Step 6: Build per-window consensus stats for this gene (N% and threshold) =========================\n"
        )

        # Store consensus stats keyed by (start, end, del_type)
        window_consensus_info: Dict[Tuple[int, int, int], dict] = {}

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

        # ---------------------- Step 7: Consensus-based classification and category adjustment ----------------------
        print(
            "\n# ========================= Step 7: Consensus-based classification and category adjustment =========================\n"
        )

        # 1) Derive a per-gene "classified_consensus_variant" by:
        #    - keeping only windows where consensus_N_pct >= consensus_N threshold
        #    - if multiple windows pass, pick the one with the largest expected_variant
        #      (i.e. the longest canonical deletion).
        consensus_candidates: List[Tuple[int, int, int, float, int]] = []
        # elements: (expected_variant, expected_start, expected_end, consensus_N_pct, consensus_N_length)

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
            # Pick the candidate with the largest expected_variant (longest canonical deletion)
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

        # 2) Adjust category based on consensus rules:
        #
        #   Rule 1: If classification is Category 1 (exact-ish high-confidence call),
        #           we do NOT let consensus change the category.
        #
        #   Rule 2: If we have a Category 2–5 from call/mpileup AND a consensus classification:
        #
        #       - If classified_consensus_variant == classified_deletion_variant:
        #           * Let expected_len = classified_deletion_variant (canonical length).
        #           * If consensus_N_length == expected_len (100% N in that window):
        #                 → "exact" consensus support → new_category = category - 2
        #           * Else (same variant but inexact N coverage):
        #                 → "inexact" consensus support → new_category = category - 1
        #           * Category is not allowed to go below 1.
        #
        #       - If classified_consensus_variant != classified_deletion_variant:
        #           → Category 7 (conflict between read-based and consensus-based deletion).
        #
        #   Rule 3:
        #       - If there is NO Category 1–5 (no call/mpileup deletion),
        #         but consensus has a classification → Category 6 (consensus-only).
        #
        # NOTE: Category 6 means "consensus-only"; Category 7 means "conflict".
        if classified_deletion_variant is not None and chosen_category is not None:
            # We have a deletion classification (Categories 1–5)
            if chosen_category == "1":
                # Category 1: high confidence from call; do not change.
                print(
                    f"[INFO] Gene '{gene}' has Category 1 from CALL; "
                    f"consensus will NOT change this classification."
                )
                # support_source_for_gene stays "call"
            else:
                # Categories 2–5: allow consensus to upgrade or set conflict.
                if classified_consensus_variant is not None:
                    if classified_consensus_variant == classified_deletion_variant:
                        # Consensus supports the SAME canonical deletion type
                        expected_len = classified_deletion_variant  # canonical bp length
                        if (
                            consensus_chosen_N_length is not None
                            and consensus_chosen_N_length == expected_len
                        ):
                            # Exact consensus support: boost by 2 categories
                            boost = 2
                            support_type = "exact"
                        else:
                            # Inexact but supportive: boost by 1 category
                            boost = 1
                            support_type = "inexact"

                        original_cat_int = int(chosen_category)
                        new_cat_int = max(1, original_cat_int - boost)
                        print(
                            f"[INFO] Consensus {support_type} support for gene '{gene}' "
                            f"(deletion_variant={classified_deletion_variant}): "
                            f"original_category={chosen_category} → new_category={new_cat_int}"
                        )
                        chosen_category = str(new_cat_int)
                        consensus_adjusted = True
                        # If the original source was call/mpileup, tag as "+consensus"
                        if chosen_source in ("call", "mpileup"):
                            support_source_for_gene = f"{chosen_source}+consensus"
                        else:
                            # very unlikely, but just in case
                            support_source_for_gene = "consensus"
                    else:
                        # Consensus calls a different canonical variant → conflict
                        print(
                            f"[INFO] Conflict for gene '{gene}': "
                            f"deletion_variant={classified_deletion_variant} from {chosen_source}, "
                            f"but consensus_variant={classified_consensus_variant}. "
                            f"Setting category=7 (conflict)."
                        )
                        chosen_category = "7"
                        support_source_for_gene = "conflict"
                else:
                    print(
                        f"[INFO] Gene '{gene}' has Category {chosen_category} from {chosen_source} "
                        f"but consensus did not classify a deletion; category unchanged."
                    )
                    # support_source_for_gene stays "call" / "mpileup"
        else:
            # No deletion from call/mpileup (no Category 1–5)
            if classified_consensus_variant is not None:
                # Consensus-only classification → Category 6
                chosen_category = "6"
                chosen_type = classified_consensus_variant
                chosen_source = "consensus"
                support_source_for_gene = "consensus_only"
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
                # support_source_for_gene stays None

        # ---------------------- Step 8: Build per-window rows for this gene ----------------------
        print(
            "\n# ========================= Step 8: Build per-window rows for this gene =========================\n"
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

            # Get precomputed consensus stats for this window
            key = (expected_start, expected_end, expected_variant)
            c_info = window_consensus_info.get(key, {})
            if c_info:
                consensus_N_pct = c_info.get("consensus_N_pct")
                consensus_N_length = c_info.get("consensus_N_length")

            # If we have a chosen canonical assignment (from call, mpileup or consensus-only),
            # assign it to exactly ONE row whose del_type == chosen_type.
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
                    and chosen_source in ("call", "mpileup")
                ):
                    # We have real coordinates from BCF
                    deletion_start = int(chosen_obs_start)
                    deletion_end = int(chosen_obs_end)
                    deletion_length = int(chosen_obs_len)
                else:
                    # consensus-only classification: keep deletion_* as None
                    pass

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
                    "expected_overlap_pct": expected_overlap_pct,
                    "observed_overlap_pct": observed_overlap_pct,
                    "consensus_N_pct": consensus_N_pct,
                    "consensus_N_length": consensus_N_length,
                    "classified_deletion_variant": classified_deletion_variant,
                    "classified_consensus_variant": classified_consensus_variant,
                    "category": category,
                    "support_source": support_source_for_gene,
                }
            )

    # Build final DataFrame
    out_df = pd.DataFrame(results)

    # Only keep rows where category != "0"
    nonzero_df = out_df[out_df["category"] != "0"].copy()

    if nonzero_df.empty:
        print("[INFO] No deletions detected (all category 0). Writing single dash row.")
        # Use first gene/species as context if possible
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
                    "expected_overlap_pct": "-",
                    "observed_overlap_pct": "-",
                    "consensus_N_pct": "-",
                    "consensus_N_length": "-",
                    "classified_deletion_variant": "-",
                    "classified_consensus_variant": "-",
                    "category": "-",
                    "support_source": "none",
                }
            ]
        )
        fallback.to_csv(output_path, sep="\t", index=False)
    else:
        # Make sure coordinate and length columns are integers (no .0)
        int_cols = [
            "expected_start",
            "expected_end",
            "expected_variant",
            "deletion_start",
            "deletion_end",
            "deletion_length",
            "consensus_N_length",
            "classified_deletion_variant",
            "classified_consensus_variant",
        ]
        for col in int_cols:
            if col in nonzero_df.columns:
                nonzero_df[col] = nonzero_df[col].astype("Int64")

        # Format percentage columns as float strings with 2 decimals
        pct_cols = [
            "expected_overlap_pct",
            "observed_overlap_pct",
            "consensus_N_pct",
        ]
        for col in pct_cols:
            if col in nonzero_df.columns:
                nonzero_df[col] = nonzero_df[col].map(
                    lambda x: f"{x:.2f}" if isinstance(x, (float, int)) else x
                )

        # Fill in missing support_source (if any) from category using the helper
        if "support_source" not in nonzero_df.columns:
            nonzero_df["support_source"] = nonzero_df["category"].map(
                get_support_source
            )
        else:
            nonzero_df["support_source"] = nonzero_df.apply(
                lambda r: r["support_source"]
                if isinstance(r["support_source"], str)
                and r["support_source"] not in ("", "none", "unknown", "nan")
                else get_support_source(r["category"]),
                axis=1,
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
                "expected_overlap_pct",
                "observed_overlap_pct",
                "consensus_N_pct",
                "consensus_N_length",
                "classified_deletion_variant",
                "classified_consensus_variant",
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
            "consensus-supported category adjustments."
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
        "-o",
        "--output",
        required=True,
        help=(
            "Output TSV path "
            "(species, contig_name, gene, expected_start, expected_end, "
            "expected_variant, deletion_start, deletion_end, deletion_length, "
            "expected_overlap_pct, observed_overlap_pct, consensus_N_pct, "
            "consensus_N_length, classified_deletion_variant, "
            "classified_consensus_variant, category, support_source)."
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
            "from call/mpileup (default 0.6)."
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
    )


if __name__ == "__main__":
    main()
