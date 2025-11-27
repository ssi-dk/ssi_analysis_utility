#!/usr/bin/env python3

from __future__ import annotations

import argparse
import re
from typing import Dict, List, Optional, Tuple

import pandas as pd
import pysam


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


# ========================= Step 5: Filter deletions by IMF/IDV/DP =========================


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


def filter_deletions_by_thresholds(
    variants: List[pysam.VariantRecord],
    gt_IMF: float,
    gt_IDV: float,
    gt_DP: int,
) -> List[pysam.VariantRecord]:
    """
    Filter deletion variants to only those that pass IMF, IDV, DP thresholds.
    """

    filtered: List[pysam.VariantRecord] = []
    print(
        f"[INFO] Filtering {len(variants)} deletion records by "
        f"IMF>={gt_IMF}, IDV>={gt_IDV}, DP>={gt_DP}"
    )
    for rec in variants:
        IMF = _get_info_float(rec, "IMF")
        IDV = _get_info_float(rec, "IDV")
        DP = _get_info_int(rec, "DP")
        print(
            f"[DEBUG] Deletion record at {rec.pos} with IMF={IMF}, "
            f"IDV={IDV}, DP={DP} vs thresholds IMF>={gt_IMF}, "
            f"IDV>={gt_IDV}, DP>={gt_DP}"
        )
        if IMF >= gt_IMF and IDV >= gt_IDV and DP >= gt_DP:
            filtered.append(rec)
            print(
                f"[DEBUG] Record at {rec.pos}: {rec.ref} -> {rec.alts[0]} PASSED thresholds"
            )
        else:
            print(
                f"[DEBUG] Record at {rec.pos}: {rec.ref} -> {rec.alts[0]} FAILED thresholds"
            )

    print(f"[INFO] Records remaining after threshold filter: {len(filtered)}")
    return filtered


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


# ========================= Step 6: Pick best deletion (length only) =========================


def pick_best_deletion(
    variants: List[pysam.VariantRecord],
) -> Optional[Tuple[int, int, int]]:
    """
    Pick the 'best' deletion from a list of deletion records.

    For now, choose the deletion with the largest length (len(REF) - len(ALT))
    considering the first deletion-like ALT in each record.

    Returns (start, end, length) of the best deletion, or None if no
    deletions could be parsed.
    """

    best_span: Optional[Tuple[int, int, int]] = None

    for rec in variants:
        span = get_deletion_span(rec)
        if span is None:
            continue
        del_start, del_end, del_len = span
        if best_span is None or del_len > best_span[2]:
            best_span = span
            print(
                f"[DEBUG] Chose deletion at record {rec.pos}: "
                f"{rec.ref} -> {rec.alts[0]} (len={del_len}) as current best."
            )

    if best_span is None:
        print("[INFO] No valid deletion length could be determined from records.")
    else:
        print(
            f"[INFO] Best deletion span selected: {best_span[0]}-{best_span[1]} "
            f"(len={best_span[2]})"
        )

    return best_span


# ========================= Step 7: Mpileup canonical assignment =========================

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

    Scoring rules:

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
            #      → more of the canonical region is deleted.
            #
            # 2) If tie, prefer higher frac_observed
            #      → more of the observed deletion lies inside the canonical region.
            #
            # 3) If tie, prefer smaller abs(len_diff) = |obs_len - exp_len|.
            #
            # 4) If still tie, prefer smaller canonical length (exp_type),
            #      e.g. 39 over 54.
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
                    f"[INFO] New best canonical candidate: type={best_type}, "
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
    fasta_path: str,
    output_path: str,
    deletion_region_buffer: int = 5,
    overlap_fraction: float = 0.6,
) -> None:
    # Step 1: load_deletion_metafile
    print(f"\n# ========================= Step 1: load_deletion_metafile =========================\n")

    meta_df = load_deletion_metafile(meta_path)

    # Step 2: Check_organism
    print(f"\n# ========================= Step 2: Check_organism =========================\n")

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
                "category",
            ]
        )
        out_df.to_csv(output_path, sep="\t", index=False)
        return

    print(
        f"[INFO] Using overlap_fraction={overlap_fraction} for canonical mpileup assignment."
    )

    print(f"\n# ========================= Check BCF files =========================\n")

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

    print(f"\n# ========================= Step 3: Match_gene_bcf using the call BCF header =========================\n")

    # Step 3: Match_gene_bcf using the call BCF header
    genes = sorted(org_meta["gene"].astype(str).unique())
    gene_to_contig = match_gene_bcf(bcf, genes)

    results: List[dict] = []

    # Process per gene so mpileup produces only one canonical call per gene
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

        gene_rows: List[dict] = []
        any_category1 = False

        # --- Call BCF per window (category 1) ---
        for _, row in meta_g.iterrows():
            species = str(row["species"])
            expected_start = int(row["del_start"])
            expected_end = int(row["del_end"])
            expected_variant = int(row["del_type"])
            gt_IMF = float(row["gt_IMF"])
            gt_IDV = float(row["gt_IDV"])
            gt_DP = int(row["gt_DP"])

            print(f"\n# ========================= Extracting bcf variants =========================\n")
            print(
                f"[INFO] Call BCF: Processing window {gene}:{expected_start}-{expected_end}, "
                f"expected_len={expected_variant}"
            )

            if contig_name == "-":
                # No contig match; forced wt
                variant = "wt"
                category = "0"
                deletion_start = None
                deletion_end = None
                deletion_length = None
            else:
                print(f"\n# ========================= Step 4: Extract deletion variants overlapping the region =========================\n")
                # Step 4: Extract deletion variants overlapping the region
                deletion_variants = extract_deletion_variants_in_region(
                    bcf=bcf,
                    contig=contig_name,
                    start=expected_start,
                    end=expected_end,
                    region_buffer=deletion_region_buffer,
                )
                print(f"\n# ========================= Step 5: Filter by IMF/IDV/DP thresholds =========================\n")
                # Step 5: Filter by IMF/IDV/DP thresholds
                deletion_variants_filtered = filter_deletions_by_thresholds(
                    variants=deletion_variants,
                    gt_IMF=gt_IMF,
                    gt_IDV=gt_IDV,
                    gt_DP=gt_DP,
                )
                print(f"\n# ========================= Step 6: Pick best deletion and set variant/category =========================\n")
                # Step 6: Pick best deletion and set variant/category
                best_span = pick_best_deletion(deletion_variants_filtered)

                if best_span is not None:
                    deletion_start, deletion_end, deletion_length = best_span
                    variant = str(deletion_length)
                    category = "1"
                    any_category1 = True
                    print(
                        f"[INFO] Call BCF: detected deletion len={deletion_length} "
                        f"for {gene}:{expected_start}-{expected_end}, category=1"
                    )
                else:
                    variant = "wt"
                    category = "0"
                    deletion_start = None
                    deletion_end = None
                    deletion_length = None
                    print(
                        f"[INFO] Call BCF: no passing deletion for {gene}:"
                        f"{expected_start}-{expected_end}, category=0"
                    )

            gene_rows.append(
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
                    "expected_overlap_pct": None,  # will be filled for mpileup mapping
                    "observed_overlap_pct": None,  # will be filled for mpileup mapping
                    "category": category,
                }
            )

        print(f"\n# ========================= Step 7: Mpileup-based canonical assignment (categories 2 & 3) =========================\n")
        # --- Mpileup-based canonical assignment (categories 2 & 3) ---
        if not any_category1 and contig_name != "-":
            print(
                f"[INFO] No category 1 deletion for gene '{gene}', "
                f"checking mpileup BCF for canonical assignment."
            )
            gene_min_start = int(meta_g["del_start"].min())
            gene_max_end = int(meta_g["del_end"].max())

            mpileup_deletion_variants = extract_deletion_variants_in_region(
                bcf=mpileup_bcf,
                contig=contig_name,
                start=gene_min_start,
                end=gene_max_end,
                region_buffer=deletion_region_buffer,
            )

            (
                best_type,
                best_len_diff,
                best_obs_start,
                best_obs_end,
                best_obs_len,
                best_frac_expected,
                best_frac_observed,
            ) = assign_best_canonical_for_mpileup(
                variants=mpileup_deletion_variants,
                meta_gene=meta_g,
                min_frac=overlap_fraction,
            )

            if best_type is not None:
                # category 2: exact canonical length (len_diff == 0)
                # category 3: mapped to canonical via scoring (len_diff > 0)
                if best_len_diff == 0:
                    chosen_category = "2"
                else:
                    chosen_category = "3"

                print(
                    f"[INFO] Mpileup canonical result for gene '{gene}': "
                    f"type={best_type}, len_diff={best_len_diff}, category={chosen_category}"
                )

                # Compute overlap percentages for the chosen canonical window
                # relative to its own expected_start/expected_end.
                overlap_pct_expected = None
                overlap_pct_observed = None

                if (
                    best_obs_start is not None
                    and best_obs_end is not None
                    and best_obs_len is not None
                    and best_frac_expected is not None
                    and best_frac_observed is not None
                ):
                    overlap_pct_expected = best_frac_expected * 100.0
                    overlap_pct_observed = best_frac_observed * 100.0

                # Assign to exactly one row with this del_type if it is still category 0
                assigned = False
                for row_dict in gene_rows:
                    if (
                        row_dict["category"] == "0"
                        and int(row_dict["expected_variant"]) == int(best_type)
                    ):
                        row_dict["variant"] = str(best_type)
                        row_dict["category"] = chosen_category
                        row_dict["deletion_start"] = best_obs_start
                        row_dict["deletion_end"] = best_obs_end
                        row_dict["deletion_length"] = best_obs_len
                        row_dict["expected_overlap_pct"] = overlap_pct_expected
                        row_dict["observed_overlap_pct"] = overlap_pct_observed
                        assigned = True
                        print(
                            f"[INFO] Assigned mpileup canonical type={best_type} to window "
                            f"{gene}:{row_dict['expected_start']}-"
                            f"{row_dict['expected_end']} "
                            f"with category={chosen_category}, "
                            f"expected_overlap_pct={overlap_pct_expected}, "
                            f"observed_overlap_pct={overlap_pct_observed}"
                        )
                        break

                if not assigned:
                    print(
                        f"[WARN] Mpileup canonical type={best_type} could not be assigned to any "
                        f"category 0 window for gene '{gene}'."
                    )
            else:
                print(
                    f"[INFO] Mpileup BCF: no canonical deletion assigned for gene '{gene}'."
                )

        elif contig_name != "-":
            print(
                f"[INFO] Gene '{gene}' already has at least one category 1 call; "
                f"mpileup canonical mapping skipped."
            )

        # Append per-gene rows to global result list
        results.extend(gene_rows)

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
                    "category": "-",
                }
            ]
        )
        fallback.to_csv(output_path, sep="\t", index=False)
    else:
        # Make sure columns are ordered / stringified nicely
        nonzero_df["expected_overlap_pct"] = nonzero_df["expected_overlap_pct"].map(
            lambda x: f"{x:.2f}" if isinstance(x, (float, int)) else x
        )
        nonzero_df["observed_overlap_pct"] = nonzero_df["observed_overlap_pct"].map(
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
                "expected_overlap_pct",
                "observed_overlap_pct",
                "category",
            ]
        ]
        nonzero_df.to_csv(output_path, sep="\t", index=False)

    print(f"[INFO] Deletion summary written to: {output_path}")


def main() -> None:
    arg = argparse.ArgumentParser(
        description="Identify configured deletions from a call BCF and a mpileup BCF using a deletion metafile (no KMA/BED)."
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
        help="Consensus FASTA (currently unused but kept for future functionality).",
    )
    arg.add_argument(
        "-o",
        "--output",
        required=True,
        help=(
            "Output TSV path "
            "(species, contig_name, gene, expected_start, expected_end, "
            "expected_variant, deletion_start, deletion_end, deletion_length, "
            "expected_overlap_pct, observed_overlap_pct, category)."
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
            "from mpileup (default 0.6)."
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
