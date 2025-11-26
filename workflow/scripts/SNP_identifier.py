#!/usr/bin/env python3

from __future__ import annotations

import argparse
import re
from typing import Dict, List, Optional

import pandas as pd
import pysam


# ========================= Step 1: load_metafile =========================

def load_metafile(meta_path: str) -> pd.DataFrame:
    """
    Load the metafile describing SNPs.

    Expected columns:
      species, gene, position, reference, alternative, gt_DP
    """
    df = pd.read_csv(meta_path, sep="\t")
    required_cols = {
        "species",
        "gene",
        "position",
        "reference",
        "alternative",
        "gt_DP",
    }

    missing_cols = required_cols - set(df.columns)
    if missing_cols:
        raise ValueError(
            f"Metafile is missing required columns: {', '.join(sorted(missing_cols))}"
        )

    # Normalize types
    df["species"] = df["species"].astype(str)
    df["gene"] = df["gene"].astype(str)
    df["position"] = df["position"].astype(int)
    df["reference"] = df["reference"].astype(str)
    df["alternative"] = df["alternative"].astype(str)
    df["gt_DP"] = df["gt_DP"].astype(int)

    return df


# ========================= Step 2: Check_organism =========================


def check_organism(meta_df: pd.DataFrame, organism: str) -> pd.DataFrame:
    """
    Filter metafile rows to the requested organism. In case of metafile containing multiple species sharing loci name

    Returns a new DataFrame containing only rows where species == organism.
    """
    organism = str(organism).strip()
    filtered_metadf = meta_df[meta_df["species"].str.strip() == organism].copy()
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

    for gene in genes:
        pattern = re.compile(re.escape(gene), flags=re.IGNORECASE)
        matched_contig: Optional[str] = None
        for contig in contig_names:
            if pattern.search(contig):
                matched_contig = contig
                print(f"gene from metafile {re.escape(gene)} was matched with contig names from bcf header {contig}")
                break

        if matched_contig is None:
            print(f"[WARN] No contig in BCF matches gene '{gene}' (substring regex).")
            continue

        gene_to_contig[gene] = matched_contig

    return gene_to_contig


# ========================= Step 4: Extract_pos_variants =========================

def extract_pos_variants(
    bcf: pysam.VariantFile,
    contig: str,
    position: int,
    snp_region_buffer: int,
) -> List[pysam.VariantRecord]:
    """
    Extract all variant records from `bcf` on `contig` within
    [position - buffer, position + buffer] (1-based, inclusive).
    """
    # Convert 1-based desired window to 0-based half-open coordinates for pysam
    start_1based = max(1, position - snp_region_buffer)
    end_1based = position + snp_region_buffer
    start0 = start_1based - 1  # inclusive
    end0 = end_1based          # pysam end is exclusive, still in 1-based

    records: List[pysam.VariantRecord] = []
    try:
        for rec in bcf.fetch(contig, start0, end0):
            print(f"Could fetch the bcf records {rec.pos} searching around {position} ± {snp_region_buffer}")
            records.append(rec)
    except ValueError as e:
        # E.g. contig not found in index
        print(f"[ERROR] Could not fetch {contig}:{start_1based}-{end_1based} from BCF: {e}")

    return records


# =========== Extract DP value for filtering between Step 4 and Steps 5/6 ===========

def get_dp_from_record(rec: pysam.VariantRecord) -> int:
    """
    Extract DP from a VariantRecord's INFO field.

    Handles DP being stored as a single int or as a tuple/list.
    Returns 0 if DP is missing or malformed.
    """
    dp = rec.info.get("DP")
    if isinstance(dp, int):
        return dp
    if isinstance(dp, (list, tuple)):
        if not dp:
            return 0
        # some BCFs store as a 1-element tuple
        val = dp[0]
        return int(val) if isinstance(val, int) else 0
    return 0


# ========================= Step 5: Filter positions by DP =========================


def filter_variants_by_dp(
    variants: List[pysam.VariantRecord],
    min_dp: int,
) -> List[pysam.VariantRecord]:
    """
    Filter variants to only those with INFO/DP >= min_dp.
    """
    filtered: List[pysam.VariantRecord] = []
    for rec in variants:
        dp_val = get_dp_from_record(rec)
        print(f"Fetch the bcf records {rec.pos} with reference {rec.ref} and alternative {rec.alts[0]} and DP value {dp_val}")
        if dp_val >= min_dp:
            filtered.append(rec)
    return filtered

# ========================= Step 6: Check deletions spanning position =========================


def check_pos_deletions(
    variants: List[pysam.VariantRecord],
    position: int,
) -> Optional[str]:
    """
    Check whether any variant in `variants` represents a deletion that spans
    `position` (1-based).

    A deletion is defined as an ALT allele shorter than REF. If any such
    deletion covers `position` (rec.pos <= position <= rec.pos + len(REF) - 1),
    returns a string like 'Δ117'. Otherwise returns None.
    """
    for rec in variants:
        ref = rec.ref
        if not rec.alts:
            continue

        # Check if any ALT is shorter than REF → deletion
        has_deletion_alt = any(
            alt is not None and len(alt) < len(ref) for alt in rec.alts
        )
        if not has_deletion_alt:
            continue

        deletion_start = rec.pos  # 1-based
        deletion_end = rec.pos + len(ref) - 1  # inclusive

        if deletion_start <= position <= deletion_end:
            print(f"The bcf records {rec.pos} with reference {rec.ref} and alternative {rec.alts[0]} is determined as a deletion")
            return f"Δ{position}"

    return None


# ========================= Step 7: Check SNP at position =========================


def check_pos_snp(
    variants: List[pysam.VariantRecord],
    position: int,
) -> Optional[str]:
    """
    Check whether there is a SNP exactly at `position` (1-based) in `variants`.

    - Only considers simple SNPs: len(REF) == 1 and len(ALT) == 1.
    - Returns variant in the format '{ref}{pos}{alt}' using the BCF orientation.

    If no SNP is present, returns None.
    """
    for rec in variants:
        if rec.pos != position:
            continue
        if not rec.alts:
            continue

        ref = rec.ref
        alt = rec.alts[0]  # first ALT only

        # Only treat pure SNPs
        if alt is None or len(ref) != 1 or len(alt) != 1:
            print(f"The bcf records {rec.pos} with reference {rec.ref} and alternative {rec.alts[0]} is determined as a SNP")
            continue

        return f"{ref}{position}{alt}"

    return None


# ========================= Main routine =========================

def run(
    organism: str,
    bcf_path: str,
    meta_path: str,
    output_path: str,
    snp_region_buffer: int = 1,
) -> None:
    # Step 1: load_metafile
    meta_df = load_metafile(meta_path)

    # Step 2: Check_organism
    org_meta = check_organism(meta_df, organism)

    if org_meta.empty:
        print(f"[WARN] No rows for organism '{organism}' in metafile. Writing empty output.")
        out_df = pd.DataFrame(columns=["species", "contig_name", "gene", "position", "variant"])
        out_df.to_csv(output_path, sep="\t", index=False)
        return

    # Open BCF once
    try:
        bcf = pysam.VariantFile(bcf_path)
    except Exception as e:
        raise RuntimeError(f"Failed to open BCF file '{bcf_path}': {e}")

    # Step 3: Match_gene_bcf
    genes = sorted(org_meta["gene"].astype(str).unique())
    gene_to_contig = match_gene_bcf(bcf, genes)

    results = []

    for _, row in org_meta.iterrows():
        species = str(row["species"])
        gene = str(row["gene"])
        position = int(row["position"])
        gt_dp = int(row["gt_DP"])

        if gene not in gene_to_contig:
            print(f"[WARN] Skipping {gene}:{position} - no contig mapping found in BCF.")
            variant = "wt"
            contig_name = "-"
        else:
            contig_name = gene_to_contig[gene]

            # Step 4: Extract_pos_variants
            variants = extract_pos_variants(
                bcf=bcf,
                contig=contig_name,
                position=position,
                snp_region_buffer=snp_region_buffer,
            )

            # Step 5: DP filtering step
            variants = filter_variants_by_dp(variants, min_dp=gt_dp)

            # Step 6: Check_pos_Deletions
            variant = check_pos_deletions(variants, position)

            # Step 7: Check_pos_SNP (only if no deletion)
            if variant is None:
                variant = check_pos_snp(variants, position)

            # If neither deletion nor SNP was found
            if variant is None:
                variant = "wt"

        results.append(
            {
                "species": species,
                "contig_name": contig_name,
                "gene": gene,
                "position": position,
                "variant": variant,
            }
        )

    out_df = pd.DataFrame(
        results,
        columns=["species", "contig_name", "gene", "position", "variant"],
    )
    out_df.to_csv(output_path, sep="\t", index=False)


def main() -> None:
    arg = argparse.ArgumentParser(
        description="Identify SNPs and deletions at configured positions from a BCF and metafile (no KMA/BED)."
    )
    arg.add_argument("--organism", required=True, help="Organism name to filter metafile by.")
    arg.add_argument("--call", required=True, help="BCF file with genotype calls.")
    arg.add_argument(
        "--metafile",
        required=True,
        help="TSV metafile with columns: species, gene, position, reference, alternative, gt_DP.",
    )
    arg.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output TSV path (species, contig_name, gene, position, variant).",
    )
    arg.add_argument(
        "--snp_region_buffer",
        type=int,
        default=1,
        help="bp to pad each expected position on each side (default 1).",
    )

    args = arg.parse_args()
    run(
        organism=args.organism,
        bcf_path=args.call,
        meta_path=args.metafile,
        output_path=args.output,
        snp_region_buffer=args.snp_region_buffer,
    )


if __name__ == "__main__":
    main()
