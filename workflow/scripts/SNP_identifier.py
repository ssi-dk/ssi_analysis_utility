#!/usr/bin/env python3

from __future__ import annotations
import argparse
import os
import sys
from typing import Dict, List, Tuple

import pandas as pd
import pysam


# ========================= Helpers =========================

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

def read_meta(meta_path: str, organism: str) -> pd.DataFrame:
    """Read and filter SNP meta TSV for the requested organism.

    Expected columns: species, gene, position, reference, alternative, strand
    """
    df = pd.read_csv(meta_path, sep="\t")
    required_col = {"species", "gene", "position", "reference", "alternative", "strand"}
    
    missing_col = required_col - set(df.columns)
    if missing_col:
        raise ValueError(f"provided snp information file to (--meta) is missing required columns: {', '.join(sorted(missing_col))}")
    
    # Normalize column types
    df["gene"] = df["gene"].astype(str)
    df["species"] = df["species"].astype(str)
    df["position"] = df["position"].astype(int)
    df["reference"] = df["reference"].astype(str)
    df["alternative"] = df["alternative"].astype(str)
    df["strand"] = df["strand"].astype(str)

    return df[df["species"] == organism].copy()

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
        raise ValueError(f"No contig matching gene '{gene}' in .res")
    return str(hits.iloc[0]["#Template"])  # first match

def gene_pos_to_contig_pos(gene: str, pos_in_gene: int, bed_info: Dict[str, Dict[str, object]]) -> int:
    if gene not in bed_info:
        raise ValueError(f"Gene '{gene}' not present in BED6")
    info = bed_info[gene]
    strand = info["strand"]
    length = int(info["length"])  # gene length
    if strand == "+":
        return pos_in_gene
    elif strand == "-":
        return length - (pos_in_gene - 1)
    else:
        raise ValueError(f"Invalid strand '{strand}' for gene '{gene}'")

def check_snp_variant(
    bcf_path: str,
    contig: str,
    pos: int,
    expected_ref: str,
    expected_alt: str,
    strand: str,
    flank: int = 20,
) -> Tuple[str, str]:
    """Return (label, detail) where label in {snp, wt, del, other, -}.

    - If a variant at exactly `pos` matches expected_ref>expected_alt → "snp".
    - If a different SNP at exactly `pos` → ("other", "pos_ref_alt").
    - If a deletion spans `pos` → "del".
    - If nothing in the window → "wt".
    - On error → ("-", "").
    """
    try:
        bcf = pysam.VariantFile(bcf_path)
        start = max(0, pos - flank)
        end = pos + flank
        
        for rec in bcf.fetch(contig, start, end):

            ref = rec.ref
            alt = rec.alts[0] if rec.alts else None
            
            ref_to_use = ref
            alt_to_use = alt

            if strand == "-":
                ref_to_use = reverse_complement(ref)
                alt_to_use = reverse_complement(alt)

            if rec.pos == pos:
                # Case 1: variation fitting the A>T
                if ref == expected_ref and alt == expected_alt:
                    print(f"\tCase 1: {ref}>{expected_ref} and {alt}>{expected_alt} SNP present at 117")
                    return f"snp", ""
            
                # Case 2: different variation
                else:
                    print("\tCase 2: Different SNP present at 117")
                    return "other", f"{rec.pos}_{ref_to_use}_{alt_to_use}"
            
            # Case 3: Deletion spanning the position
            elif rec.pos < pos:
                deletion_end = rec.pos + len(ref) - 1
                if deletion_end >= pos and any(len(a) < len(ref) for a in rec.alts if a is not None):
                    print(f"\tCase 3: A deletion spans the position at {pos} from {rec.pos} to {deletion_end}")
                    return f"del", ""

        return ("wt", "")
    except Exception as e:
        print(f"[ERROR] BCF read failed at {contig}:{pos}: {e}")
        return ("-", "")


# ========================= Main routine =========================

def run(organism: str,
        res_path: str,
        bcf_path: str,
        bed_path: str,
        meta_path: str,
        output_path: str) -> None:

    # Inputs
    res_df = process_res_file(res_path)
    bed_info = load_bed6(bed_path)
    meta_df = read_meta(meta_path, organism)

    if meta_df.empty:
        # No SNPs defined for this organism - create empty file '-'
        out = pd.DataFrame([{"organism": organism, "snp_info": "-"}])
        out.to_csv(output_path, sep="\t", index=False)
        return None

    # Group by gene
    snp_calls: List[str] = []
    contig_cache: Dict[str, str] = {}

    for idx, row in meta_df.iterrows():
        gene = str(row["gene"])  # e.g., tcdC
        pos = int(row["position"])  # 1-based in gene
        ref = str(row["reference"]).upper()
        alt = str(row["alternative"]).upper()
        strand = str(row["strand"]).strip()

        # Locate contig for this gene (from .res)
        try:
            if gene not in contig_cache:
                contig_cache[gene] = find_contig_for_gene(res_df, gene)
            contig = contig_cache[gene]
        except Exception as e:
            print(f"[WARN] {gene}: contig not found in .res → marking as '-' ({e})")
            snp_calls.append(f"{pos}_-")
            continue

        # Map gene pos → contig pos (gene-sized contigs)
        try:
            contig_pos = gene_pos_to_contig_pos(gene, pos, bed_info)
        except Exception as e:
            print(f"[WARN] {gene}: coordinate transform failed → '-' ({e})")
            snp_calls.append(f"{pos}_-")
            continue

        # For '-' strand, compare against reverse-complement expectations
        ref_cmp, alt_cmp = (ref, alt) if strand == "+" else (reverse_complement(ref), reverse_complement(alt))

        label, detail = check_snp_variant(
            bcf_path=bcf_path,
            contig=contig,
            pos=contig_pos,
            expected_ref=ref_cmp,
            expected_alt=alt_cmp,
            strand=strand,
        )

        if label == "snp":
            snp_calls.append(f"{pos}_{ref}>{alt}")
        elif label == "wt":
            snp_calls.append(f"{pos}_wt")
        elif label == "del":
            snp_calls.append(f"{pos}_del")
        elif label == "other":
            snp_calls.append(f"{pos}_other:{detail}")
        else:
            snp_calls.append(f"-")

    row = {
        "organism": organism,
        "snp_info": ";".join(snp_calls) if snp_calls else "-",
    }

    pd.DataFrame([row]).to_csv(output_path, sep="\t", index=False)


def main():
    ap = argparse.ArgumentParser(description="Identify configured SNPs from BCF using meta TSV.")
    ap.add_argument("--organism", required=True)
    ap.add_argument("--res", required=True, help="KMA .res file for contig lookup")
    ap.add_argument("--call", required=True, help="Main genotype BCF")
    ap.add_argument("--bed", required=True, help="BED6 with gene coordinates/strand")
    ap.add_argument("--metafile", required=True, help="TSV: species gene position reference alternative strand")
    ap.add_argument("-o", "--output", required=True, help="Output TSV path")
    args = ap.parse_args()

    run(
        organism=args.organism,
        res_path=args.res,
        bcf_path=args.call,
        bed_path=args.bed,
        meta_path=args.metafile,
        output_path=args.output,
    )


if __name__ == "__main__":
    main()

"""
SNP_identifier.py

Extracts SNP status per sample from a meta TSV describing expected SNPs.

Input files (per sample):
  --res     : KMA .res file (tab-separated) used to locate contig(s) per gene
  --call    : main genotype BCF (SNPs/indels)
  --bed     : BED6 for genes (contig, start, end, gene, score, strand)
  --meta    : TSV with columns: species, gene, position, reference, alternative, strand

CLI example:
  python SNP_identifier.py \
    --organism "Clostridioides difficile" \
    --res examples/Results/SRR2915763/Cdiff_KMA_Toxin/SRR2915763.res \
    --call examples/Results/SRR2915763/GenotypeCalls/SRR2915763.Cdiff_KMA_Toxin.calls.bcf \
    --bed refs/toxin_genes.bed6 \
    --meta configs/snp_meta.tsv \
    -o snp_summary.tsv

Output (TSV):
  organism\tsnp_info
  SRR2915763\tClostridioides difficile\t117_del;184_wt

Notes:
- Assumes per-gene contigs in the BCF match names present in the .res file (#Template contains gene name).
- Coordinates are gene-relative on those contigs: for plus strand, pos stays pos; for minus, pos' = gene_len - (pos-1).
- If a deletion spans the SNP position, the label is "<pos>_del" (1bp deletions treated as SNP-like events).
- If a different SNP is present at the exact position, label "<pos>_other:<pos_ref_alt>"; if none, "<pos>_wt".
"""