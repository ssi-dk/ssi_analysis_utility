#!/usr/bin/env python3
import pandas as pd
import argparse
from typing import List, Dict, Tuple
import os
import sys

# ------------------------- Metafile & Thresholds ------------------------- #

def load_thresholds_from_metafile(metafile: str, organism: str) -> Tuple[Dict[str, List[float]], List[str]]:
    """
    Read a TSV metafile with columns:
      organism  gene  coverage_min  identity_min  depth_min
    Returns:
      thresholds: dict[gene] -> [coverage_min, identity_min, depth_min] for the given organism
      gene_columns: list of genes (excluding 'other') for presence/absence columns in output
    """
    if not os.path.exists(metafile):
        raise FileNotFoundError(f"Metafile not found: {metafile}")

    meta = pd.read_csv(metafile, sep="\t", dtype={"organism": str, "gene": str})
    req = {"organism", "gene", "coverage_min", "identity_min", "depth_min"}
    if not req.issubset(meta.columns):
        raise ValueError(f"Metafile is missing required columns: {sorted(req)}")

    sub = meta[meta["organism"].astype(str) == organism]
    if sub.empty:
        raise ValueError(f"No rows for organism '{organism}' in metafile {metafile}")

    thresholds: Dict[str, List[float]] = {}
    gene_columns: List[str] = []
    for _, r in sub.iterrows():
        g = str(r["gene"]).strip()
        thresholds[g] = [float(r["coverage_min"]), float(r["identity_min"]), float(r["depth_min"])]
        if g.lower() != "other":
            gene_columns.append(g)

    # keep order as in metafile (remove dupes)
    gene_columns = list(dict.fromkeys(gene_columns))
    return thresholds, gene_columns

# ------------------------- Parsing & Threshold resolution ------------------------- #

def parse_gene_from_template(template: str) -> Tuple[str, str]:
    """
    Extract (gene, allele) from KMA '#Template'.
    Supports '__' as primary delimiter, '_' as fallback.
    """
    NA = "-"
    parts = template.split("__") if "__" in template else template.split("_")
    gene = parts[1] if len(parts) >= 2 else parts[0]
    allele = parts[2] if len(parts) >= 3 else NA
    return gene, allele

def resolve_threshold_for_gene(gene: str, thresholds: Dict[str, List[float]]) -> List[float]:
    """
    Priority:
      1) exact key match
      2) case-insensitive exact
      3) prefix match (e.g., 'stx1' uses 'stx' thresholds)
      4) 'other' (case-insensitive)
    Returns a list: [coverage_min, identity_min, depth_min]
    """
    if gene in thresholds:
        return thresholds[gene]
    for k in thresholds:
        if k.lower() == gene.lower():
            return thresholds[k]
    for k in thresholds:
        if k.lower() != "other" and gene.lower().startswith(k.lower()):
            return thresholds[k]
    for k in thresholds:
        if k.lower() == "other":
            return thresholds[k]
    # if we get here, metafile lacks both a matching entry and 'other'
    raise ValueError(f"No threshold for gene '{gene}' and no 'other' fallback in metafile.")

# ------------------------- KMA Processing (kept like your original) ------------------------- #

def get_threshold(template_name: str, thresholds: Dict[str, List[float]]) -> List[float]:
    """
    Your original signature — now parses the gene from #Template,
    resolves it via the thresholds dict, and returns a 3-item list.
    """
    gene, _allele = parse_gene_from_template(template_name)
    return resolve_threshold_for_gene(gene, thresholds)

def process_kma_res(res_file: str, thresholds: Dict[str, List[float]]) -> pd.DataFrame:
    df = pd.read_csv(res_file, sep="\t")

    required_cols = {"#Template", "Template_length", "Template_Coverage", "Template_Identity", "Depth"}
    if not required_cols.issubset(df.columns):
        raise ValueError(f"Missing required columns in {res_file}: {df.columns.tolist()}")

    # exactly like your original: attach 'threshold' column of lists, filter via indices
    df["threshold"] = df["#Template"].apply(lambda x: get_threshold(x, thresholds))
    filtered_df = df[
        (df["Template_Coverage"] >= df["threshold"].apply(lambda x: x[0])) &
        (df["Template_Identity"]   >= df["threshold"].apply(lambda x: x[1])) &
        (df["Depth"]               >= df["threshold"].apply(lambda x: x[2]))
    ].copy()

    print("Applied thresholds:")
    for _, row in df.iterrows():
        t = row["threshold"]
        print(f"\t {row['#Template']}: Coverage>={t[0]}, Identity>={t[1]}, Depth>={t[2]}")

    # keep parsed gene/allele for summarizers
    parsed = df["#Template"].apply(parse_gene_from_template)
    df["__gene"] = parsed.apply(lambda x: x[0])
    df["__allele"] = parsed.apply(lambda x: x[1])
    parsed_f = filtered_df["#Template"].apply(parse_gene_from_template)
    filtered_df["__gene"] = parsed_f.apply(lambda x: x[0])
    filtered_df["__allele"] = parsed_f.apply(lambda x: x[1])

    return filtered_df

# ------------------------- Summaries ------------------------- #

def summarize_generic(sample_id: str, organism: str, filtered_df: pd.DataFrame,
                      gene_columns: List[str], verbose_flag: int = 1) -> Dict[str, str]:
    NA = "-"
    result = {g: "Negative" for g in gene_columns}
    result.update({"Other": NA, "verbose": NA, "sample_id": sample_id, "organism": organism})

    if filtered_df.empty:
        return result

    other_genes = set()
    verbose_parts = []

    for _, row in filtered_df.iterrows():
        gene = row["__gene"]
        allele = row["__allele"]

        matched = None
        for g in gene_columns:
            if g == gene or g.lower() == gene.lower() or gene.lower().startswith(g.lower()):
                matched = g
                break
        if matched:
            result[matched] = "Positive"
        else:
            other_genes.add(gene)

        if verbose_flag:
            verbose_parts.append(
                f"{gene}_{allele}_{row['Template_length']}_{min(row['Template_Coverage'], 100.0):.2f}_{row['Template_Identity']:.2f}_{row['Depth']:.2f}"
            )

    if other_genes:
        result["Other"] = ";".join(sorted(other_genes))
    if verbose_parts:
        result["verbose"] = ";".join(verbose_parts)
    return result

def summarize_ecoli(sample_id: str, organism: str, filtered_df: pd.DataFrame,
                    gene_columns: List[str], verbose_flag: int = 1) -> Dict[str, str]:
    NA = "-"
    output = {g: "Negative" for g in gene_columns}
    output.update({"sample_id": sample_id, "organism": organism, "OH": NA, "stx": NA, "Other": NA, "verbose": NA})

    if filtered_df.empty:
        return output

    stx_alleles = set()
    other_genes = set()
    verbose_parts = []

    fli = NA
    fliC = NA
    wzx_allele = wzy_allele = wzt_allele = wzm_allele = NA

    for _, row in filtered_df.iterrows():
        gene = row["__gene"]
        allele = row["__allele"]

        for g in gene_columns:
            if g == gene or g.lower() == gene.lower() or gene.lower().startswith(g.lower()):
                output[g] = "Positive"
                break
        else:
            other_genes.add(gene)

        if gene == "wzx": wzx_allele = allele
        elif gene == "wzy": wzy_allele = allele
        elif gene == "wzt": wzt_allele = allele
        elif gene == "wzm": wzm_allele = allele
        elif gene == "fli": fli = allele
        elif gene == "fliC": fliC = allele

        if gene.lower().startswith("stx") and allele != NA:
            stx_alleles.add(allele)

        if verbose_flag:
            verbose_parts.append(
                f"{gene}_{allele}_{row['Template_length']}_{min(row['Template_Coverage'], 100.0):.2f}_{row['Template_Identity']:.2f}_{row['Depth']:.2f}"
            )

    # O-type logic
    Otype = NA

    # identify gene pairs used to determine O-type
    wzx_wzy_present = [allele for allele in (wzx_allele, wzy_allele) if allele != NA] # require both wzx and wzy to infer O type

    wzm_wzt_present = [allele for allele in (wzt_allele, wzm_allele) if allele != NA] # require both wzt and wzm to infer O type for a different pair
    
    print("pair1_present:", wzx_wzy_present)
    print("pair2_present:", wzm_wzt_present)

    wzx_wzy_set = set(wzx_wzy_present)
    wzm_wzt_set = set(wzm_wzt_present)
    print("pair1_set:", wzx_wzy_set)
    print("pair2_set:", wzm_wzt_set)

    Otype = NA

    # rare occasions where all four genes are present and identical - probably only for O8 - https://pmc.ncbi.nlm.nih.gov/articles/PMC4508402/
    if len(wzx_wzy_present) == 2 and len(wzm_wzt_present) == 2:
        union_set = set(wzx_wzy_present + wzm_wzt_present)
        print("union_set:", union_set)
        if len(union_set) == 1:
            Otype = next(iter(union_set))

    # If Otype is not set by rare example above then we have a strong call for remaining O types - if both members of the given "gene pairs" are present and identical
    if Otype == NA:
        if len(wzx_wzy_present) == 2 and len(wzx_wzy_set) == 1:
            Otype = next(iter(wzx_wzy_set))
        elif len(wzm_wzt_present) == 2 and len(wzm_wzt_set) == 1:
            Otype = next(iter(wzm_wzt_set))

    Htype = fli if fli != NA else fliC
    if Otype != NA or Htype != NA:
        output["OH"] = f"{Otype if Otype != NA else '-'};{Htype if Htype != NA else '-'}"

    if stx_alleles:
        output["stx"] = ";".join(sorted(stx_alleles))

    if other_genes:
        output["Other"] = ";".join(sorted(other_genes))
    if verbose_parts:
        output["verbose"] = ";".join(verbose_parts)
    return output

# ------------------------- Main ------------------------- #

def main(args):
    try:
        thresholds, gene_columns = load_thresholds_from_metafile(args.metafile, args.organism)

        print(f"Processing KMA .res file: {args.KMA_res}")
        filtered_df = process_kma_res(args.KMA_res, thresholds)

        if args.organism in {"Escherichia coli", "E. coli", "E.coli"}:
            result_dict = summarize_ecoli(args.sample_id, args.organism, filtered_df, gene_columns, verbose_flag=args.verbose)
        else:
            result_dict = summarize_generic(args.sample_id, args.organism, filtered_df, gene_columns, verbose_flag=args.verbose)

    except Exception as e:
        print(f"Error processing {args.sample_id}: {e}")
        # fallback on error: try to at least create all gene columns (Negative)
        try:
            thresholds, gene_columns = load_thresholds_from_metafile(args.metafile, args.organism)
            base_cols = {g: "Negative" for g in gene_columns}
        except Exception:
            base_cols = {}
        result_dict = {**base_cols, "Other": "-", "verbose": "-", "sample_id": args.sample_id, "organism": args.organism}

    output_df = pd.DataFrame([result_dict])

    # ensure all metafile genes appear as columns (Negative if missing)
    for g in gene_columns:
        if g not in output_df.columns:
            output_df[g] = "Negative"

    first_cols = ["sample_id", "organism"]
    gene_cols_in_order = [g for g in gene_columns if g in output_df.columns]
    remaining = [c for c in output_df.columns if c not in first_cols + gene_cols_in_order]
    output_df = output_df[first_cols + gene_cols_in_order + remaining]

    output_df.to_csv(args.output, sep="\t", index=False)
    print(f"Summary written to: {args.output}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Filter KMA .res using thresholds from a TSV metafile (organism + gene)."
    )
    parser.add_argument("--KMA_res", required=True, help="Input KMA .res file")
    parser.add_argument("--metafile", required=True, help="TSV with: organism,gene,coverage_min,identity_min,depth_min")
    parser.add_argument("--organism", required=True, help="Organism to match in the metafile")
    parser.add_argument("--sample_id", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--verbose", type=int, choices=[0, 1], default=1)
    args = parser.parse_args()
    main(args)

# python KMAfilter.py --KMA_res ../../examples/Results/SRR10518319/Cdiff_KMA_Toxin/SRR10518319.res --Gene_list tcdA tcdB tcdC cdtAB --organism "Clostridioides difficile" --sample_id SRR10518319 --output ../../examples/Results/SRR10518319/KMA.tsv