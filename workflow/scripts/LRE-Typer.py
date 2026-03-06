#!/usr/bin/env python3

import sys
import argparse
import gzip
import pandas as pd
import os
import warnings


"""
LRE-Finder
--------------------------------------

This script is a reimplementation of the LRE-finder tool 
(https://bitbucket.org/genomicepidemiology/lre-finder/src/master/) 
It parses KMA `.mat.gz` and `.res` files, filters relevant genes,
and calculates nucleotide frequencies and mutation percentages at 
specific positions (e.g., 2505 and 2576) for Enterococcus 23S genes.

Author: Simone Scrima
Date: [2025-10-22]
"""

def read_res_file(res_path):
    """
    Parse a .res file to extract a list of valid gene names (templates).

    Parameters
    ----------
    res_path : str
        Path to the .res file produced by KMA.

    Returns
    -------
    list
        A list of gene/template names found in the "#Template" column.
    """

    try:
        df = pd.read_csv(res_path, sep="\t")
        return df["#Template"].dropna().tolist()
    except Exception as e:
        warnings.warn(f"[ERROR] Could not parse {res_path}: {e}")
        sys.exit("Quitting...")



def parse_mat_file(mat_path, 
                   wanted_genes):
    """
    Parse a KMA .mat.gz file and build a dictionary of gene matrices.

    Parameters
    ----------
    mat_path : str
        Path to the .mat.gz file.
    wanted_genes : list
        Gene names to retain.

    Returns
    -------
    dict
        {gene_name: pandas.DataFrame} with columns ["RefBase", "A", "C", "G", "T", "N", "Gap"]
    """

    # Initialize storage
    gene_dict = {}
    current_gene = None
    current_rows = []


    # Open the file (gzip handles both compressed and uncompressed)
    with gzip.open(mat_path, "rt") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                # Skip blank lines
                continue

            if line.startswith("#"):
                # Save previous gene if one is being processed
                if current_gene and current_rows:
                    df = pd.DataFrame(
                        current_rows,
                        columns=["RefBase", "A", "C", "G", "T", "N", "Gap"]
                    )
                    gene_dict[current_gene] = df

                # Start a new gene section
                current_gene = line
                current_rows = []

            else:
                # Data line: RefBase + 6 counts
                parts = line.split()
                if len(parts) == 7:
                    current_rows.append(parts)
            # Store the last gene
            df = pd.DataFrame(
                current_rows,
                columns=["RefBase", "A", "C", "G", "T", "N", "Gap"]
            )
            gene_dict[current_gene] = df

    for gene in wanted_genes:
        if gene in gene_dict:
            continue
        else:
            gene_dict.pop(gene, None) # Remove useless matrices

    return gene_dict




def count_occurences(gene_dict):
    """
    Calculate base frequencies and mutation percentages at positions 2505 and 2576.

    Parameters
    ----------
    gene_dict : dict
        Dictionary of gene matrices produced by `parse_mat_file`.

    Returns
    -------
    pandas.DataFrame
        Summary DataFrame with raw counts and per-base percentages.
    """

    e_faecium = { "G2505A" : 0.10,
                  "G2576T" : 0.10}
    e_faecalis = { "G2505A" : 0.10,
                   "G2576T" : 0.10}

    results = []
    cols = ["A", "C", "G", "T", "N", "Gap"]

    for gene, matrix in gene_dict.items():
        if gene not in {"#23S_Enterococcus_faecium", "#23S_Enterococcus_faecalis"}:
            continue

        # Remove gap rows and reset indexing
        matrix = matrix[matrix["RefBase"] != "-"].reset_index(drop=True)

        # Check if target positions exist
        if len(matrix) <= 2575:
            print(f"[WARNING] {gene} too short — skipping.")
            continue

        # Ensure both positions are G
        if matrix.loc[2504, "RefBase"] != "G" or matrix.loc[2575, "RefBase"] != "G":
            print(f"[INFO] {gene} has non-G reference bases — skipping.")
            continue

        # Extract the two target positions
        subset = matrix.loc[[2504, 2575]].copy()
        subset[cols] = subset[cols].apply(pd.to_numeric, errors="coerce")

        # Compute per-base percentages (relative to total coverage)
        subset_perc = subset[cols].div(subset[cols].sum(axis=1), axis=0) * 100
        subset_perc = subset_perc.round(2)
        subset_perc.columns = ["A[%]", "C[%]", "G[%]", "T[%]", "N[%]", "-[%]"]

        # Merge raw + percentage data
        gene_final = pd.concat([subset, subset_perc], axis=1)
        gene_final["POS"] = gene_final.index + 1
        gene_final.index = [gene] * len(gene_final)
        results.append(gene_final)

    if not results:
        warnings.warn("No valid genes processed — all skipped.")
        return pd.DataFrame()

    df_final = pd.concat(results)
    return df_final
       




def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-ires',
                        '--res_file',
                        dest='res',
                        type=str,
                        required=True,
                        metavar='',
                        help='.res file output from KMA'
                        )

    parser.add_argument('-imat',
                        '--matrix',
                        dest='mat',
                        type=str,
                        required=True,
                        metavar='',
                        help='.mat.gz file output from KMA'
                        )
    
    parser.add_argument('-o',
                        '--pos_file',
                        dest='pos',
                        type=str,
                        required=True,
                        metavar='',
                        help='Position description file (tab-delimited)'
                        )

    args = parser.parse_args()

    # Define Flags
    res_file = os.path.abspath(args.res)
    mat = os.path.abspath(args.mat)
    pos = args.pos



    # Read .res file for valid genes
    wanted_genes = read_res_file(res_file)
    
    if not wanted_genes:
        sys.exit(f"Error: Could not read valid templates from {res_file}")
 
    d = parse_mat_file(mat, 
                       wanted_genes)
    df_final = count_occurences(d)
    df_final = df_final[["POS",
                         "RefBase", 
                         "A",
                         "C",
                         "G",
                         "T",
                         "N",
                         "Gap",
                         "A[%]",
                         "C[%]",
                         "G[%]",
                         "T[%]",
                         "N[%]",
                         "-[%]"]]
    base, ext = os.path.splitext(pos)
    pos = pos if pos.endswith('.tsv') else pos + '.tsv'
    
    df_final.to_csv(pos,
                    sep="\t") # write final tsv 


if __name__ == "__main__":
    main()

