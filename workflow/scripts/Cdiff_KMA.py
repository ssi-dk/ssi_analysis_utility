import pandas as pd
import argparse
import os
from typing import List, Dict
import pysam
import logging
from logging_utils import setup_logging
from thresholds import get_kma_thresholds_for_species, get_threshold

# ========================= KMA FILE HANDLING =========================
def process_res_file(res_file_path: str, thresholds: Dict[str, List[int]]) -> pd.DataFrame:
    try:
        res_df = pd.read_csv(res_file_path, sep="\t")
    except FileNotFoundError:
        raise FileNotFoundError(f"File not found: {res_file_path}")
    except pd.errors.EmptyDataError:
        raise ValueError(f"File is empty or not properly formatted: {res_file_path}")

    required_columns = {"#Template", "Template_Coverage", "Query_Identity", "Depth"}
    if not required_columns.issubset(res_df.columns):
        raise ValueError(f"Missing expected columns in {res_file_path}. Found: {', '.join(res_df.columns)}")

    res_df["threshold"] = res_df["#Template"].apply(lambda x: get_threshold(x, thresholds))
    res_df_filtered = res_df[
        (res_df["Template_Coverage"] >= res_df["threshold"].apply(lambda x: x[0])) &
        (res_df["Query_Identity"] >= res_df["threshold"].apply(lambda x: x[1])) &
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
                depth = row["Depth"]
                coverage = row["Template_Coverage"]
                identity = row["Query_Identity"]
                verbose_parts.append(f"{gene}_{depth:.2f}_{coverage:.2f}_{identity:.2f}")
        output_data["verbose"] = ";".join(verbose_parts)

    logging.info(f"Successfully processed sample: {sample_name}")
    return output_data, filtered_df

# ========================= TOXIN VARIATION FILE HANDLING =========================

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
        print(f"start {start} converted start {converted_start} - end {end} converted end {converted_end}")
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

def check_deletions_in_region(
    bcf_path: str,
    contig: str,
    gene_name: str,
    coord_dict: Dict[str, Dict[str, str | int]],
    target_regions: dict[int, tuple[int, int]],  # already in contig-relative forward coordinates
    region_buffer: int = 5,
    length_tolerance: int = 1
) -> tuple[str, str]:
    """
    Scan for deletions overlapping BCF regions (with ±region_buffer) and matching length (±tolerance).
    
    Args:
        bcf_path (str): Path to BCF file.
        contig (str): Contig name in BCF (e.g. AM180355.1_tcdC_...).
        gene_name (str): Gene symbol (e.g., tcdC).
        coord_dict (dict): Gene coordinate dictionary.
        target_regions (dict): Pre-converted region coords: {length: (start, end)} in forward strand.
        region_buffer (int): Buffer size (in nt) to apply symmetrically to region boundaries.
        length_tolerance (int): Allowed deviation for matching deletion lengths.

    Returns:
        Tuple[str, str]: (matched deletion keys like "18;36", detailed string like "pos_REF_ALT_len")
    """
    print(f"Scanning deletions in {bcf_path} for {gene_name}")
    try:
        bcf = pysam.VariantFile(bcf_path)
        matched_lengths = []
        deletion_info = []

        # Pad each region explicitly by ±region_buffer
        padded_regions = {
            expected_len: (
                max(1, region_start - region_buffer),
                region_end + region_buffer
            )
            for expected_len, (region_start, region_end) in target_regions.items()
        }

        for rec in bcf.fetch(contig):
            ref = rec.ref
            alt = rec.alts[0] if rec.alts else None
            pos = rec.pos

            if alt and len(ref) > len(alt):  # deletion
                del_len = len(ref) - len(alt)
                del_start = pos
                del_end = pos + del_len - 1
                deletion_info.append(f"{pos}_{ref}_{alt}_{del_len}")

                for expected_len, (region_start, region_end) in padded_regions.items():
                    if abs(del_len - expected_len) > length_tolerance:
                        continue

                    # Compute overlap
                    overlap_start = max(del_start, region_start)
                    overlap_end = min(del_end, region_end)
                    overlap = max(0, overlap_end - overlap_start + 1)
                    target_len = region_end - region_start + 1

                    print(f"{del_len}bp deletion at {del_start}-{del_end} vs region {region_start}-{region_end} (overlap={overlap})")

                    if overlap / target_len >= 0.6:
                        print(f"Matched deletion {expected_len}bp at {del_start}-{del_end}")
                        matched_lengths.append(str(expected_len))
                        break

        return (
            ";".join(sorted(set(matched_lengths))) if matched_lengths else "-",
            ";".join(deletion_info) if deletion_info else "-"
        )

    except Exception as e:
        print(f"Error while scanning deletions: {e}")
        return "-", "-"

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
    variation_columns = ["tcdC117", "tcdCdel", "deletion_details"]
    desired_columns = input_columns + summary_columns + variation_columns

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
        res_path = f"examples/Results/{sample}/Cdiff_KMA_Toxin/{sample}.res"
        bcf_path = f"examples/Results/{sample}/GenotypeCalls/{sample}.Cdiff_KMA_Toxin.calls.bcf"

        print(f"\n================= Processing sample: {sample} =================")

        # Threshold loading
        try:
            thresholds = get_kma_thresholds_for_species(row_dict["organism"])
        except ValueError as e:
            logging.error(f"Sample '{sample}' skipped: {e}")
            continue

        # Summarize result
        result_dict, filtered_df = summarize_single_sample(sample, res_path, verbose_flag=args.verbose, thresholds=thresholds)
        row_dict.update(result_dict)

        if filtered_df.empty:
            row_dict["tcdC117"] = "-"
            row_dict["tcdCdel"] = "-"
        else:
            try:
                contig = find_matching_contig(filtered_df, "tcdC")

                # --- SNP check: tcdC position 117 (gene-relative) ---
                _, genomic_pos = gene_pos_to_genomic("tcdC", 117, coord_dict)
                print(f"📍 Position 117 in tcdC maps to contig pos: {genomic_pos}")

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
                # --- Deletion check: convert coordinates BEFORE checking ---
                gene_length = coord_dict["tcdC"]["length"]
                converted_regions = convert_reverse_strand_regions(tcdC_deletion_regions, gene_length)

                del_status, del_details = check_deletions_in_region(
                    bcf_path=bcf_path,
                    contig=contig,
                    gene_name="tcdC",
                    coord_dict=coord_dict,
                    target_regions=converted_regions,
                    region_buffer=5,
                    length_tolerance=1
                )
                row_dict["tcdCdel"] = del_status
                row_dict["deletion_details"] = del_details

            except Exception as e:
                logging.warning(f"Failed tcdC deletion check for {sample}: {e}")
                row_dict["tcdCdel"] = "-"

        # Final formatting
        for col in desired_columns:
            row_dict.setdefault(col, "-")
        result_rows.append(row_dict)

        # Write individual sample output
        if args.outputfile is None:
            output_path = os.path.join("examples", "Results", sample, "kma", f"{sample}_kma.{args.suffix}")
            os.makedirs(os.path.dirname(output_path), exist_ok=True)
            pd.DataFrame([row_dict])[desired_columns].to_csv(output_path, sep="\t" if args.suffix == "tsv" else ",", index=False)
            logging.info(f"Wrote per-sample summary: {output_path}")

    # Write merged output if requested
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




def check_variant_at_position(bcf_path: str, contig: str, pos: int) -> str:
    print("inside check_variant_at_position")
    try:
        bcf = pysam.VariantFile(bcf_path)
        print(bcf)
        for rec in bcf.fetch(contig, pos - 1, pos):
            print(f"RECORD IS {rec}")
            if rec.pos == pos:
                return "Present"
        return "-"
    except Exception as e:
        logging.warning(f"Could not read or query BCF {bcf_path}: {e}")
        return "-"
    
def find_contig_for_gene(filtered_df: pd.DataFrame, gene_name: str) -> str:
    """
    Extracts the actual contig name from a KMA template in the .res file,
    based on a known gene (e.g., 'tcdC').

    Args:
        filtered_df (pd.DataFrame): Filtered KMA result DataFrame.
        gene_name (str): Target gene to match (e.g., 'tcdC').

    Returns:
        str: Extracted contig name (e.g., 'AM180355.1')
    """
    for template in filtered_df["#Template"]:
        if gene_name in template:
            return template.split("_")[0]
    raise ValueError(f"No contig found for gene {gene_name} in templates.")
