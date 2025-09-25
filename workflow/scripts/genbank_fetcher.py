import argparse
import sys
from typing import List, Tuple, Dict, Optional, Any
from datetime import datetime
import os
import pandas as pd
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
from Bio.Seq import Seq 
from collections import defaultdict  # add this

# ---------- OUTPUT FILES ----------
# note these files are not necessary as output - as such no raise error ro similar is needed if the paths are None

def genbank_records_report(path: Optional[str], info_lines: List[str], append: bool) -> None:
    """Write the genbank information report."""
    if not path:
        return
    mode = 'a' if append else 'w'
    with open(path, mode) as f:
        f.write("\n".join(info_lines) + "\n")


def create_locus_bed(path: Optional[str], chrom: str,
                 regions: List[Tuple[int, int, str, str]],
                 append: bool) -> None:
    """
    Extract the regional information for desired locus to create a .bed6 formation.
    """
    if not path or not regions:
        return

    # create bed header
    header = "# chrom\tchromStart\tchromEnd\tname\tscore\tstrand\n"
    mode = 'a' if append else 'w'
    need_header = True # i need a boolean value here - for appended regions each region will have their own header so only write it first time
    if append and os.path.exists(path) and os.path.getsize(path) > 0:
        need_header = False
    
    # extract inferred bed file information and store values
    bed_lines = [f"{chrom}\t{chromStart}\t{chromEnd}\t{name}\t0\t{strand}"
                 for (chromStart, chromEnd, strand, name) in regions]

    with open(path, mode) as bed6file:
        if need_header:
            bed6file.write(header)
        bed6file.write("\n".join(bed_lines) + "\n")

def create_locus_fasta(path: Optional[str], record,  # Biopython SeqRecord
                       regions: List[Tuple[int, int, str, str]],
                       merge_distance: Optional[int],
                       append: bool,
                       strand_correct: bool) -> None:
    """
    Build and write FASTA entries for one GenBank record from a list of genomic regions.

    path : str | None
        Output FASTA filepath. If None or empty, this function returns immediately.
    record : Bio.SeqRecord.SeqRecord
        Biopython SeqRecord parsed from GenBank (contains .id and .seq).
    regions : list[tuple[int, int, str, str]]
        A list of regions to export. Each item is:
            (start, end, strand, name)
        where:
            - start, end are genomic, 0-based, half-open coordinates [start, end).
            - strand is '+', '-' or '.' (unknown).
            - name is the gene name for the corresponding region (e.g. gene name).
    merge_distance : int | None
        If set, adjacent regions on the SAME strand are merged when the gap between
        them is <= merge_distance. If None/0, only directly touching regions merge.
    append : bool
        If True, append to the FASTA file; otherwise overwrite.
    """
    # do nothing for no paths or regions
    if not path or not regions:
        return

    # Sort by regions by strand, start
    locus_regions = sorted(regions, key=lambda x: (x[2], x[0]))
    print(f"locus_regions : {locus_regions}")

    # Merge adjacent regions per strand when gap <= merge_distance
    merged_groups: List[List[Tuple[int, int, str, str]]] = [] # e.g. (795842, 797843, '+', 'tcdA_795843-797843')
    current = [locus_regions[0]]

    for reg_info in locus_regions[1:]:
        print(f"reg_info: {reg_info} for current {current} ")
        prev = current[-1]
        same_strand = (reg_info[2] == prev[2])

        # gap between previous end and current start
        vicinity = ((reg_info[0] - prev[1]) <= (merge_distance or 0))

        if same_strand and vicinity:
            current.append(reg_info)
        else:
            merged_groups.append(current)
            current = [reg_info]

    merged_groups.append(current)
    print(f"merged_groups: {merged_groups}")

    # Build FASTA entries (genomic orientation; NOT strand-corrected)
    fasta_entries = []
    for group in merged_groups:
        group_start = group[0][0]
        group_end = group[-1][1]
        group_strand = group[0][2]  # all items in a group share strand
        group_name = "+".join(g[3] for g in group)

        print(f"name: {group_name} \t start {group_start} \t end {group_end}")

        seq = record.seq[group_start:group_end]

        #reverse-complement or not?
        rc_applied = False
        if strand_correct and group_strand == "-":
            print(f"corrects strand orientation for record {group_name} belonging to strand {group_strand}")
            # convert to gene/CDS 5'->3' orientation
            seq = seq.reverse_complement()
            rc_applied = True

        header = f">{record.id}_{group_name}_{group_start}_{group_end}_strand_{group_strand}"
        if rc_applied:
            header += "_rc"
        
        fasta_entries.append(f"{header}\n{str(seq)}")
        
    # write fasta file either as append or overwrite
    mode = 'a' if append else 'w'
    with open(path, mode) as fasta_out:
        fasta_out.write("\n".join(fasta_entries) + "\n")

def process_metafile(
    meta_path: str,
    email: str,
    record_file: Optional[str],
    bed_file: Optional[str],
    fasta_file: Optional[str],
    merge_distance: Optional[int],
    append: bool,
    strand_correct: bool,
) -> None:
    rows = read_metafile_rows(meta_path)
    if not rows:
        print(f"Metafile '{meta_path}' is empty or unreadable.")
        return

    # group so we fetch each accession once
    by_acc: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
    for r in rows:
        acc = r.get("accession") or ""
        if not acc:
            print("Warning: row without accession; skipping.")
            continue
        by_acc[acc].append(r)

    Entrez.email = email

    for accession, group in by_acc.items():
        print(f"Fetching GenBank record for {accession} from NCBI (metafile)...")
        try:
            with Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text") as handle:
                record = SeqIO.read(handle, "genbank")
        except Exception as e:
            print(f"Failed to fetch or parse {accession}: {e}")
            continue

        ts = datetime.now().strftime("%Y-%m-%d_%H:%M:%S")
        info_lines: List[str] = [
            f"Timestamp: {ts}\n",
            f"Running command: {' '.join(sys.argv)}\n",
            f"Record ID: {record.id}",
            f"Description: {record.description}",
            f"Sequence Length: {len(record.seq)} bp",
            "Mode: --metafile",
            ""
        ]

        # index CDS once
        cds_list, key_map = collect_cds_index(record)

        # organism text (species + lineage) for tolerant matching
        organism_text = (
            (record.annotations.get("organism", "") + " " +
             " ".join(record.annotations.get("taxonomy", [])))
            .lower()
        )

        regions: List[Tuple[int, int, str, str]] = []

        for r in group:
            gene = (r.get("gene") or "").strip()
            org  = (r.get("organism") or "").strip()
            s: Optional[int] = r.get("start")
            e: Optional[int] = r.get("end")

            # 1) organism check if provided
            if org:
                if organism_text and (org.lower() not in organism_text):
                    info_lines.append(f" - Organism mismatch for row ({gene or 'NA'}): '{org}' not found in record organism/lineage; skipping.")
                    continue

            # 2) gene must exist
            if not gene:
                info_lines.append(" - Missing 'gene' in metafile row; skipping.")
                continue

            feature = match_feature_by_name(gene, key_map, cds_list)
            if not feature:
                info_lines.append(f" - '{gene}' not matched to any CDS; skipping.")
                continue

            start = int(feature.location.start)
            end   = int(feature.location.end)
            strand = strand_char(feature)
            cds_len = end - start
            base = CDS_labels(feature)

            # 3) decide which mode to emulate
            if s is None or e is None:
                # full CDS (locus)
                info_lines.append(f" - {base} (metafile locus): {start+1}..{end} ({strand}) length={cds_len}")
                regions.append((start, end, strand, base))
                continue

            # both s and e present:
            if max(s, e) <= cds_len:
                # treat as cds_region (CDS-relative 1-based)
                cds_start = max(1, min(int(s), int(cds_len)))
                cds_end   = max(1, min(int(e), int(cds_len)))
                if (cds_start, cds_end) != (s, e):
                    info_lines.append(f"\tWarning: {gene} CDS subrange {s}-{e} clamped to {cds_start}-{cds_end} (CDS length {cds_len}).")

                if strand == "+":
                    genomic_start = start + (cds_start - 1)
                    g_end         = start + cds_end
                elif strand == "-":
                    genomic_start = end - cds_end
                    g_end         = end - (cds_start - 1)
                else:
                    genomic_start = start + (cds_start - 1)
                    g_end         = start + cds_end

                label = f"{base}_{cds_start}-{cds_end}"
                info_lines.append(f" - {base} (metafile cds_region): CDS {cds_start}-{cds_end} -> genomic {genomic_start+1}..{g_end} ({strand})")
                regions.append((genomic_start, g_end, strand, label))

            else:
                # treat as acc_region (genomic 1-based inclusive)
                record_len = len(record.seq)
                acc_start_sub = max(0, min((s - 1), record_len))
                acc_end_sub   = max(0, min(e, record_len))

                # clamp to CDS bounds (like your acc_region)
                genomic_start = max(start, min(acc_start_sub, end))
                g_end = max(start, min(acc_end_sub, end))

                if g_end <= genomic_start:
                    info_lines.append(f"   Warning: empty slice for {gene}; skipping.")
                    print(f"warning empty regions for {gene}")
                    continue

                label = f"{base}_{genomic_start+1}-{g_end}"
                info_lines.append(f" - {base} (metafile acc_region): genomic {genomic_start+1}..{g_end} ({strand})")
                regions.append((genomic_start, g_end, strand, label))

        # finalize this accession
        unknown_count = sum(1 for r in regions if r[2] == ".")
        if unknown_count:
            info_lines.append(f"Warning: {unknown_count} region(s) have unknown strand '.'; they will be merged together (same '.') and never reverse-complemented.")
        genbank_records_report(record_file, info_lines, append)
        create_locus_bed(bed_file, record.id, regions, append)
        create_locus_fasta(fasta_file, record, regions, merge_distance, append, strand_correct)

# ---------- HELPERS ----------

def parse_range(text: str) -> Tuple[int, int]:
    """
    Parse a 1-based inclusive range 'START-END' into a (start, end) tuple of ints.
    
    return (start, end) : tuple[int, int]
    """
    try:
        start, end = text.replace(',', '').split('-')
        start, end = int(start), int(end)
        if start <= 0 or end <= 0 or end < start:
            raise ValueError
        return start, end
    except Exception:
        # ensure the argument input is an valid range 
        raise argparse.ArgumentTypeError(f"Invalid range '{text}'. Use START-END with positive integers.")

def parse_name_ranges_list(items: Optional[List[str]]) -> List[Tuple[str, Tuple[int, int]]]:
    """
    parse a list of strings (chrom:start-end) to accomodate multiple regions specificed in the arguments 
    ['tcdA:1-200', 'tcdA:400-600', 'tcdB:10-50'] converted to [('tcda', (1,200)), ('tcda', (400,600)), ('tcdb', (10,50))]
    """
    converted_coord: List[Tuple[str, Tuple[int, int]]] = []
    for item in items or []:
        try:
            name, region = item.split(':', 1)
            converted_coord.append((name.strip().lower(), parse_range(region)))
        except Exception:
            raise argparse.ArgumentTypeError(f"Invalid item '{item}'. Use NAME:START-END")
    return converted_coord

def collect_cds_index(record: SeqRecord) -> Tuple[List[SeqFeature], Dict[str, SeqFeature]]:
    """
    Build an index of CDS features for a GenBank record 
    """
    identifier_to_feature: Dict[str, SeqFeature] = {}
    cds_features: List[SeqFeature] = []

    for feature in record.features:
        if feature.type != "CDS":
            continue

        cds_features.append(feature)

        # The feature identifiers for genes to search for.
        for qualifier_name in ("gene", "locus_tag", "old_locus_tag"):
            for identifier in feature.qualifiers.get(qualifier_name, []):
                key = identifier.strip().lower()
                if key and key not in identifier_to_feature:
                    identifier_to_feature[key] = feature

    return cds_features, identifier_to_feature

def match_feature_by_name(
    query_name: str,
    identifier_to_feature: Dict[str, SeqFeature],
    cds_features: List[SeqFeature],
) -> Optional[SeqFeature]:
    """
    Match gene/locus/old_locus label to CDS features
    """
    
    # query name (e.g., "tcdA", "CD0660")
    key = query_name.lower()

    # identify match by gene / locus_tag / old_locus_tag
    if key in identifier_to_feature:
        return identifier_to_feature[key]

    # Fallback: product substring match (convenience; may be ambiguous)
    for feature in cds_features:
        product_text = " ".join(feature.qualifiers.get("product", []))
        if product_text and key in product_text.lower():
            return feature

    return None

def strand_char(feature) -> str:
    if feature.location.strand == 1:
        return "+"
    if feature.location.strand == -1:
        return "-"
    return "."

def CDS_labels(feature) -> str:
    """ extract the different GenBank record labels for the CDS """

    gene = feature.qualifiers.get("gene", [""])[0].strip()
    old_tag = feature.qualifiers.get("old_locus_tag", [""])[0]
    locus_tag = feature.qualifiers.get("locus_tag", [""])[0]

    return gene or old_tag or locus_tag or "unknown"

def normalize_header(name: str) -> str:
    return name.strip().lower().replace(" ", "_")

def read_metafile_rows(path: str) -> List[Dict[str, Any]]:
    """
    Read a TSV meta file with headers:
      accession  gene  start  end  organism
    Returns a list of dicts with normalized keys; start/end as Optional[int].
    """
    # strict TSV; if you want auto-detect, use sep=None, engine="python"
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)

    # normalize headers
    df.columns = [normalize_header(c) for c in df.columns]

    required = {"accession", "gene", "start", "end", "organism"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Metafile missing required columns: {', '.join(sorted(missing))}")

    # trim whitespace in text columns
    for c in ["accession", "gene", "organism", "start", "end"]:
        df[c] = df[c].astype(str).str.strip()

    # convert start/end to Optional[int]
    def to_opt_int(s: str):
        s = (s or "").replace(",", "").strip()
        return int(s) if s.isdigit() else None

    df["start"] = df["start"].map(to_opt_int)
    df["end"]   = df["end"].map(to_opt_int)

    # return list of dicts (start/end are Python ints or None)
    return df.to_dict(orient="records")

# ---------- Core fetch genbank features ----------

def fetch_and_store_genbank_features(
    email: str,
    accession: str,
    record_file: Optional[str],
    bed_file: Optional[str],
    fasta_file: Optional[str],
    merge_distance: Optional[int],
    append: bool,
    mode: str,                              # 'locus' | 'cds_region' | 'acc_region'
    locus_list: Optional[List[str]] = None,
    cds_regions: Optional[List[Tuple[str, Tuple[int,int]]]] = None,
    acc_regions: Optional[List[Tuple[str, Tuple[int,int]]]] = None,
    strand_correct: bool = True,
) -> None:
    Entrez.email = email
    print(f"Fetching GenBank record for {accession} from NCBI...")

    try:
        with Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text") as handle:
            record = SeqIO.read(handle, "genbank")
    except Exception as e:
        print(f"Failed to fetch or parse {accession}: {e}")
        return

    # 1) record_info
    record_info: List[str] = []
    ts = datetime.now().strftime("%Y-%m-%d_%H:%M:%S")
    record_info.append(f"Timestamp: {ts}\n")
    record_info.append(f"Running command: {' '.join(sys.argv)}\n")
    record_info.append(f"Record ID: {record.id}")
    record_info.append(f"Description: {record.description}")
    record_info.append(f"Sequence Length: {len(record.seq)} bp")
    record_info.append("")

    # 2) regions for BED/FASTA (0-based)
    regions: List[Tuple[int, int, str, str]] = []

    # Build CDS index
    cds_list, key_map = collect_cds_index(record)

    if mode == "locus":
        record_info.append("Mode: --locus (full CDS per locus)")
        requested = [x.strip() for x in (locus_list or [])]
        for req in requested:
            feature = match_feature_by_name(req, key_map, cds_list)
            if not feature:
                record_info.append(f" - '{req}' not matched to any CDS.")
                continue
            start = int(feature.location.start)
            end = int(feature.location.end)
            label = CDS_labels(feature)
            strand = strand_char(feature)
            record_info.append(f" - {label} (matched '{req}'): {start+1}..{end} ({strand}) length={end-start}")
            regions.append((start, end, strand, label))

    elif mode == "cds_region":
        record_info.append("Mode: --cds_region (CDS-relative sub-slices)")
        for req_name, (sub_start, sub_end) in (cds_regions or []):
            feature = match_feature_by_name(req_name, key_map, cds_list)
            if not feature:
                record_info.append(f" - '{req_name}' not matched to any CDS.")
                continue
            start = int(feature.location.start)
            end = int(feature.location.end)
            strand = strand_char(feature)
            cds_len = end - start

            # CDS-relative 1-based range to CDS length
            cds_start = max(1, min(int(sub_start), int(cds_len)))
            cds_end = max(1, min(int(sub_end), int(cds_len)))
            if (cds_start, cds_end) != (sub_start, sub_end):
                record_info.append(f"\tWarning: {req_name} CDS subrange {sub_start}-{sub_end} clamped to {cds_start}-{cds_end} (CDS length {cds_len}).")

            # Map CDS 1-based inclusive -> genomic 0-based half-open
            if strand == "+":
                genomic_start = start + (cds_start - 1)
                g_end = start + cds_end
            elif strand == "-":
                genomic_start = end - cds_end
                g_end   = end - (cds_start - 1)
            else:
                genomic_start = start + (cds_start - 1)
                g_end   = start + cds_end

            base = CDS_labels(feature)
            label = f"{base}_{cds_start}-{cds_end}"
            record_info.append(f" - {base} (matched '{req_name}'): CDS {cds_start}-{cds_end} -> genomic {genomic_start+1}..{g_end} ({strand})")
            regions.append((genomic_start, g_end, strand, label))

    elif mode == "acc_region":
        record_info.append("Mode: --acc_region (genomic/Accession sub-slices within CDS)")
        record_len = len(record.seq)

        for req_name, (acc_start, acc_end) in (acc_regions or []):
            feature = match_feature_by_name(req_name, key_map, cds_list)
            if not feature:
                record_info.append(f" - '{req_name}' not matched to any CDS.")
                continue

            if acc_start < 1 or acc_end > record_len:
                record_info.append(f"Warning: {req_name} accession range {acc_start}-{acc_end} lies outside record bounds 1-{record_len}; it will be clamped.")

            start = int(feature.location.start)  # CDS genomic bounds (0-based)
            end = int(feature.location.end)

            # Convert provided 1-based inclusive -> 0-based half-open
            acc_start_sub = int(acc_start) - 1
            acc_end_sub   = int(acc_end)

            # Clamp to CDS bounds
            genomic_start = max(start, min(acc_start_sub, end))
            g_end   = max(start, min(acc_end_sub,end))
            if (genomic_start, g_end) != (acc_start_sub, acc_end_sub):
                record_info.append(
                    f"   Warning: {req_name} accession range {acc_start}-{acc_end} clamped to {genomic_start+1}-{g_end} "
                    f"within CDS {start+1}-{end}."
                )
            if g_end <= genomic_start:
                record_info.append(f"   Warning: empty slice for {req_name} after clamping; skipping.")
                continue

            base = CDS_labels(feature)
            strand = strand_char(feature)
            label = f"{base}_{genomic_start+1}-{g_end}"
            record_info.append(f" - {base} (matched '{req_name}'): genomic {genomic_start+1}..{g_end} ({strand})")
            regions.append((genomic_start, g_end, strand, label))

    else:
        record_info.append("Error: no valid mode selected.")
        genbank_records_report(record_file, record_info, append)
        return

    # Warn if any unknown-strand regions present (merged together; never reverse-complemented)
    unknown_count = sum(1 for r in regions if r[2] == ".")
    if unknown_count:
        record_info.append(f"Warning: {unknown_count} region(s) have unknown strand '.'; they will be merged together (same '.') and never reverse-complemented.")

    # 3) finalize outputs
    genbank_records_report(record_file, record_info, append)
    create_locus_bed(bed_file, record.id, regions, append)
    create_locus_fasta(fasta_file, record, regions, merge_distance, append,strand_correct)


# ---------- CLI ----------

def main():
    parser = argparse.ArgumentParser(description="Fetch GenBank CDS features and export full or partial regions.")

    parser.add_argument("-e", "--email", default="Your.Name.Here@example.org",
                        help="Your email address (required by NCBI).")
    parser.add_argument("-a", "--accession", nargs="+",
                        help="One or more GenBank accession numbers.")

    # Exactly one of these three must be used
    parser.add_argument("-l", "--locus", nargs="*",
                        help="List of gene/locus names (e.g., tcdA tcdB) → export FULL CDS for each.")
    parser.add_argument("--cds_region", nargs="*", metavar="NAME:START-END",
                        help="CDS-relative slices (1-based, inclusive), e.g., tcdA:1-200 tcdB:300-400.")
    parser.add_argument("--acc_region", nargs="*", metavar="NAME:START-END",
                        help="Accession/genomic slices (1-based, inclusive) within each CDS, "
                             "e.g., tcdA:795843-797843.")

    parser.add_argument("--metafile", help="TSV with columns: accession, gene, start, end, organism.")

    # outputs
    parser.add_argument("-r", "--records", default="genbank_record.txt",
                        help="Text report file (default: genbank_record.txt)")
    parser.add_argument("--bed", default="genbank_coord.bed6",
                        help="BED6 output file (default: genbank_coord.bed6)")
    parser.add_argument("--fasta", default="genbank_seq.fasta",
                        help="FASTA output file with extracted slices (default: genbank_seq.fasta)")

    parser.add_argument("--merge", type=int,
                        help="Merge FASTA regions if gap ≤ this value (same strand).")
    parser.add_argument("--append", action="store_true",
                        help="Append to outputs instead of overwriting.")

    # Strand-correction toggle: default ON; disable with -ns / --no_strand_correct
    parser.add_argument(
        "-ns", "--no_strand_correct",
        action="store_true",
        help="Do NOT reverse-complement '-' strand slices; keep genomic/reference orientation."
    )

    args = parser.parse_args()

    # metafile exclusivity
    if args.metafile:
        if args.accession or args.locus or args.cds_region or args.acc_region:
            parser.error("--metafile cannot be combined with --accession/--locus/--cds_region/--acc_region.")

        # default behavior = strand correction ON
        strand_correct = not args.no_strand_correct

        # Run metafile path and exit
        process_metafile(
            meta_path=args.metafile,
            email=args.email,
            record_file=args.records,
            bed_file=args.bed,
            fasta_file=args.fasta,
            merge_distance=args.merge,
            append=args.append,
            strand_correct=strand_correct,
        )
        return

    # -------- normal (non-metafile) path below ----------

    # Validate required inputs for the 3 classic modes
    if not args.accession:
        parser.error("--accession is required unless you use --metafile.")

    # Warn if multiple accessions without append
    if args.accession and len(args.accession) > 1 and not args.append:
        print("Warning: multiple accessions provided without --append; outputs will be overwritten for each accession.", file=sys.stderr)

    # Exactly one of the three modes must be selected
    chosen = sum([bool(args.locus), bool(args.cds_region), bool(args.acc_region)])
    if chosen != 1:
        parser.error("You must specify exactly ONE of: --locus, --cds_region, or --acc_region (unless using --metafile).")

    if args.locus:
        mode = "locus"
        locus_list = args.locus
        cds_regions = None
        acc_regions = None
    elif args.cds_region:
        mode = "cds_region"
        locus_list = None
        cds_regions = parse_name_ranges_list(args.cds_region)
        acc_regions = None
    else:
        mode = "acc_region"
        locus_list = None
        cds_regions = None
        acc_regions = parse_name_ranges_list(args.acc_region)

    # default behavior = strand correction ON
    strand_correct = not args.no_strand_correct

    for acc in args.accession:
        fetch_and_store_genbank_features(
            email=args.email,
            accession=acc,
            record_file=args.records,
            bed_file=args.bed,
            fasta_file=args.fasta,
            merge_distance=args.merge,
            append=args.append,
            mode=mode,
            locus_list=locus_list,
            cds_regions=cds_regions,
            acc_regions=acc_regions,
            strand_correct=strand_correct
        )

if __name__ == "__main__":
    main()

#python scripts/genbank_fetcher.py --accession AM180355.1 --cds_region tcdA:100-200 tcdB:300-400 tcdC:20-175 --bed genbank_coord_cds.bed6 --records genbank_records_cds.txt --fasta genbank_seq_cds.fasta --merge 200
#python scripts/genbank_fetcher.py --meta meta_cds.tsv --bed genbank_coord_cds_meta.bed6 --records genbank_records_cds_meta.txt --fasta genbank_seq_cds_meta.fasta --merge 200

#python scripts/genbank_fetcher.py --accession AM180355.1 --acc_region tcdA:795942-803175 tcdB:787492-794393 tcdC:804409-804908 --bed genbank_coord_acc.bed6 --records genbank_records_acc.txt --fasta genbank_seq_acc.fasta --merge 200
#python scripts/genbank_fetcher.py --meta meta_acc.tsv --bed genbank_coord_acc_meta.bed6 --records genbank_records_acc_meta.txt --fasta genbank_seq_acc_meta.fasta --merge 200

#python scripts/genbank_fetcher.py --accession AM180355.1 --locus tcdA tcdB tcdC --bed genbank_coord_locus.bed6 --records genbank_records_locus.txt --fasta genbank_seq_locus.fasta --merge 200
#python scripts/genbank_fetcher.py --meta meta_locus.tsv --bed genbank_coord_locus_meta.bed6 --records genbank_records_locus_meta.txt --fasta genbank_seq_locus_meta.fasta --merge 200
