import argparse
import sys
from typing import List, Tuple, Dict, Optional
from Bio import Entrez, SeqIO
from datetime import datetime
import os


# ---------- FINALIZERS ----------

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

    bed_lines = [f"{chrom}\t{s}\t{e}\t{name}\t0\t{strand}"
                 for (s, e, strand, name) in regions]

    header = "# chrom\tchromStart\tchromEnd\tname\tscore\tstrand\n"
    mode = 'a' if append else 'w'
    need_header = True
    if append and os.path.exists(path) and os.path.getsize(path) > 0:
        need_header = False

    with open(path, mode) as bedf:
        if need_header:
            bedf.write(header)
        bedf.write("\n".join(bed_lines) + "\n")


def create_locus_fasta(path: Optional[str], record,  # Biopython SeqRecord
                   regions: List[Tuple[int, int, str, str]],
                   merge_distance: Optional[int],
                   append: bool) -> None:
    """
    Extract the DNA sequence for the desired locus based on the identified regions, and merges adjacent regions belonging to same CDS and strand within a merge_distance.
    The DNA sequence within the .fasta file is within the strand orientation, such that sequences originating from the negative strand is not reverse-complemented.
    """
    if not path or not regions:
        return

    regs = sorted(regions, key=lambda x: (x[2], x[0]))  # (strand, start)

    merged_groups: List[List[Tuple[int, int, str, str]]] = []
    current = [regs[0]]
    for r in regs[1:]:
        prev = current[-1]
        same_strand = (r[2] == prev[2])
        close_enough = ((r[0] - prev[1]) <= (merge_distance or 0))
        if same_strand and close_enough:
            current.append(r)
        else:
            merged_groups.append(current)
            current = [r]
    merged_groups.append(current)

    entries = []
    for group in merged_groups:
        gs = group[0][0]
        ge = group[-1][1]
        gname = "+".join(g[3] for g in group)
        seq = record.seq[gs:ge]
        entries.append(f">{record.id}_{gname}_{gs}_{ge}\n{str(seq)}")

    mode = 'a' if append else 'w'
    with open(path, mode) as fasta_out:
        fasta_out.write("\n".join(entries) + "\n")


# ---------- HELPERS ----------

def parse_range(text: str) -> Tuple[int, int]:
    """Parse 'START-END' (1-based inclusive) into (start, end)."""
    try:
        s, e = text.replace(',', '').split('-')
        s, e = int(s), int(e)
        if s <= 0 or e <= 0 or e < s:
            raise ValueError
        return s, e
    except Exception:
        raise argparse.ArgumentTypeError(f"Invalid range '{text}'. Use START-END with positive integers.")


def parse_name_ranges(items: Optional[List[str]]) -> Dict[str, Tuple[int, int]]:
    """
    Parse ['name:1-200', 'other:300-400'] -> {'name': (1,200), 'other': (300,400)} (lowercased keys).
    """
    out: Dict[str, Tuple[int, int]] = {}
    for item in items or []:
        try:
            name, r = item.split(':', 1)
            out[name.strip().lower()] = parse_range(r)
        except Exception:
            raise argparse.ArgumentTypeError(f"Invalid item '{item}'. Use NAME:START-END")
    return out


def collect_cds_index(record) -> Tuple[List[object], Dict[str, object]]:
    """
    Build an index of CDS features keyed by lowercased gene/locus identifiers.
    Returns: (cds_list, key_map)
      key_map: dict key -> feature (gene, locus_tag, old_locus_tag)
    """
    key_map: Dict[str, object] = {}
    cds_list: List[object] = []
    for f in record.features:
        if f.type != "CDS":
            continue
        cds_list.append(f)
        for qname in ("gene", "locus_tag", "old_locus_tag"):
            for k in f.qualifiers.get(qname, []):
                k = k.strip().lower()
                if k and k not in key_map:
                    key_map[k] = f
    return cds_list, key_map


def match_feature_by_name(name: str, key_map: Dict[str, object], cds_list: List[object]) -> Optional[object]:
    """
    Try exact match on gene/locus/old_locus tag; else fallback to product substring match.
    """
    nm = name.lower()
    if nm in key_map:
        return key_map[nm]
    for f in cds_list:
        product = " ".join(f.qualifiers.get("product", []))
        if product and nm in product.lower():
            return f
    return None


def strand_char(feature) -> str:
    if feature.location.strand == 1:
        return "+"
    if feature.location.strand == -1:
        return "-"
    return "."


def safe_label(feature) -> str:
    gene = feature.qualifiers.get("gene", [""])[0].strip()
    old_tag = feature.qualifiers.get("old_locus_tag", [""])[0]
    locus_tag = feature.qualifiers.get("locus_tag", [""])[0]
    return gene or old_tag or locus_tag or "unknown"


# ---------- CORE ----------

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
    cds_regions: Optional[Dict[str, Tuple[int, int]]] = None,
    acc_regions: Optional[Dict[str, Tuple[int, int]]] = None,
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

    # 2) regions for BED/FASTA (0-based half-open)
    regions: List[Tuple[int, int, str, str]] = []

    # Build CDS index
    cds_list, key_map = collect_cds_index(record)

    if mode == "locus":
        record_info.append("Mode: --locus (full CDS per locus)")
        requested = [x.strip() for x in (locus_list or [])]
        for req in requested:
            f = match_feature_by_name(req, key_map, cds_list)
            if not f:
                record_info.append(f" - '{req}' not matched to any CDS.")
                continue
            start = int(f.location.start)
            end = int(f.location.end)
            label = safe_label(f)
            strand = strand_char(f)
            record_info.append(f" - {label} (matched '{req}'): {start+1}..{end} ({strand}) length={end-start}")
            regions.append((start, end, strand, label))

    elif mode == "cds_region":
        record_info.append("Mode: --cds_region (CDS-relative sub-slices)")
        for req_name, (sub_s, sub_e) in (cds_regions or {}).items():
            f = match_feature_by_name(req_name, key_map, cds_list)
            if not f:
                record_info.append(f" - '{req_name}' not matched to any CDS.")
                continue
            start = int(f.location.start)
            end = int(f.location.end)
            strand = strand_char(f)
            cds_len = end - start

            # Clamp CDS-relative 1-based range to CDS length
            cs = max(1, min(sub_s, cds_len))
            ce = max(1, min(sub_e, cds_len))
            if (cs, ce) != (sub_s, sub_e):
                record_info.append(
                    f"   Warning: {req_name} CDS subrange {sub_s}-{sub_e} clamped to {cs}-{ce} (CDS length {cds_len})."
                )

            # Map CDS 1-based inclusive -> genomic 0-based half-open
            if strand == "+":
                g_start = start + (cs - 1)
                g_end   = start + ce
            elif strand == "-":
                g_start = end - ce
                g_end   = end - (cs - 1)
            else:
                g_start = start + (cs - 1)
                g_end   = start + ce

            base = safe_label(f)
            label = f"{base}_{cs}-{ce}"
            record_info.append(f" - {base} (matched '{req_name}'): CDS {cs}-{ce} -> genomic {g_start+1}..{g_end} ({strand})")
            regions.append((g_start, g_end, strand, label))

    elif mode == "acc_region":
        record_info.append("Mode: --acc_region (genomic/Accession sub-slices within CDS)")
        for req_name, (abs_s, abs_e) in (acc_regions or {}).items():
            f = match_feature_by_name(req_name, key_map, cds_list)
            if not f:
                record_info.append(f" - '{req_name}' not matched to any CDS.")
                continue
            start = int(f.location.start)  # CDS genomic bounds (0-based)
            end = int(f.location.end)

            # Convert provided 1-based inclusive -> 0-based half-open
            g_start_wish = abs_s - 1
            g_end_wish   = abs_e

            # Clamp to CDS bounds
            g_start = max(start, min(g_start_wish, end))
            g_end   = max(start, min(g_end_wish,   end))
            if (g_start, g_end) != (g_start_wish, g_end_wish):
                record_info.append(
                    f"   Warning: {req_name} accession range {abs_s}-{abs_e} clamped to {g_start+1}-{g_end} "
                    f"within CDS {start+1}-{end}."
                )
            if g_end <= g_start:
                record_info.append(f"   Warning: empty slice for {req_name} after clamping; skipping.")
                continue

            base = safe_label(f)
            strand = strand_char(f)
            label = f"{base}_{g_start+1}-{g_end}"
            record_info.append(f" - {base} (matched '{req_name}'): genomic {g_start+1}..{g_end} ({strand})")
            regions.append((g_start, g_end, strand, label))

    else:
        record_info.append("Error: no valid mode selected.")
        finalize_record_info(record_file, record_info, append)
        return

    # 3) finalize outputs
    finalize_record_info(record_file, record_info, append)
    finalize_bed(bed_file, record.id, regions, append)
    finalize_fasta(fasta_file, record, regions, merge_distance, append)


# ---------- CLI ----------

def main():
    parser = argparse.ArgumentParser(description="Fetch GenBank CDS features and export full or partial regions.")

    parser.add_argument("-e", "--email", default="Your.Name.Here@example.org",
                        help="Your email address (required by NCBI).")
    parser.add_argument("-a", "--accession", required=True, nargs="+",
                        help="One or more GenBank accession numbers.")

    # Exactly one of these three must be used
    parser.add_argument("-l", "--locus", nargs="*",
                        help="List of gene/locus names (e.g., tcdA tcdB) → export FULL CDS for each.")
    parser.add_argument("--cds_region", nargs="*", metavar="NAME:START-END",
                        help="CDS-relative slices (1-based, inclusive), e.g., tcdA:1-200 tcdB:300-400.")
    parser.add_argument("--acc_region", nargs="*", metavar="NAME:START-END",
                        help="Accession/genomic slices (1-based, inclusive) within each CDS, "
                             "e.g., tcdA:795843-797843.")

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

    args = parser.parse_args()

    chosen = sum([bool(args.locus), bool(args.cds_region), bool(args.acc_region)])
    if chosen != 1:
        parser.error("You must specify exactly ONE of: --locus, --cds_region, or --acc_region.")

    if args.locus:
        mode = "locus"
        locus_list = args.locus
        cds_regions = None
        acc_regions = None
    elif args.cds_region:
        mode = "cds_region"
        locus_list = None
        cds_regions = parse_name_ranges(args.cds_region)
        acc_regions = None
    else:
        mode = "acc_region"
        locus_list = None
        cds_regions = None
        acc_regions = parse_name_ranges(args.acc_region)

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
            acc_regions=acc_regions
        )


if __name__ == "__main__":
    main()
