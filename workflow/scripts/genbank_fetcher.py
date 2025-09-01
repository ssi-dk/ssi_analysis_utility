import argparse
import sys
from Bio import Entrez, SeqIO
from datetime import datetime
import os

def fetch_and_store_genbank_features(email: str, accession: str, locus_filter=None, output_file=None, bed_file=None, fasta_file=None, merge_distance=None, append=False):
    Entrez.email = email
    print(f"Fetching GenBank record for {accession} from NCBI...")

    try:
        with Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text") as handle:
            record = SeqIO.read(handle, "genbank")
    except Exception as e:
        print(f"Failed to fetch or parse {accession}: {e}")
        return

    output = []
    bed_lines = []
    fasta_regions = []  # Store features as (start, end, strand, name)

    timestamp = datetime.now().strftime("%Y-%m-%d_%H:%M:%S")
    output.append(f"Timestamp: {timestamp}\n")
    output.append(f"Running command: {' '.join(sys.argv)}\n")

    output.append(f"\nRecord ID: {record.id}")
    output.append(f"Description: {record.description}")
    output.append(f"Sequence Length: {len(record.seq)} bp")

    output.append("\nCoding sequence (CDS) coordinates:")
    matched_loci = set()

    for feature in record.features:
        if feature.type != "CDS":
            continue

        gene_name = feature.qualifiers.get("gene", [""])[0].strip()
        product = feature.qualifiers.get("product", [""])[0]
        locus_tag = feature.qualifiers.get("locus_tag", [""])[0]
        old_locus_tag = feature.qualifiers.get("old_locus_tag", [""])[0]
        gene_synonyms = feature.qualifiers.get("gene_synonym", [""])[0]
        note = feature.qualifiers.get("note", [""])[0]

        matched = False
        comments = []
        match_sources = {}

        if locus_filter:
            for target in locus_filter:
                matched_fields = []

                if gene_name == target:
                    matched_fields.append(f"Matched by gene ({target})")
                    match_sources[target] = match_sources.get(target, []) + ["gene"]
                    matched_loci.add(target)

                if target.lower() in product.lower():
                    matched_fields.append(f"Matched by product (partial) ({target})")
                    match_sources[target] = match_sources.get(target, []) + ["product"]
                    matched_loci.add(target)

                if locus_tag == target:
                    matched_fields.append(f"Matched by locus_tag ({target})")
                    match_sources[target] = match_sources.get(target, []) + ["locus_tag"]
                    matched_loci.add(target)

                if old_locus_tag == target:
                    matched_fields.append(f"Matched by old_locus_tag ({target})")
                    match_sources[target] = match_sources.get(target, []) + ["old_locus_tag"]
                    matched_loci.add(target)

                if matched_fields:
                    matched = True
                    comments.extend(matched_fields)

            for target, fields in match_sources.items():
                if len(fields) > 1:
                    output.append(f"Warning: Multiple matches for CDS '{gene_name}' on locus '{target}': {', '.join(fields)}")
        else:
            matched = True

        if not matched:
            continue

        start = int(feature.location.start)
        end = int(feature.location.end)
        strand = "+" if feature.location.strand == 1 else "-"
        gene_label = gene_name or old_locus_tag or "unknown"

        output.append(f" - {gene_label}:")
        output.append(f"     coordinates: {start + 1}..{end}")
        output.append(f"     length: {end - start}")
        output.append(f"     strand: ({strand})")
        output.append(f"     Locus tag: {locus_tag}")
        output.append(f"     Old locus tag: {old_locus_tag}")
        output.append(f"     Gene synonym(s): {gene_synonyms}")
        output.append(f"     Product: {product}")
        output.append(f"     Note: {note}")
        if comments:
            output.append(f"     Comments: {', '.join(comments)}")

        if bed_file:
            bed_lines.append(f"{record.id}\t{start}\t{end}\t{gene_label}\t0\t{strand}")

        if fasta_file:
            fasta_regions.append((start, end, strand, gene_label))

    if locus_filter:
        unmatched = set(locus_filter) - matched_loci
        for missing in sorted(unmatched):
            output.append(f"\n - Gene name '{missing}' not found in any CDS feature.")

    if output_file:
        mode = 'a' if append else 'w'
        with open(output_file, mode) as f:
            f.writelines("\n".join(output) + "\n")

    if bed_file and bed_lines:
        header = "# chrom\tchromStart\tchromEnd\tname\tscore\tstrand\n"
        write_header = not os.path.exists(bed_file) or (not append and os.path.getsize(bed_file) == 0)
        mode = 'a' if append else 'w'
        with open(bed_file, mode) as bedf:
            if write_header:
                bedf.write(header)
            bedf.writelines("\n".join(bed_lines) + "\n")

    if fasta_file and fasta_regions:
        fasta_entries = []

        # Sort by start position
        fasta_regions.sort(key=lambda x: (x[2], x[0]))  # group by strand, then start

        merged = []
        current = [fasta_regions[0]]

        for region in fasta_regions[1:]:
            prev = current[-1]
            if region[2] == prev[2] and (region[0] - prev[1]) <= (merge_distance or 0):
                current.append(region)
            else:
                merged.append(current)
                current = [region]
        merged.append(current)

        for group in merged:
            chrom = record.id
            start = group[0][0]
            end = group[-1][1]
            strand = group[0][2]
            name = "+".join([r[3] for r in group])
            seq = record.seq[start:end]
            fasta_entries.append(f">{chrom}_{name}_{start}_{end}\n{str(seq)}")

        mode = 'a' if append else 'w'
        with open(fasta_file, mode) as fasta_out:
            fasta_out.writelines("\n".join(fasta_entries) + "\n")

def main():
    parser = argparse.ArgumentParser(description="Fetch and print CDS features from a GenBank record.")
    parser.add_argument("-e", "--email", default="Your.Name.Here@example.org",
                        help="Your email address (required by NCBI, default: Your.Name.Here@example.org)")
    parser.add_argument("-a", "--accession", required=True, nargs="+",
                        help="One or more GenBank accession numbers")
    parser.add_argument("-l", "--locus", nargs="*",
                        help="Filter by gene/locus names (e.g., tcdA CD0660 cdtA)")
    parser.add_argument("-o", "--output", required=True,
                        help="Output file to store the results")
    parser.add_argument("--bed", help="Optional BED6 file to write matched gene coordinates")
    parser.add_argument("--fasta", help="Optional FASTA file to write raw genomic sequences of CDS (not strand-corrected)")
    parser.add_argument("--merge", type=int,
                        help="Merge CDS sequences if they are within this number of nucleotides (FASTA only)")
    parser.add_argument("--append", action="store_true",
                    help="Append to output/bed/fasta files instead of overwriting them.")

    args = parser.parse_args()

    for acc in args.accession:
        fetch_and_store_genbank_features(
            email=args.email,
            accession=acc,
            locus_filter=args.locus,
            output_file=args.output,
            bed_file=args.bed,
            fasta_file=args.fasta,
            merge_distance=args.merge,
            append=args.append
        )

if __name__ == "__main__":
    main()

# Example usage:
# python genbank_fetcher.py -a AM180355.1 --locus tcdA tcdC cdtA tcdB CD630_06600 CD0660 cdtB cdtS -o output.txt --bed coord.bed6 --fasta cdifftoxin.fasta --merge 500


