# 🧬 genbank_fetcher.py

**Fetch and extract genomic or CDS regions from GenBank accessions**

This script is part of the [SSI Analysis Utility](https://github.com/ssi-dk/ssi_analysis_utility) workflow.  
It retrieves GenBank records via NCBI Entrez, identifies coding sequences (CDS), and outputs specified loci or genomic slices as BED and FASTA files.  
It supports multiple operating modes and batch processing through a metadata file.

---

## 🧩 Overview

`genbank_fetcher.py` enables dynamic and generalized fetching of gene loci from GenBank, removing the need for hardcoded gene lists in Snakemake workflows.  
The script can extract:

- **Full loci** (complete CDS regions)  
- **CDS-relative slices** (subregions within a coding sequence)  
- **Genomic slices** (absolute coordinates)  
- **Metadata-driven batch extraction** (via a TSV metafile)

Developed under PR [#73](https://github.com/ssi-dk/ssi_analysis_utility/pull/73).

---

## ⚙️ Dependencies

- Python ≥ 3.8  
- [Biopython](https://biopython.org/) (Entrez & SeqIO)  
- [pandas](https://pandas.pydata.org/)  
- Standard library: `argparse`, `csv`, `logging`, `pathlib`, etc.

---

## 🚀 Usage

### Basic CLI syntax

```bash
python scripts/genbank_fetcher.py     -e <email>     --accession <ACC>     (--locus | --cds_region | --acc_region) <arguments>     --bed out.bed --fasta out.fasta --records out.txt
```

### Metafile mode (batch processing)

Instead of passing arguments manually, provide a **TSV metafile**:

```bash
python scripts/genbank_fetcher.py     -e <email>     --metafile data/Cdiff_Toxins_genbank_metafile.tsv     --bed out.bed --fasta out.fasta --records out.txt
```

Each row in the metafile describes one target locus or region.

---

## 📄 Metafile Format

| Column      | Description |
|-------------|-------------|
| `accession` | GenBank accession (required) |
| `gene`      | Target gene or locus name |
| `start`     | Start coordinate (CDS-relative or genomic) |
| `end`       | End coordinate (CDS-relative or genomic) |
| `organism`  | (Optional) Organism name for validation |

Depending on which columns are filled, the script infers whether the row is **locus**, **CDS-relative**, or **accession-region** mode.

---

## 🧠 Modes of Operation

| Mode | Flag | Behavior |
|------|------|----------|
| **Locus mode** | `--locus gene1 gene2` | Extract full CDS regions of specified genes |
| **CDS-relative** | `--cds_region gene:START-END` | Extract a subsequence within a CDS (1-based inclusive) |
| **Accession-region** | `--acc_region gene:START-END` | Extract absolute genomic coordinates (1-based inclusive) |
| **Metafile mode** | `--metafile <file.tsv>` | Batch process mixed operations from a TSV |

> Exactly one of `--locus`, `--cds_region`, or `--acc_region` must be provided unless using `--metafile`.

---

## 🧾 Key Arguments

| Flag | Description |
|------|-------------|
| `-e, --email` | Required by NCBI Entrez |
| `-a, --accession` | One or more accession IDs |
| `-l, --locus` | Gene/locus names for extraction |
| `--cds_region` | CDS-relative region, e.g. `toxA:1-200` |
| `--acc_region` | Absolute region, e.g. `toxA:795942-803175` |
| `--metafile` | TSV specifying multiple accessions and genes |
| `--bed` | Output BED6 file path |
| `--fasta` | Output FASTA file path |
| `--records` | Text summary report path |
| `--merge` | Merge nearby slices (gap ≤ N bp) |
| `--append` | Append to existing output files |
| `--no_strand_correct` | Disable strand-based reverse-complementing |

Important constraints:

- `--metafile` **cannot** be combined with `--accession`, `--locus`, `--cds_region`, or `--acc_region`.  
- If multiple accessions are processed and `--append` is not used, outputs may be overwritten.

---

## 🧬 Workflow Summary

1. **Parse inputs** — CLI arguments or metafile rows.  
2. **Fetch GenBank record** — via `Entrez.efetch` and `SeqIO.read`.  
3. **Build CDS index** — map gene names, `locus_tag`, and `old_locus_tag` to CDS features.  
4. **Match features** — exact match on gene/locus tags, fallback to substring match in `product` qualifier.  
5. **Derive genomic coordinates** — based on mode; map CDS-relative coordinates to genomic positions considering strand.  
6. **Generate outputs** — BED6, FASTA, and a textual records summary.  
7. **Log results & warnings** — missing features, clamping, organism mismatches, Entrez errors.

---

## 📦 Outputs

### BED file (`--bed`)
- Standard 6-column BED:  
  `chrom  start  end  name  .  strand`  
- Adjacent slices on the same strand may be merged if gap ≤ `--merge`.

### FASTA file (`--fasta`)
- Extracted sequence slices per locus.  
- Header format example: `>accession_gene_start_end_strand_(rc)`  
- Negative-strand sequences are reverse-complemented by default (unless `--no_strand_correct`).

### Records file (`--records`)
- Human-readable summary per accession:
  - Timestamp & command used  
  - Accession, description, sequence length  
  - Per-gene slice details (requested vs. actual coordinates, clamping notes)  
  - Warnings or errors (e.g., no matching CDS, organism mismatch)

---

## ⚠️ Edge Cases & Notes

- Rows missing `accession` or `gene` are skipped (with a warning).  
- If `organism` is provided in the metafile but does not match the GenBank record (case-insensitive substring check), that row is skipped.  
- Requested coordinates outside CDS/sequence boundaries are clamped and logged.  
- Regions on unknown strand `.` are not reverse-complemented and are merged as-is.  
- Ambiguous gene names use fallback matching on `product` text—this may produce false positives for short or generic names.

---

## 💡 Example Commands

Extract full loci:

```bash
python genbank_fetcher.py   -e user@domain.com   --accession AM180355.1   --locus tcdA tcdB tcdC   --bed toxins.bed --fasta toxins.fasta --records toxins.txt
```

Extract CDS subregions:

```bash
python genbank_fetcher.py   -e user@domain.com   --accession AM180355.1   --cds_region tcdA:100-200 tcdB:1-150   --bed regions.bed --fasta regions.fasta --records regions.txt
```

Batch via metafile:

```bash
python genbank_fetcher.py   -e user@domain.com   --metafile Cdiff_Toxins_genbank_metafile.tsv   --bed meta.bed --fasta meta.fasta --records meta.txt
```

---

## 🔗 Integration with Snakemake

This script is intended to be used by a generalized `fetch_genbank` rule (see PR [#73](https://github.com/ssi-dk/ssi_analysis_utility/pull/73)).  
Example Snakemake snippet:

```python
rule fetch_genbank:
    input:  "data/{species}_genbank_metafile.tsv"
    output: "results/{species}/genbank/{species}.bed",
            "results/{species}/genbank/{species}.fasta",
            "results/{species}/genbank/{species}_records.txt"
    params: email=config["entrez_email"]
    script: "workflow/scripts/genbank_fetcher.py"
```

Using a metafile lets you avoid hardcoding loci in the workflow and supports species-specific configurations.

---

## 🧭 Future Improvements

- Smarter gene matching (synonyms, fuzzy/regex matching).  
- Retry/backoff and caching for Entrez fetches.  
- Configurable metafile path via species config to avoid rule edits.  
- Improved reporting (JSON output for downstream parsing).

---

**Author(s):** SSI Analysis Utility contributors  
**Source:** [`workflow/scripts/genbank_fetcher.py`](https://github.com/ssi-dk/ssi_analysis_utility/blob/dev/workflow/scripts/genbank_fetcher.py)  
**Pull Request:** [#73](https://github.com/ssi-dk/ssi_analysis_utility/pull/73)

---

*This documentation was generated with **ChatGPT (GPT-5)**.*
