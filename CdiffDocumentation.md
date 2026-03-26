# *Clostridioides difficile*

This species-specific analysis pipeline is used to analyze samples belonging to the *Clostridioides difficile* species and determine the presence of predefined canonical sequence variants, including Tandem Repeat Sequence Typing (TRST) markers, canonical SNPs, and canonical deletions such as 18 bp, 39 bp, and 54 bp deletion events in target loci.

The workflow combines alignment-based read evidence, pileup-based indel evidence, consensus-sequence ambiguity, and optional assembly support to classify variants in a biologically interpretable way.

## Requirements
### Database
#### Download databases
Specific loci are necessary to determine the canonical sequence variants relevant for future serotyping of *Clostridioides difficile*, and can be downloaded using a customized [python script](https://github.com/ssi-dk/MMASeq/blob/dev/workflow/scripts/genbank_fetcher.py).
- The desired regions of one or multiple loci/locus can be downloaded directly from [genbank](https://www.ncbi.nlm.nih.gov/genbank/) based on their accession numbers. The sequences of interest can be downloaded in several ways, depending on the information provided.
  1. The loci/locus can be downloaded in its entirety by simply providing the accession number(s) within the [metafile](https://github.com/ssi-dk/MMASeq/blob/dev/config/Metadata/meta_locus.tsv).
  2. Any coding sequence within a particular locus can be downloaded based on the locus-specific coordinates, within the [metafile](https://github.com/ssi-dk/MMASeq/blob/dev/config/Metadata/meta_cds.tsv)
  3. The loci/locus/genomic region can be downloaded based on the accession number and the chromosomal coordinates within the [metafile](https://github.com/ssi-dk/MMASeq/blob/dev/config/Metadata/meta_acc.tsv).
#### Curated databases
Several distinct databases are required to determine repeat-associated signals, SNPs, and deletions. 
- The relevant loci for the preceding analysis to determine the canonical variants have been curated into a [database](https://github.com/ssi-dk/ssi_analysis_utility_db/blob/main/clostridioides_difficile/Cdiff_Toxins.fasta)
- The TRST markers are derived from combining repeat-sequence patterns from two loci, TR6 and TR10. The sequence of each locus (e.g. [TR10](https://github.com/ssi-dk/ssi_analysis_utility_db/blob/main/clostridioides_difficile/type_repeats/TR10.fasta)) and the [repeat type](https://github.com/ssi-dk/ssi_analysis_utility_db/blob/main/clostridioides_difficile/type_repeats/TR10.txt) based on the locus-specific sequence is defined in the currated database. The databse contain similar information for the TR6 locus, and the [TRST](https://github.com/ssi-dk/ssi_analysis_utility_db/blob/main/clostridioides_difficile/type_repeats/TRST.txt) markers derived from the combined TR6/TR10 pattern.
### Bioinformatical tools
- The component is alignment-based using the tool [KMA](https://bitbucket.org/genomicepidemiology/kma).
- The pileup-based evidence is generated from the aligned files ([SAM/BAM](https://samtools.github.io/hts-specs/SAMv1.pdf)) and processed using the [htslib](https://www.htslib.org/) library involving a suite of programs, i.e [samtools](https://github.com/samtools/samtools) and [bcftools](https://samtools.github.io/bcftools/)
- The alignment and pileup-based evidence is supported by *de novo* generated assemblies obtained either from [SKESA](https://github.com/ncbi/SKESA), [spades](https://github.com/ablab/spades) or [shovill](https://github.com/tseemann/shovill)
### Snakemake environments
- The alignments creating the Sequence Alignment Map format and the consensus sequence utilize the [KMA](https://github.com/ssi-dk/MMASeq/blob/dev/workflow/envs/kmeraligner.yaml)
- The pileup-based functionality for filtering, sorting and indexing of the alignment files, followed by pileup generation, genotype calling and indel filtering are all based on the [htslib](https://github.com/ssi-dk/MMASeq/blob/dev/workflow/envs/htslib.yaml) library environment. 
- The *de novo* assembly tools are all based on an environment for [shovill](https://github.com/ssi-dk/MMASeq/blob/dev/workflow/envs/shovill.yaml), which includes the SKESA and spades internally.
- The following downstream analysis utilizes customized Python scripts wrangling all of these generated results (described later) and is all based on the same [python environment](https://github.com/ssi-dk/MMASeq/blob/dev/workflow/envs/python_functions.yaml)
### Configuration files
To accurately filter the generated alignment files to identify the genes of interest, deleterious regions, and specific SNPs, several configuration files used by the customized Python scripts are required.
- The *KMA* generated alignment results are filtered on the columns *"Template_Coverage", "Template_Identity", "Depth"* using the thresholds defined in the species-specific entries within the [metafile](https://github.com/ssi-dk/MMASeq/blob/dev/config/Metadata/kma_filter.tsv), keeping those genes from the *Template* column for additional downstream analysis. The generated consensus sequence which is used for potential deletion identification, is likewise filtered based on a percentage of undefined nucleotide *N* within the region (column *[consensus_N](https://github.com/ssi-dk/MMASeq/blob/dev/config/Metadata/deletion_metafiles.tsv)*) 
- The SNPs and Deletions are filtered on the genotype information created by the *mpileup, calling* procedure. The [SNPs](https://github.com/ssi-dk/MMASeq/blob/dev/config/Metadata/SNP_metafile.tsv) are solely filtered on the genotype depth, defined in the column *gt_DP*, while the [deletions](https://github.com/ssi-dk/MMASeq/blob/dev/config/Metadata/deletion_metafiles.tsv) are filtered on the bcftools values IDV (*Integer Maximum number of reads supporting an indel*, threshold defined in column *gt_IDV*), IMF (*Float Maximum fraction of reads supporting an indel*, threshold defined in column *gt_IMF*) and the read depth value (column *gt_DP*)

## Analysis
The Clostridioides difficile-specific analysis is performed in separate steps, including detection of locus presence/absence, TRST markers, canonical SNPs, and deletions. In combination, these features can be indicative, but not conclusive, of a genome-based sequence type ([ST](https://github.com/ssi-dk/ssi_analysis_utility_db/blob/main/clostridioides_difficile/Cdiff_ST_table.tsv)).

## Repeat identification
The *cdiff_repeat_identifier* rule in the [pipeline](https://github.com/ssi-dk/MMASeq/blob/dev/workflow/rules/finders.smk) uses a custom [script](https://github.com/ssi-dk/MMASeq/blob/dev/workflow/scripts/Repeat_Identifier.py) to identify repeat types from a genome assembly by scanning contigs for known TR6 and TR10 repeat types, which are defined in metadata tables and the reference repeat fragments are provided as FASTA records, where each header (for example R001) defines a short repeat fragment sequence.

e.g. TR6 repeat [sequence](https://github.com/ssi-dk/ssi_analysis_utility_db/blob/main/clostridioides_difficile/type_repeats/TR6.fasta) sequence
```bash
>R001
CTTGCATACCACTAATAGTGC
>R002
CTTGCATATCACTAATAGTAC
.
.
.
>R100
ATCCATACTACTGATACTGC
>R101
CTTGCATATTGATGATAGTAC
```
The rule reconstructs the full nucleotide sequence for each known TR6 (Axxx) and TR10 (Bxxx) type by concatenating the corresponding fragment sequences from the reference FASTA files (e.g., R044 with R011 ...). The reconstructed repeat-type sequences are then searched against the assembly contigs using exact forward and reverse-complement matching. Any detected TR6 and TR10 types are reported for the sample. 

With the various combinations of the concatenated sequences being defined in the metadata,

e.g. TR6 [types](https://github.com/ssi-dk/ssi_analysis_utility_db/blob/main/clostridioides_difficile/type_repeats/TR6.txt)
```bash
A001,	R044-R011-R051-R007-R023-R024-R052-R053-R004-R019-R018-R019-R018-R022-R035
A002,	R044-R011-R051-R007-R020-R007-R015-R008-R009-R022-R035
.
.
.
A134,	R044-R011-R051-R007-R016-R011-R022-R035
A135,	R024-R015-R019-R084-R095-R007-R008-R028-R082-R024-R008-R007-R008-R007-R008-R007-R008-R028-R026-R028-R026-R028-R038-R071-R022-R035
```

Finally, the identified TR6 and TR10 types are matched against the [TRST](https://github.com/ssi-dk/ssi_analysis_utility_db/blob/main/clostridioides_difficile/type_repeats/TRST.txt) lookup table to assign a combined TRST type (trxxx) where possible.

 e.g.
```bash
ST_COMB	ST_A	ST_B
tr001	A001	B001
tr002	A002	B002
tr003	A003	B003
tr004	A005	B004
tr005	A001	B005
.
.
.
tr287	A134	B123
tr288	A003	B009
tr289	A095	B059
tr290	A133	B032
tr291	A132	B057
```

## SNP identification
The snp_identifier rule in the [pipeline](https://github.com/ssi-dk/MMASeq/blob/dev/workflow/rules/finders.smk) uses a custom [script](https://github.com/ssi-dk/MMASeq/blob/dev/workflow/scripts/SNP_identifier.py) to identify predefined SNPs and deletion events from a sample BCF file using a metadata [table](https://github.com/ssi-dk/MMASeq/blob/dev/config/Metadata/SNP_metafile.tsv) of expected variant positions.

The metadata file specifies the organism, target gene, genomic position, expected reference and alternative allele, and a minimum read-depth threshold (*gt_DP*) for each marker.
```bash
species	gene	position	reference	alternative	gt_DP
Clostridioides difficile	tcdC	117	A	T	10
Clostridioides difficile	tcdC	184	C	T	10
```

For each entry in the metadata table, the script first filters rows to the requested organism, then matches each target gene name to a contig in the BCF header. Then it extracts variant records at the expected position, with a small configurable buffer around the site, and filters these records by minimum depth using the INFO/DP field and *gt_DP* as the threshold.

The script classifies each target position as one of the following:

- the expected SNP, reported in the form *"A117T"*
- a deletion spanning the position, reported as *"Δ117"*
- another SNP at the same position that does not match the expected allele change, reported as *"other"*
- no qualifying variant at the site, reported as *"wt"*

Results are written as a tab-delimited table containing the species, matched contig name, gene, position, and inferred variant state for each configured marker.

## Deletion identification

The deletion_identifier in the [pipeline](https://github.com/ssi-dk/MMASeq/blob/dev/workflow/rules/finders.smk) uses a custom [script](https://github.com/ssi-dk/MMASeq/blob/dev/workflow/scripts/deletion_identifier.py) to identify predefined deletion variants from a sample using multiple analysis results as evidence-based sources: genotype calls in a filtered call BCF, indel calls in a mpileup-derived BCF, ambiguous sequence in a consensus FASTA, and optional CIGAR-based deletion events from an assembly alignment SAM file.

Deletion targets are defined in a metadata [table](https://github.com/ssi-dk/MMASeq/blob/dev/config/Metadata/deletion_metafiles.tsv) containing the organism, gene, expected deletion start and end coordinates, expected deletion length (del_type), minimum genotype-support thresholds (*gt_IMF*, *gt_IDV*, *gt_DP*), and a consensus ambiguity threshold (*consensus_N*). 

For *Clostridioides difficile*, for example:
```bash
species	gene	del_start	del_end	del_type	gt_IMF	gt_IDV	gt_DP	consensus_N
Clostridioides difficile	tcdC	330	347	18	0.10	5	10	70
Clostridioides difficile	tcdC	341	379	39	0.10	5	10	70
Clostridioides difficile	tcdC	313	366	54	0.10	5	10	70
```

For each gene, the script first matches the gene name to a contig in the call BCF header. If the gene is present, the script searches for deletion-like variants (in the bcf file) overlapping the expected region, as defined in the metadata table, and attempts to assign the best-matching canonical deletion type based on overlap with the expected deletion windows in the metadata table.

Canonical deletion assignment requires that the observed deletion overlaps at least a minimum fraction of the expected deletion interval (*--overlap_fraction* parameter with a default 0.6; see rule [deletion_identifier](https://github.com/ssi-dk/MMASeq/blob/dev/workflow/rules/finders.smk)). The coordinates of the defined deletions within the metadata all overlap with one another, so when the deletion variants are identified (**call** or **mpileup** based) and multiple of the investigated deletions could match, the best candidate is selected by prioritizing:
1. The highest fraction of the expected deletion interval covered
2. The highest fraction of the observed deletion is contained within the expected interval
3. Smallest difference between observed and expected deletion length
4. Smaller expected deletion type in the event of a remaining tie

The rule evaluates evidence in a stepwise manner:
1. Call BCF: primary source for deletion assignment, using thresholded genotype evidence (IMF, IDV, DP)
2. Mpileup BCF: fallback source if no call-based deletion is assigned
3. Consensus FASTA: used to measure the proportion of N bases across each expected deletion window as indirect support for a deletion
4. Assembly SAM/CIGAR: optional supporting evidence from deletions observed in read or assembly alignments

If a high-confidence call-based assignment is made, lower-priority evidence is not used to alter it. Otherwise, consensus and assembly evidence may strengthen an existing assignment, provide weaker standalone support, or indicate conflict. Based on the different amounts of evidence, the identified deletions are assigned to a category. 

### Possible categories
The output reports one row for the selected deletion type per gene and assigns one of the following categories describing confidence and support.

### Category 1 — high-confidence canonical deletion
A deletion in the **call BCF** matches a configured canonical deletion and passes the IMF, IDV, and DP thresholds, with an observed deletion length differing by at most 1 bp from the expected length.

<ins>Biological interpretation:<ins>
- strong evidence for the expected deletion
- likely reflects a true biological deletion with good support in the variant calls
- suitable as the most reliable deletion classification in this scheme

### Category 2 — good-confidence canonical deletion
A deletion in the **call BCF** matches a configured canonical deletion and passes thresholds, but the observed deletion length differs by more than 1 bp from the expected canonical length.

<ins>Biological interpretation:<ins>
- Likely a real deletion in the expected region
- May represent a less exact breakpoint call, alignment ambiguity, or a nearby related small deletion event
- Useful evidence, but less precise than Category 1

### Category 3 — threshold-failing call-based deletion
A deletion in the **call BCF** overlaps and matches a configured canonical deletion, but does not pass one or more of the IMF, IDV, or DP thresholds.

<ins>Biological interpretation:<ins>
- A possible deletion signal is present
- May represent a weaker deletion signal, due to low depth or low allele support
- Defined thresholds used for filtering might be too conservative and could be reevaluated
- Reflect a lower-confidence deletion than Category 1 & Category 2

### Category 4 — exact mpileup-supported deletion
No qualifying call-based deletion was found, but the fallback **mpileup BCF** deletion matches a canonical deletion with exact expected length.

<ins>Biological interpretation:<ins>
- There is evidence for the deletion from indel-sensitive pileup calls
- Useful fallback when the main caller did not report a confident deletion
- typically weaker than a good call-based detection, but still consistent with the expected event

### Category 5 — inexact mpileup-supported deletion
No qualifying call-based deletion was found, but the fallback **mpileup BCF** deletion overlaps a canonical deletion with a non-exact length match.

<ins>Biological interpretation:<ins>
- Suggests a deletion in the expected region
- Breakpoint precision is lower, or the observed event may only partially correspond to the canonical deletion
- Should be treated as weaker than Categories 1 - 4

### Category 6 — indirect or assembly-only support
No canonical deletion was assigned from **call BCF** or **mpileup BCF** evidence, but there is supporting evidence from either:

- High N content across an expected deletion interval in the consensus FASTA, generated from a [KMA](https://github.com/genomicepidemiology/kma) alignment. 
- A matching deletion detected from SAM/CIGAR parsing of the assembly alignment, generated from one of multiple [assemblies](https://github.com/ssi-dk/MMASeq/blob/dev/workflow/rules/assemblers.smk) tools.

<ins>Biological interpretation:<ins>
- The region may be disrupted or unresolved in a way consistent with a deletion
- A high N percentage can indicate poor consensus resolution, local assembly uncertainty, low coverage, or ambiguity caused by a true deletion
- Assembly-only support may indicate a structural event that was not confidently represented in the variant calls
- This category is suggestive, but less direct than explicit variant-call evidence

### Category 7 — conflicting evidence
A canonical deletion was assigned from **call BCF** or **mpileup BCF** evidence, but the consensus-based classification supports a different configured deletion type.

<ins>Biological interpretation:<ins>
- Different evidence sources disagree about which deletion is present
- may reflect low-quality data, misalignment, imprecise breakpoint assignment, or ambiguity between overlapping canonical deletion types
- should be interpreted cautiously and not treated as a clean deletion call - and perhaps further investigated

### Upgrade scenarios
Consensus and assembly evidence can strengthen an existing call or mpileup classification when they support the same deletion type, thus increasing the confidence in the assigned deletion type.

### Consensus-based upgrades
If the consensus FASTA supports the same deletion type as the call/mpileup classification:

- The category is improved by 1 level if the N-supported interval is consistent but not exact in length
- The category is improved by 2 levels if the number of N bases exactly matches the expected deletion length

**Examples:**
- Category 3 plus good consensus support may become Category 2 or Category 1
- Category 5 plus exact consensus support may improve to Category 3

<ins>Biological interpretation:<ins>
- A high fraction of N bases in the expected deletion window can indicate that the consensus could not be resolved across the deleted region
- exact-length N tracts provide stronger indirect support than a general ambiguous interval
- However, Ns are still indirect evidence and may also result from poor local sequence quality or assembly uncertainty

### Assembly-based upgrades
If the SAM/CIGAR evidence contains a deletion with the same length as the chosen canonical deletion type, the category is improved by 2 levels, up to a maximum of Category 1.

Examples:
- Category 5 with exact assembly-length support may improve to Category 3
- Category 3 with exact assembly-length support may improve to Category 1

<ins>Biological interpretation:<ins>
- alignment-derived deletions provide an additional signal that the structural event is real
- exact-length support is particularly useful when standard variant calling is weak or ambiguous
- still, the signal depends on alignment quality and may be affected by alignment bias, repetitive or difficult regions
