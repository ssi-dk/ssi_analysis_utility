# Tools and Databases


## Supported Tools

This pipeline integrates a broad collection of well-established bioinformatics tools and curated reference databases to enable comprehensive bacterial genome analysis across multiple species and sequencing technologies. 


**AMRFinder**

:   Detects antimicrobial resistance genes, resistance-associated point mutations, 
    and stress response genes using NCBI-curated databases.  

---

**ResFinder**

:   Identifies acquired antimicrobial resistance genes by sequence similarity against the ResFinder database.

---

**Kleborate**

:   Performs species confirmation, sequence typing, virulence profiling, and resistance detection 
    specifically for *Klebsiella pneumoniae* species complexes.

---

**MLST**

:   Assigns multilocus sequence types based on allele profiles from housekeeping genes using PubMLST schemes.

---

**MeningoType**

:   Provides genogrouping and fine typing of *Neisseria meningitidis* based on capsule and antigen genes. 

---

**Shovill**

:   A wrapper for short-read genome assembly that optimizes parameters and automates preprocessing, typically using SPAdes or SKESA.

---

**SerotypeFinder**

:   Predicts bacterial serotypes by identifying serotype-specific genetic markers.

---

**spaTyper**

:   Determines *Staphylococcus aureus* spa types by detecting and classifying repeat patterns in the protein A (spa) gene.

---

**VirulenceFinder**

:   Identifies known virulence-associated genes using curated virulence gene databases.

---

**BLAST**

:   Performs sequence similarity searches for gene detection, confirmation, or custom database screening.

---

**Bowtie2**

:   Aligns short sequencing reads to reference sequences, commonly used for targeted mapping and presence/absence analysis.

---

**Minimap2**

:   Aligns long-read or mixed-read sequencing data to reference genomes or gene databases.

---

**HTSlib**

:   Provides core utilities for handling and processing BAM, SAM, and CRAM alignment files.

---

**SeqSero2**

:   Predicts *Salmonella* serotypes using genomic data, supporting both read- and assembly-based workflows.

---

**SISTR**

:   Performs *Salmonella* in silico serotyping and provides cgMLST-based phylogenetic context.

---

**PlasmidFinder**

:   Detects plasmid replicon sequences to infer plasmid content and incompatibility groups.

---

**kma**

:   Uses k-mer–based matching to rapidly assign genotypes, lineages, or profiles against predefined reference sets.

---



All tools are executed within isolated Conda environments, and their associated databases are versioned and either bundled or automatically downloaded and indexed by the pipeline to ensure reproducibility and consistent results across analyses.

For citations of included tools please see the [references](references.md).