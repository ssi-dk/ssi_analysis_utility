# 🧬 SSI Analysis Utility

[![Snakemake](https://img.shields.io/badge/snakemake-v7.32-blue?logo=snakemake)](https://snakemake.readthedocs.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](./LICENSE)

**ssi_analysis_utility** is a modular, reproducible **bioinformatics pipeline** built for bacterial genome analysis. 

It integrates various tools for antimicrobial resistance profiling, sequence typing, virulence detection, and more. Designed for flexibility and transparency, it leverages **Snakemake** for workflow management and **nbdev** for literate programming.

---

## 📚 Table of Contents

- [🧬 SSI Analysis Utility](#-ssi-analysis-utility)
- [📚 Table of Contents](#-table-of-contents)
- [📁 Project Structure](#-project-structure)
- [🚀 Getting Started](#-getting-started)
- [🧫 Supported Databases](#-supported-databases)
- [📊 Running Pipeline](#-running-pipeline)
---

## 📁 Project Structure

### 📁 `ssi_analysis_utility/`
The core Python module implements the logic for sample management, configuration parsing, and pipeline execution. This package is designed to be importable and modular.

### 📁 `workflow/`
Contains all Snakemake-related files: rule definitions, species-specific configs, conda environments, and helper scripts.

### 📁 `examples/`
Includes test data, example results, and template samplesheets to demonstrate pipeline functionality and structure.

### 📁 `resources/`
Stores the important external datasets and pre-built KMA or tool-specific databases (e.g., for AMRFinder, kmeraligner, etc.).

### 📁 `nbs/`
Notebook-based documentation and exploratory development environment built using `nbdev`. These include tutorials, demos, and tests.

---

<details>
<summary><strong>📂 CLICK TO EXPAND: FULL PROJECT TREE STRUCTURE WITH EXAMPLES FOR SPECIES SAMPLE</strong></summary>


The structure below shows all information within the directory after running two snakemake pipelines

1. Updating all databases
2. Running species-specific pipeline for one sample [SRR10518319](https://www.ebi.ac.uk/ena/browser/view/SRR10518319) belonging to the species Clostridioides difficile 
```text
.
├── LICENSE
├── Logs
│   ├── Databases
│   │   ├── setup_AMRFinder.log
│   │   ├── setup_CdiffTRST.log
│   │   ├── setup_CdiffToxin.log
│   │   ├── setup_DisinFinder.log
│   │   ├── setup_EcoliKmerAligner.log
│   │   ├── setup_PlasmidFinder.log
│   │   ├── setup_PointFinder.log
│   │   ├── setup_ResFinder.log
│   │   └── setup_VirulenceFinder.log
│   ├── SRR10518319
│   │   ├── Cdiff_KMA_Toxin.log
├── MANIFEST.in
├── README.md
├── conda.dev.env.yaml
├── config
│   └── config.yaml
├── examples
│   ├── Dataset
│   │   ├── assemblies
│   │   │   ├── SRR10518319.fasta
│   │   └── reads
│   │       ├── SRR10518319_1.fastq.gz
│   │       ├── SRR10518319_2.fastq.gz
│   │       └── dl_script.sh
│   ├── Log
│   │   ├── SRR10518319_kma_cdiff.log
│   ├── Results
│   │   ├── SRR10518319
│   │   │   ├── Cdiff_KMA_Toxin
│   │   │   ├── GenotypeCalls
│   │   │   └── skesa
│   ├── samplesheet.tsv
├── input
├── nbs
│   ├── 00_core.ipynb
├── pyproject.toml
├── resources
│   ├── Clostridioides_difficile_db
│   │   ├── TRST
│   │   │   ├── TR10_repeat_sequences.fa
│   │   │   ├── TR10_repeat_types.txt
│   │   │   ├── TR6_repeat_sequences.fa
│   │   │   ├── TR6_repeat_types.txt
│   │   │   └── TRST_repeat_types.txt
│   │   └── Toxin
│   │       ├── Cdiff_Toxin.bed6
│   │       ├── Cdiff_Toxin.comp.b
│   │       ├── Cdiff_Toxin.fasta
│   │       ├── Cdiff_Toxin.fasta.fai
│   │       ├── Cdiff_Toxin.length.b
│   │       ├── Cdiff_Toxin.name
│   │       ├── Cdiff_Toxin.seq.b
│   │       └── Cdiff_Toxin.txt
│   ├── amrfinderplus
│   ├── disinfinder_db
│   ├── kmeraligner
│   ├── plasmidfinder_db
│   ├── pointfinder_db
│   ├── resfinder_db
│   └── virulencefinder_db
├── settings.ini
├── setup.py
├── ssi_analysis_utility
│   ├── analysis_utility.py
│   ├── config
│   │   ├── config.default.env
│   │   └── config.default.yaml
│   ├── core.py
│   ├── pipeline_runner.py
│   └── sample_manager.py
└── workflow
    ├── Snakefile
    ├── configs_species
    │   ├── C.diff.yaml
    ├── envs
    │   ├── DatabaseFetch.yaml
    │   ├── skesa.yaml
    ├── rules
    │   ├── characterizers.smk
    │   ├── db_setups.smk
    │   ├── finders.smk
    │   └── others.smk
    └── scripts
        ├── Cdiff_KMA.py
        ├── genbank_fetcher.py
        ├── logging_utils.py
        └── thresholds.py
```
</details>

---

## 🚀 Getting Started

### 📦 Requirements:
- 🧪 [Conda (Miniconda or Anaconda)](https://docs.conda.io/en/latest/miniconda.html) >=24.7.1
- 🧬 [Snakemake](https://snakemake.readthedocs.io/en/stable/)
  
### 📝 Prepare a sample sheet:  
For inspiration, inspect the example sheet found in [`examples/samplesheet.tsv`](examples/samplesheet.tsv)

| sample_name   | Illumina_read_files                                                                                      | Nanopore_read_file | assembly_file | organism              | variant | notes |
|---------------|-----------------------------------------------------------------------------------------------------------|---------------------|----------------|------------------------|---------|-------|
| SRR10518319   | examples/Dataset/reads/SRR10518319_1.fastq.gz,examples/Dataset/reads/SRR10518319_2.fastq.gz              | Na                  | SRR10518319.fasta| Clostridioides difficile | Na      | ST2   |
---

### 📦 Installation
 
1. **Clone the repository:**
   ```bash
   git clone https://github.com/your-org/ssi_analysis_utility.git
   cd ssi_analysis_utility
   ```

2. **Create and activate the conda environment:**
   ```bash
    conda env create -f conda.dev.env.yaml
    conda activate ssi_analysis_dev
   ```

### 📦 Configuration files

**Input manager configuration:**

The most predominant configuration file ([`ssi_analysis_utility/config/config.yaml`](ssi_analysis_utility/config/config.yaml)) determines the default input files and the desired databases to be updated using a specific environment when running the pipeline initially

```yaml
####################### INPUT MANAGER #######################################
input_manager:
    path: examples/samplesheet.tsv       # Path to the input samplesheet (CSV/TSV file)
    database_path: resources
    config_species: workflow/configs_species/
    out_folder: examples/Results         # Folder where the analysis results will be saved

###################### ANALYSIS & TOOLS SETTINGS ############################
analysis_settings:

    kmeraligner:
        yaml: ../envs/kmeraligner.yaml
        database: kmeraligner

    Clostridioides_difficile_db:
        yaml: ../envs/DatabaseFetch.yaml
        database: Clostridioides_difficile_db

.
.
.
   skesa:
        yaml : ../envs/skesa.yaml
```

**Species-specific configurations:**

Specific species might require unique tools when running the pipeline or different parameters for tools shared across numerous species. The options which distinguish them are defined in the species-specific configuration files, such as *Clostridioides difficile* specific ([`ssi_analysis_utility/workflow/configs_species/C.diff.yaml`](ssi_analysis_utility/workflow/configs_species/C.diff.yaml))


```yaml
analyses_to_run: 

    kmeraligner:
        status : False
        Title : Kmer Aligner on two pair reads
        ID: 90
        additional_option : -matrix 
        database : resources/plasmidfinder_db/enterobacteriales
        wrangler: workflow/scripts/KMA_wrangler.py

    Cdiff_KMA_Toxin:
        status : True
        Title : Kmer Aligner for Clostridium difficile Toxin detection
        ID: 91
        additional_option : -ref_fsa -nf -sam 4 -vcf 2 
        database : resources/Clostridioides_difficile_db/Toxin/
.
.
.
    bcftools_view_filter:
        status: True
        Title: Filter INDELs from mpileup
        ID: 91
        additional_option: -v indels -i 'INFO/DP>10'
        region: AM180355.1_tcdC_804309_805008
```


## 🧫 Supported Databases

When cloning the repository, the databases have not been downloaded. The source for each database is defined in the database setup snakemake file (ssi_analysis_utility/workflow/rules/db_setups.smk). The current database includes:

**AMR Profiling**
- `setup_AMRFinder`: Sets up AMRFinderPlus for AMR gene and mutation detection (multi-species).
- `setup_ResFinder`: Prepares ResFinder for acquired AMR gene identification.
- `setup_PointFinder`: Indexes PointFinder for chromosomal resistance mutations (e.g., *E. coli*, *Salmonella*).
- `setup_DisinFinder`: Installs DisinFinder for disinfectant resistance gene detection.

**Virulence & Plasmids**
- `setup_PlasmidFinder`: Sets up PlasmidFinder for plasmid replicon detection (e.g., *E. coli*).
- `setup_VirulenceFinder`: Sets up VirulenceFinder for virulence gene detection (e.g., *E. coli*).

**Typing & Subtyping**
- `setup_EcoliKmerAligner`: Downloads and indexes k-mer gene set for *E. coli* typing.
- `setup_CdiffToxin`: Builds *Clostridioides difficile* toxin gene database from GenBank.
- `setup_CdiffTRST`: Downloads TRST repeat files for *C. difficile* strain typing.


### 📦 Update databases:

   1. When running the entire pipeline, the species-specific pipeline (based on the organism provided in the sample sheet) is determined, and based on the configuration files, it automatically downloads the required databases for the species-specific analysis.
      
   2. It is also possible for the user to manually download or update all databases
      ```bash
      snakemake --use-conda --cores all setup_all_databases
      ```
    
   3. It is also possible for the user to manually download or update specific databases

      **Clostridioides difficile toxin database**
      ```bash
      snakemake --use-conda --cores 1 setup_CdiffToxin
      ```

## 📊 Running pipeline

### 📦 Run the pipeline using Snakemake:
   
   To execute the full workflow using all available cores:
   ```bash
   snakemake --use-conda --cores all
   ```

   To re-run incomplete or failed steps:
   ```bash
   snakemake --use-conda --cores all --rerun-incomplete
   ```

   To run a specific rule (e.g., cdiff_kma):
   ```bash
   snakemake --use-conda --cores 1 cdiff_kma
   ```

   To preview what will run without executing (dry-run):
   ```bash
   snakemake -n
   ```
### 📂 Output folder structure

Results are organized inside the `out_folder` variable specified within the `config/config.yaml` file, (e.g. Results). Analysis output are stored according to {Sample}/{Tool}:
```text
    Results/
    ├── sample1/
    │   ├── analysis_1/
    │   ├── analysis_2/
    │   └── analysis_3/
    └── sample2/
    |   ├── analysis_1/
    |   ├── analysis_2/
    |   └──...
    └──...
```

### 📝 Species-specific results - extending samplesheet

Depending on the chosen tools used for analysis, individual species require specific information extracted from similar output data files, defined using the species-specific config files. This final species-specific output functions extracts from the Output folder structure the necessary information and extends the original [`examples/samplesheet.tsv`](examples/samplesheet.tsv) file with specific information

   **Clostridioides difficile specific final output**
   ```bash
   python workflow/scripts/Cdiff_KMA.py --samplesheet examples/samplesheet.tsv --outputfile Results/Cdiff_results.tsv
   ```

   **Extended samplesheet output file**

| sample_name   | Illumina_read_files| Nanopore_read_file | assembly_file | organism| variant | notes | tcdA | tcdB | tcdC | cdtAB | verbose | tcdC117 | tcdCdel | deletion_details | TRST | TR6 | TR10 
|---------------|-----------------------------------------------------------------------------------------------------------|---------------------|----------------|------------------------|---------|-------|-------|-------|-------|-------|-------|-------|-------|-------|-------|-------|-------|
| SRR10518319   | examples/Dataset/reads/SRR10518319_1.fastq.gz,examples/Dataset/reads/SRR10518319_2.fastq.gz| Na| SRR10518319.fasta| Clostridioides difficile | Na| ST2 | Positive | Positive | Positive | - | tcdB_283.94_100.00_99.97;tcdA_268.76_100.00_99.85;tcdC_256.69_100.00_99.86 | wt | - | - | tr046 | A001 | B038
---

                                        


