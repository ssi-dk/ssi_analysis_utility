# 🧬 SSI Analysis Utility

**ssi_analysis_utility** is a modular, reproducible **bioinformatics pipeline** built for bacterial genome analysis. 

It integrates various tools for antimicrobial resistance profiling, sequence typing, virulence detection, and more. Designed for flexibility and transparency, it leverages **Snakemake** for workflow management and **nbdev** for literate programming.

---

## 📚 Table of Contents

- [🧬 SSI Analysis Utility](#-ssi-analysis-utility)
- [📚 Table of Contents](#-table-of-contents)
- [📁 Project Structure](#-project-structure)
- [🚀 Getting Started](#-getting-started)

---

## 📁 Project Structure

### 📁 `ssi_analysis_utility/`
The core Python module implementing the logic for sample management, configuration parsing, and pipeline execution. This package is designed to be importable and modular.

### 📁 `workflow/`
Contains all Snakemake-related files: rule definitions, species-specific configs, conda environments, and helper scripts. This is the heart of the pipeline's execution logic.

### 📁 `examples/`
Includes test data, example results, and template samplesheets to demonstrate pipeline functionality and structure.

### 📁 `resources/`
Houses important external datasets and pre-built KMA or tool-specific databases (e.g., for AMRFinder, kmeraligner, etc.).

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
│   ├── 01_sample_manager.ipynb
│   ├── 02_analysis_utility.ipynb
│   ├── 03_pipeline_runner.ipynb
│   ├── _quarto.yml
│   ├── index.ipynb
│   ├── nbdev.yml
│   └── styles.css
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
│   ├── __init__.py
│   ├── _modidx.py
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
    │   ├── E.coli.yaml
    │   └── example_species.yaml
    ├── envs
    │   ├── CHtyper.yaml
    │   ├── DatabaseFetch.yaml
    │   ├── LREfinder.yaml
    │   ├── amrfinder.yaml
    │   ├── assembly_lineage_determination.yaml
    │   ├── cgMLSTfinder.yaml
    │   ├── emm_typing.yaml
    │   ├── htslib.yaml
    │   ├── kleborate.yaml
    │   ├── kmeraligner.yaml
    │   ├── kmerfinder.yaml
    │   ├── mlst.yaml
    │   ├── plasmidfinder.yaml
    │   ├── resfinder.yaml
    │   ├── resistence_gene_detection.yaml
    │   ├── serotypefinder.yaml
    │   ├── skesa.yaml
    │   └── virulencefinder.yaml
    ├── rules
    │   ├── characterizers.smk
    │   ├── db_setups.smk
    │   ├── finders.smk
    │   └── others.smk
    └── scripts
        ├── CHTyper-1.0.py
        ├── Cdiff_KMA.py
        ├── KMA_wrangler.py
        ├── LRE-Finder.py
        ├── __pycache__
        ├── blaster.py
        ├── cgMLST.py
        ├── convert_external_genome.py
        ├── genbank_fetcher.py
        ├── logging_utils.py
        └── thresholds.py
```
</details>

---

## 🚀 Getting Started

### 📦 Requirements
**📝 **:  
- 🧪 [Conda (Miniconda or Anaconda)](https://docs.conda.io/en/latest/miniconda.html) >=24.7.1
- 🧬 [Snakemake](https://snakemake.readthedocs.io/en/stable/) (automatically installed via conda env)

**📝 Prepare a sample sheet**:  
   For inspiration, inspect the example sheet found in [`examples/samplesheet.tsv`](examples/samplesheet.t


---

### 📦 Requirements

 
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
3. **Run the pipeline using Snakemake:**
   
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



# Species specific pipeline

## Overview

This Snakemake pipeline automates genomic analyses for antimcrobial surveillance. It integrates tools used in the Danish antimicrobial surveilance and includes detection of resistance genes, virulence factors, and plasmid replicons, in addition to several characterizations tools such as MLST and SerotypeFinder.



# Setup
## Requirements
* snakemake
* conda>=24.7.1

## Installation

1. (Recommendation) Use Micromamba or the like to set up an environment with the above requirement.

2.  Clone the repository to desired location (e.g. ~/repos/ssi_analysis_utility):

    ``` bash
    git clone https://github.com/ssi-dk/ssi_analysis_utility.git ~/repos/ssi_analysis_utility
    ```

# Execution

## Quick start
For first time use and testing, execute **Step 1** and **Step 2**
1. Navigate into the repository (e.g. `cd ~/repos/ssi_analysis_utility`)

. Run `snakemake` with the `--use-conda` flag and dedicate the desired amount of threads (e.g. `4`) using `--cores 4` 

    ``` bash
    snakemake --use-conda --cores 4
    ```

## Running the Pipeline


1.  Prepare a sample sheet, for inspiration inspect the example sheet found in `examples/samplesheet.tsv` of the repository:

**Required fields:**
* sample_name # Name of the sample, will be used to name result files
* Illumina_read_files # Expects Paired-end read file locations, comma-separated
* assembly_file # Path to sample assembly file, used for blast based analysis
* organism # Used for specifying the sample specific configurations, must follow the same naming conventions as seen for the configuration files located in the repository folder `workflow/configs_species` 

*Future planned fields:*
* Nanopore_read_file

2.  Inspect the species specific configuration files located in the repository folder `workflow/configs_species`.
For each analysis module determine whether to include in the pipe by setting its `status` to `True` for includion or `False` for dismissing

3. Set the overall configuration file located in `config/config.yaml`, for more details inspect the below **Configuration** section

4. Execute Snakemake:
    ```bash
    snakemake --use-conda --cores <N cores>
    ```
5.  Inspect the output results files located at the location specified in the `out_folder` filed of the `config/config.yaml` file. For more details on the output, see **Output Structure** section below 

## Configuration
To configure the pipeline and determine input and output, update the input_manager section of the configuration file found in `config/config.yaml` of the repository.
* Set path to the samplesheet by updating the `path` variable.
* Define the location for where to store and locate databases for the analysis tools using the `database_path` variable.
* Leave the `config_species` variable, only there for debugging and development reasons.
* Determine the location on where to provide output for each analysis by updating the `out_folder` variable.

### Output Structure

Results are organized inside the `out_folder` variable specified within the `config/config.yaml` file, (e.g. Results). Analysis output are stored according to {Sample}/{Tool}:

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

## Running status
Each rule of the species specific pipeline are determined inside the `workflow/rules` folder. Logs are generated according to the rule names, so if any step fails, inspect the log file for further details. Log files are organised in a similar manner to the results and can be located in the `Logs` folder of the reporistory.

# Analysis step
Currently supported tools with automatic setup include

* ResFinder, DisinFinder, PointFinder
* PlasmidFinder
* VirulenceFinder
* AMRFinder
* MLST
* KmerAligner for E. coli OH typing and more

Planned future additions include
* CHTyper
* Kleborate
* EMM-typing & Assembly_lineage_determination   # Currently included but not validated for automatic setup.
