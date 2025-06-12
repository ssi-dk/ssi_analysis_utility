# 🧬 SSI Analysis Utility

[![Snakemake](https://img.shields.io/badge/snakemake-v7.32-blue?logo=snakemake)](https://snakemake.readthedocs.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](./LICENSE)

**ssi_analysis_utility** is a modular, reproducible **bioinformatics pipeline** built for bacterial genome analysis. 

It integrates various tools for antimicrobial resistance profiling, sequence typing, virulence detection, and more. Designed for flexibility and transparency, it leverages **Snakemake** for workflow management and **nbdev** for literate programming.

---

## 📚 Table of Contents

- [🧬 SSI Analysis Utility](#-ssi-analysis-utility)
- [📚 Table of Contents](#-table-of-contents)
- [🚀 Getting Started](#-getting-started)
- [📁 Project Structure](#-project-structure)
- [🧫 Supported Databases](#-supported-databases)
- [📊 Running Pipeline](#-running-pipeline)
- [🧬 Coding: Species Additions](#-coding-species-additions)
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

Results are organized inside the `out_folder` variable specified within the `config/config.yaml` file (e.g., Results). Analysis output are stored according to {Sample}/{Tool}:
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

Depending on the chosen tools used for analysis, individual species require specific information extracted from similar output data files, defined using the species-specific config files. This final species-specific output function extracts from the Output folder structure the necessary information and extends the original [`examples/samplesheet.tsv`](examples/samplesheet.tsv) file with specific information

   **Clostridioides difficile specific final output**
   ```bash
   python workflow/scripts/Cdiff_KMA.py --samplesheet examples/samplesheet.tsv --outputfile Results/Cdiff_results.tsv
   ```

   **Extended samplesheet output file**

| sample_name   | Illumina_read_files| Nanopore_read_file | assembly_file | organism| variant | notes | tcdA | tcdB | tcdC | cdtAB | verbose | tcdC117 | tcdCdel | deletion_details | TRST | TR6 | TR10 
|---------------|-----------------------------------------------------------------------------------------------------------|---------------------|----------------|------------------------|---------|-------|-------|-------|-------|-------|-------|-------|-------|-------|-------|-------|-------|
| SRR10518319   | examples/Dataset/reads/SRR10518319_1.fastq.gz,examples/Dataset/reads/SRR10518319_2.fastq.gz| Na| SRR10518319.fasta| Clostridioides difficile | Na| ST2 | Positive | Positive | Positive | - | tcdB_283.94_100.00_99.97;tcdA_268.76_100.00_99.85;tcdC_256.69_100.00_99.86 | wt | - | - | tr046 | A001 | B038
---


## 🧬 Coding: Species Additions

The snakemake pipeline will always depend on the config files and overall structure, supporting a limited number of species. If a user needs to include a new species, several crucial files must be altered, and specific steps must be taken to update existing files.

This section describes how to extend the pipeline to accommodate new bacterial species through coding with the example of "Clostridioides difficile". 

### 📌 Quick Navigation
- [🧪 Species Databases](#-species-databases)
- [⚙️ Species Configurations](#-species-configurations)
- [🐍 Species Workflow](#-species-workflow)
- [🧾 Species Data Wrangling](#-species-data-wrangling)

---

### 🧪 Species Databases

If a species requires a new database, the following files need to be altered:
1. `config/config.yaml` : add the database 
   * a yaml entry name - new unique name
   * a yaml environment - add dependencies to the `workflow/envs/DatabaseFetch.yaml` or create a separate yaml file 
   * A database name used during the pipeline. The current convention is that the yaml entry name and database name remain identical.

```yaml
Clostridioides_difficile_db:
        yaml: ../envs/DatabaseFetch.yaml
        database: Clostridioides_difficile_db 
```
2. `workflow/rules/db_setups.smk` : Add a new Snakemake rule which downloads the database, which refers to the database specifications defined in step 1
```yaml
rule setup_CdiffToxin:
    conda:
        config["analysis_settings"]["Clostridioides_difficile_db"]["yaml"]
    output:
        database = directory(f'{database_path}/{config["analysis_settings"]["Clostridioides_difficile_db"]["database"]}/Toxin')
    params:
        accession_loci = "AM180355.1:tcdA,tcdB,tcdC AF271719.1:cdtA,cdtB",
        db_toxin = "Cdiff_Toxin"
    log:
        stdout = f'Logs/Databases/setup_CdiffToxin.log'
    message:
        "[setup_CdiffToxin]: Setting up C. difficile toxin database"
    shell:
        r"""
        code to download the database files and create index files if necessary
        """
```
* The "conda:", "output:", "log:", "message:" parts of the Snakemake rule are the conventions for all database download rules. 
3. Add the database setup rule name from step 2 to as an input in the "rule setup_all_databases:" of `workflow/rules/db_setups.smk`.

```yaml
rule setup_all_databases:
    input:
        rules.setup_EcoliKmerAligner.output.database,
        .
        .
        rules.setup_CdiffToxin.output.database,
```

---
### ⚙️ Species Configurations

During development, test data might be required:
1. `examples/Dataset/reads/dl_script.sh` : add the SRA link for simpler download of multiple samples from same or different species
```yaml
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR142/ERR142064/ERR142064_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR142/ERR142064/ERR142064_2.fastq.gz
```
2. `examples/samplesheet.tsv` : add the new sample information as described in [`examples/samplesheet.tsv`](examples/samplesheet.tsv).

**Species-specific information**

3. `workflow/Snakefile` : add the various naming conventions of the species to the "species_name_map" dictionary, with the key representing the "organism" column from the [`examples/samplesheet.tsv`](examples/samplesheet.tsv), and the value representing the prefix of the species-specific configuration file. By defining various naming conventions, the workflow allows for greater flexibility of different spelling and older organism names.
```yaml
species_name_map = {
    "Clostridioides difficile": "C.diff",
    "Clostridium difficile": "C.diff",
    "C. difficile": "C.diff",
    "C difficile": "C.diff",
    "C. diff": "C.diff",
    "C.diff": "C.diff",
    "Escherichia coli": "E.coli",
    "E coli": "E.coli",
    "E.coli": "E.coli",
    # Add more mappings here
}
```
4. `workflow/config_species/{species_name_map.values}` : Create a species-specific configuration file, with the naming convention following the defined values from the "species_name_map" dictionary (described in step 3), e.g. `workflow/configs_species/C.diff.yaml`. Once created, the conventions are as follows:
   * All species specific configuration files must **</u>contain the same</u>** rules
   * Thus, when adding a new unique species-specific snakemake rule to the *C.diff.yaml*, the status is set as *status : True*. The rule should also be copied to the *E. coli.yaml* (and other species) with *status : False*

**</u>C.diff.yaml</u>**
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
```

**</u>E.coli.yaml</u>**
```yaml
analyses_to_run: 

    kmeraligner:
        status : True
        Title : Kmer Aligner on two pair reads
        ID: 90
        additional_option : -matrix 
        database : resources/plasmidfinder_db/enterobacteriales
        wrangler: workflow/scripts/KMA_wrangler.py

    Cdiff_KMA_Toxin:
        status : False
        Title : Kmer Aligner for Clostridium difficile Toxin detection
        ID: 91
        additional_option : -nc -nf -sam 4 -vcf 2
        database : resources/Clostridioides_difficile_db/Toxin/
```

### 🐍 Species Workflow

5. `workflow/rules/others.smk` : Once the status of the species-specific rules have been set, the corresponding rule and expected output is defined.
```yaml
rule Cdiff_KMA_Toxin:
    input:
        R1 = lambda wildcards: sample_to_illumina[wildcards.sample][0],
        R2 = lambda wildcards: sample_to_illumina[wildcards.sample][1],
        database = rules.setup_CdiffToxin.output.database
    params:
        db_prefix = "Cdiff_Toxin",
        add_opt = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["Cdiff_KMA_Toxin"]["additional_option"],
        prefix = "%s/{sample}/Cdiff_KMA_Toxin/{sample}" % OUT_FOLDER
    output:
        aln = temp("%s/{sample}/Cdiff_KMA_Toxin/{sample}.aln" % OUT_FOLDER),
        res = "%s/{sample}/Cdiff_KMA_Toxin/{sample}.res" % OUT_FOLDER,
        fsa = "%s/{sample}/Cdiff_KMA_Toxin/{sample}.fsa" % OUT_FOLDER,
        sam = temp("%s/{sample}/Cdiff_KMA_Toxin/{sample}.sam" % OUT_FOLDER),
        vcf_gz = temp("%s/{sample}/Cdiff_KMA_Toxin/{sample}.vcf.gz" % OUT_FOLDER),
    conda:
        config["analysis_settings"]["Clostridioides_difficile_db"]["yaml"]
    log:
        stdout = 'Logs/{sample}/Cdiff_KMA_Toxin.log'
    message:
        "[Cdiff_KMA_Toxin]: Running Clostridium difficile kmer aligner on {wildcards.sample}"
    shell:
        """
        mkdir -p $(dirname {output.sam})
        cmd="kma -ipe {input.R1} {input.R2} -o {params.prefix} {params.add_opt} -t_db {input.database}/{params.db_prefix} > {output.sam}"
        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """
```
6. `workflow/Snakefile` : Once these species-specific rules and their status have been defined, the expected output from those rules should be added to the **rule all** in the Snakefile, ensuring the correct config rules are initiated. 
```yaml
rule all:
    input:
        expand(
            "%s/{sample}/kmeraligner/" %OUT_FOLDER, 
            sample=[
                s for s in SAMPLES 
                if species_configs[sample_to_organism[s]]["analyses_to_run"]["kmeraligner"]["status"] is True
            ]
        ),
        expand(
            "%s/{sample}/Cdiff_KMA_Toxin/{sample}.sam" % OUT_FOLDER,
            sample=[
                s for s in SAMPLES 
                if species_configs[sample_to_organism[s]]["analyses_to_run"]["Cdiff_KMA_Toxin"]["status"] is True
            ]
        ),
.
.
.
include : "rules/db_setups.smk"
include : "rules/finders.smk"
include : "rules/characterizers.smk"
include : "rules/others.smk"
```
### 🧾 Species Data Wrangling

7. Once the snakemake rules have been completed, the results is present in the defined sample specific folders. Depending on the species, different output information is relevant. As such, each species has its own data wrangling script. e.g `workflow/scripts/Cdiff_wrangler.py` or `workflow/scripts/Ecoli_wrangler.py`.
   * When creating a new data wrangling script for a species, the convention is that it should at least take as input the [`examples/samplesheet.tsv`](examples/samplesheet.tsv)
   * The output should extend the columns of the [`examples/samplesheet.tsv`](examples/samplesheet.tsv) as explained in section [Species-specific results – extending samplesheet](#species-specific-results---extending-samplesheet)

8. `workflow/scripts/thresholds.py` : contains the thresholds used for all species to filter on the input data, such as removing low-quality alignment information etc. This is to ensure that if two species need to apply the same threshold to output files created from the same Snakemake rule, it can easily be imported as opposed to redefined.
```yaml

### KMAfinder thresholds
ecoli_kma_threshold = {
    "stx": [98, 98],
    "wzx": [98, 98],
    "wzy": [98, 98],
    "wzt": [98, 98],
    "wzm": [98, 98],
    "fliC": [90, 90],
    "fli": [90, 90],
    "eae": [95, 95],
    "ehxA": [95, 95],
    "other": [98, 98]
}

cdiff_kma_threshold = {
    "tcdA": [90,90,10],
    "tcdB": [90,90,10],
    "tcdC": [90,90,10],
    "cdtAB": [90,90,10],
    "other": [98, 98,10]
}
.
.
```
