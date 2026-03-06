# Workflow description


## Input data

**MMASeq (Mixed Microbial Analysis on Sequencing data)** accepts as input either raw sequencing reads (FASTQ format) or assembled genomes (FASTA format). At the present state, the pipeline is designed to handle short-read. However, long-read support is planned for future releases. 
The input data is specified through a samplesheet, which allows for batch processing of multiple samples in a single run. Each sample can be associated with its own set of reads and/or assemblies, enabling flexible analysis across different species or strains. The pipeline also supports the use of species-specific configuration files, which can be linked to each sample in the samplesheet to tailor the analysis according to the organism being studied. This modular approach ensures that the pipeline can be easily adapted to a wide range of bacterial species and genomic contexts, making it a versatile tool for microbial genomics research and surveillance. 


## Repository Structure

Within the `mmaseq` folder repository, the project is organized into several key directories and files that facilitate the workflow and its management:

### `workflow/`

:   Contains all Snakemake-related files: rule definitions, species-specific configs, conda environments, and helper scripts.

### `data/`
:   Includes test data and template samplesheets to demonstrate pipeline functionality and structure.

### `config/`
:   Stores the main configuration file, metadata and the results catalogue.

<details>
<summary><strong>Click to expand: Full project tree structure </strong></summary>

```
.
├── LICENSE
├── README.md
├── mkdocs.yml
├── pyproject.toml
├── mmaseq/
    ├── config/
        |   ├── target_screening/              # Metadata Folder
        |   ├── species_configs/               # Folder containing all species configuration file
        |   ├── results_catalogue.yaml         # File containing all the information regarding the                  
        |   ├── config.yaml                    # Main pipeline default configuration file                 
        ├── data/
        │   ├── assemblies/                    # Example assemblies  
        |   ├── reads/                         # Reads folder
        |   |    └── dl_script.sh              # Script to download test reads
        │   └── samplesheet.tsv                # Example input sheet
        │                          
        ├── workflow/
        │   ├── Snakefile                      # Main Snakemake workflow definition
        │   ├── envs/                          # Folder containing all environments configuration files 
        │   ├── rules/                         # Folder containing all snakemake rules 
        │   └── scripts/                       # Folder containing all helper scripts
        ├── helper_functions.py                # Helper functions used in the launcher
        └── mmaseq.py                          # Main launcher of the pipeline
        
``` 

</details>
