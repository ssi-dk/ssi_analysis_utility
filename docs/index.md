# Home

## Overview

**SSI Analysis utility** is a fully modular and  reproducible bioinformatics pipeline for bacterial isolate characterization developed at [Statens Serum Institute](https://www.ssi.dk/) as part of the  [MicrobeSeq Denmark project](https://microbeseq.ssi.dk/) 
The pipeline is designed for flexibility and transparency, and it leverages Snakemake for workflow management. Indeed, it integrates various tools for antimicrobial resistance profiling, sequence typing, virulence detection, and more. 

The main objective is to provide a reproducible, modular framework for bacterial Whole Genome Sequencing (WGS) analysis within national surveillance. Its species-specific design and automated dependency management ensure transparent, scalable, and standardized workflows that can be easily adapted to emerging pathogens and implemented in other laboratory settings.


## Description

Loreum Ipsum


## Repository Structure

### `workflow/`

:   Contains all Snakemake-related files: rule definitions, species-specific configs, conda environments, and helper scripts.

### `examples/`
:   Includes test data and template samplesheets to demonstrate pipeline functionality and structure.

### `config/`
:   Stores the main configuration file, metadata and the results catalogue.

<details>
<summary><strong>Click to expand: Full project tree structure </strong></summary>

```
.
├── LICENSE
├── config/
|   ├── Metadata                      # Metadata Folder
│   │   ├── Cdiff_Toxins_genbank_metafile.tsv
│   │   ├── SNP_metafile.tsv
│   │   ├── deletion_metafiles.tsv
│   │   ├── kma_filter.tsv
│   │   ├── meta_acc.tsv
│   │   ├── meta_cds.tsv
│   │   └── meta_locus.tsv
|   ├── config.yaml                    # Main pipeline config
|   └── results_catalogue.yaml
│   └── config.yaml                    
├── examples/
|   ├── Dataset 
│   |   ├── assemblies                 # Example assemblies
│   │   |   ├── ERR14229029.fasta
│   │   |   └── ERR3528110.fasta
|   │   ├── databases
│   │   |   └── OXAndm.fasta           # OXAndm junction database
|   │   ├── reads    
|   |   |   └── dl_script.sh           # Script to download reads
│   └── samplesheet.tsv                # Example input sheet
│                          
├── workflow/
    ├── Snakefile
    ├── configs_species                # Folder containing all species configuration file 
    │   ├── C_difficile.yaml
    │   ├── E_Faecalis.yaml
    │   ├── E_Faecium.yaml
    │   ├── E_coli.yaml
    │   ├── K_pneumoniae.yaml
    │   ├── N_meningitidis.yaml
    │   ├── S_aureus.yaml
    │   ├── S_enterica.yaml
    │   └── default.yaml
    ├── envs                           # Folder containing all environment configuration file 
    │   ├── amrfinder.yaml
    │   ├── assembly_lineage_determination.yaml
    │   ├── blast.yaml
    │   ├── bowtie2.yaml
    │   ├── emm_typing.yaml
    │   ├── fetch.yaml
    │   ├── htslib.yaml
    │   ├── kleborate.yaml
    │   ├── kmeraligner.yaml
    │   ├── meningotype.yaml
    │   ├── minimap2.yaml
    │   ├── mlst.yaml
    │   ├── plasmidfinder.yaml
    │   ├── python_functions.yaml
    │   ├── resfinder.yaml
    │   ├── seqsero2.yaml
    │   ├── serotypefinder.yaml
    │   ├── shovill.yaml
    │   ├── sistr.yaml
    │   ├── spatyper.yaml
    │   └── virulencefinder.yaml
    ├── rules                          # Folder containing all snakemake rules 
    │   ├── assemblers.smk
    │   ├── characterizers.smk
    │   ├── custom_wranglers.smk
    │   ├── fetch_custom_dbs.smk
    │   ├── finders.smk
    │   ├── htslib_utils.smk
    │   ├── mappers.smk
    │   └── setup_dbs.smk
    └── scripts                        # Folder containing all helper scripts
        ├── KMA_Filter.py
        ├── LRE-Typer.py
        ├── LongTable.py
        ├── Repeat_Identifier.py
        ├── SNP_identifier.py
        ├── SPATyper_V2.py
        ├── blaster.py
        ├── convert_external_genome.py
        ├── deletion_identifier.py
        ├── genbank_fetcher.py
        ├── helper_functions.py
        ├── logging_utils.py
        └── thresholds.py
``` 

</details>


