# Workflow description


## Input data





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
|   ├── Metadata/                      # Metadata Folder
|   ├── config.yaml                    # Main pipeline config
|   └── results_catalogue.yaml
│   └── config.yaml                    
├── examples/
|   ├── Dataset/ 
│   |   ├── assemblies/                # Example assemblies
|   │   ├── databases/                 # OXAndm junction database   
|   │   ├── reads/    
|   |   |   └── dl_script.sh           # Script to download test reads
│   └── samplesheet.tsv                # Example input sheet
│                          
├── workflow/
    ├── Snakefile/
    ├── configs_species/               # Folder containing all species configuration file 
    ├── envs/                          # Folder containing all environments configuration files 
    ├── rules/                         # Folder containing all snakemake rules 
    └── scripts/                       # Folder containing all helper scripts
``` 

</details>
