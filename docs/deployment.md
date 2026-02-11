
## Installation Guide

The latest version of **SSI Analysis Utility** can be installed using [Conda](#conda-installation), [pip](#pip), or directly [from source](#install-from-source).

!!! info
    SSI Analysis Utility currently supports **Linux** and **macOS** only.

---

### <a name="conda-installation"></a> Conda Installation

Because the workflow is implemented in **Snakemake**, and Snakemake does **not** fully support micromamba, installation via **Conda** is the recommended approach.

To create a clean, isolated environment:

```bash
conda create -n ssi_analysis_utility -c conda-forge -c bioconda ssi_analysis_utility
```

Activate the enviroment

```bash
conda activate ssi_analysis_utility
```

### <a name="pip"></a> pip

You can install the latest stable release of SSI Analysis Utility via pip:

```bash
pip install ssi_analysis_utility 
```

### <a name="Install from source"></a> Install from source  

This section is for user who wants to install **SSI Analysis utility** directly from source (e.g., for reproducibility or to debug a release).
However, the user must ensure to have the **Snakemake >=9.12.0** dependency is installed in the environment. 

```bash
conda create -n ssi_analysis_utility -c conda-forge -c bioconda "snakemake>=9.12.0"
git clone https://github.com/your-org/ssi_analysis_utility.git
```

If the user wants to work on a development branch:

```bash
conda create -n ssi_analysis_utility -c conda-forge -c bioconda "snakemake>=9.12.0"
git clone https://github.com/your-org/ssi_analysis_utility.git
cd ssi_analysis_utility
git switch dev    # work on the develop branch
git checkout -b new-feature-branch  # Create a new branch from dev to work on new changes
``` 


## Configure workflow

This section is dedicated on how to configure and set up the different workflow components. 
The pipeline is designed to work using configurations and samplesheet files which are called during its execution.  

### Input manager configuration

We refer the main configuration file as "input manager" and contains all the informations regarding input, metadata, output etc...
It is formatted as the following:

TO BE CHANGED

```yaml
samplesheet: path/to/samplesheet.tsv
database_path: path/to/databases 
config_species: path/to/configs_species
metadata_path: path/to/Metadata
illumina_read_path: path/to/read
nanopore_read_path: path/to/read
assembly_path: path/to/assembly
output_folder: path/where/to/store/results  
temp_storage: examples/Dataset/databases
result_catalogue: config/results_catalogue.yaml
```

### Samplesheet configuration

The pipeline employs a samplesheet in the .tsv format to run all the analysis  in batch for the samples.
The format is the following:


| sample_name   |      Illumina_mate1     |    Illumina_mate2       |        Assembly        |          config             |
|---------------|-------------------------|-------------------------|------------------------|-----------------------------|
| Sample1       | path/to/read2.fastq.gz  | path/to/read2.fastq.gz  | path/to/assembly.fasta | path/to/species_config.yaml |
| Sample2       | path/to/read2.fastq.gz  | path/to/read2.fastq.gz  | path/to/assembly.fasta | path/to/species_config.yaml |
---

We included a "samplesheet.tsv" in the "examples" folder. This can be modified by the user to suit its need.  


### Species-specific configuration

Certain species may require different tools than those currently defined or tool-specific parameters that differ from those used in other species. The options which distinguish them are defined in the species-specific configuration files.
The format is the following:

```yaml
mlst:
    assemblers: shovill  
    
resfinder:
    options: --species 'Other'

plasmidfinder:
    options: ""

virulencefinder:
    options: ""

serotypefinder:
    options: ""
          
amrfinder:
    assemblers: shovill
    options: ""
    
meningotype:
    assemblers: shovill

kmeraligner:
    database: elmDB
    options : -ID 80 -1t1 -cge
```

The analysis supported and their option are the following:

<details>
<summary><strong>Click to expand: Analysis supported </strong></summary>

```
mlst:
    assemblers: shovill  
    
resfinder:
    options: --species 'Other'

plasmidfinder:
    options: ""

virulencefinder:
    options: ""

serotypefinder:
    options: ""
          
amrfinder:
    assemblers: shovill
    options: ""
    
meningotype:
    assemblers: shovill

kmeraligner:
    database: elmDB
    options : -ID 80 -1t1 -cge

### Typical of C.Difficile
kma_filter:
    options: --organism 'Clostridioides difficile'
    database: Cdiff_Toxins

snp_identifier:
    options: --organism 'Clostridioides difficile' --snp_region_buffer 5
    database: Cdiff_Toxins

deletion_identifier:
    options: --organism 'Clostridioides difficile' --deletion_region_buffer 5 --overlap_fraction 0.5
    database: Cdiff_Toxins
    assemblers: [spades]        # <--- NEW: matches the assembler wildcard

assembly_minimap2:
    options: -ax asm5 --cs
    database: Cdiff_Toxins
    assemblers: [spades]        # used to configure which assembler(s) to align

cdiff_repeat_identifier:
    repeats: TR6 TR10
    combos: TRST
    assemblers: [spades, skesa]

samtools:
    view_options: -q 20
    sort_options: -n

bcftools:
    view_options: -v indels -i 'INFO/DP>10'

# Typical of E.Faecalis and E.Faecium
LREfinder:
    database : [elmDB]

# Typical of E.Coli
chtyper:
    database: fumCH_db 

custom_blaster:
    options: "-perc_identity 90.0"
    assemblers: shovill
    database : OXAndm

# Typical of K.Pneumoniae
kleborate:
    options: --preset kpsc
    assemblers: shovill  

#  Typical of N.Meningitis
meningotype:
    assemblers: shovill

# Typical of S.Aureus
spatyper:
  assemblers: [spades]    

# S.Enterica
seqsero2: true

sistr:
  assemblers: shovill  

``` 

</details>



## Setup 

Generally, the pipeline will setup the [conda](#conda) environments and the different [databases](#databases) automatically in the ".snakemake" folder located in the "ssi_utility/" , however they may be situation in which the user/s would like to set up those in different folders.


### <a name="conda"></a> Setup Conda environments 

To setup the conda environment in a folder of choice run:

```bash
PLACEHOLDER
```


### <a name="databases"></a> Setup Databases

To setup the databases in a folder of choice run:

```bash
PLACEHOLDER
```

## Running the workflow

Once input, conda environments and databases are ready just run the pipeline with the following command:

```bash
PLACEHOLDER
```


## Final report

All the workflow results are stored in "xxx". 
The final output is a file named "xxx". This file is a concatenated file of all the collected output that can be further processed by the user.
It is in the following format:

```bash
PLACEHOLDER
```