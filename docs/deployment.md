# Deployment

## Installation Guide

The latest version of **MMASeq** can be installed using [Conda](#conda-installation), [pip](#pip), or directly [from source](#install-from-source).

!!! info
    MMASeq currently supports **Linux** and **macOS** only.

---

### <a name="conda-installation"></a> Conda Installation

Because the workflow is implemented in **Snakemake**, and Snakemake does **not** fully support micromamba, installation via **Conda** is the recommended approach.

To create a clean, isolated environment:

```bash
conda create -n mmaseq -c conda-forge -c bioconda mmaseq
```

Activate the enviroment

```bash
conda activate mmaseq
```

### <a name="pip"></a> pip

You can install the latest stable release of MMASeq via pip:

```bash
pip install mmaseq 
```

### <a name="Install from source"></a> Install from source  

This section is for user who wants to install **MMASeq** directly from source (e.g., for reproducibility or to debug a release).
However, the user must ensure to have the **Snakemake >=9.12.0** dependency is installed in the environment. 

```bash
conda create -n mmaseq -c conda-forge -c bioconda "snakemake>=9.12.0"
pip install pandas
git clone https://github.com/ssi-dk/ssi_analysis_utility
pip install . 
```

If the user wants to work on a development branch:

```bash
conda create -n mmaseq -c conda-forge -c bioconda "snakemake>=9.12.0"
pip install pandas
git clone https://github.com/ssi-dk/ssi_analysis_utility
pip install -e . 
cd ssi_analysis_utility
git switch dev    # work on the develop branch
git checkout -b new-feature-branch  # Create a new branch from dev to work on new changes
``` 



## Running the workflow

The pipeline is designed to be run using the executable `mmaseq.py`, which is the main launcher of the pipeline. The launcher is responsible for parsing the input configurations, setting up the necessary environments and databases, and executing the Snakemake workflow according to the specified parameters. When run for the first time the pipeline  will setup the [conda](#conda) environments and the different [databases](#databases) automatically in the specified deployment directory. The deployment directory is a folder where all the conda environments and databases used during the pipeline execution are stored. By default, it is set to be in the current working directory, but it can be specified by the user when running the pipeline. If the specified deployment directory does not exist, it will be created during execution. If it already exists and correctly specified in the launcher command, the pipeline will check for the presence of the required databases and environments before execution, and it will use them if they are found. This approach allows for efficient management of resources and ensures that the necessary components are readily available for the pipeline to run smoothly.

### Options 

```bash
SYNOPSIS
  Configure and execute MMASeq pipeline
USAGE
  usage: mmaseq [-h] 
                [--samplesheet SAMPLESHEET_FILE] 
                [--input_dir INPUT_DIR] 
                [--deploy_dir DEPLOY_DIR] 
                [--outdir OUTDIR] 
                [--threads THREADS]
                [--config CONFIG] 
                [--test] 
                [--debug]
GENERAL
  --help             This help

  --debug            Add debug messages during execution. Mostly used for development and debugging purposes
INPUT
  --samplesheet XXX  Path to samplesheet TSV used by the pipeline. If the samplesheet doesn't exist, `input_dir` must be specified to create one.
  --input_dir XXX    To be specified if the samplesheet does not yet exist. Input directory will be screened for `.fasta` and `fastq.gz` files, 
                     sample_names will be infered from the detected files, and used to populate a samplesheet. After samplesheet creation, 
                     the pipeline will be executed in dry-run mode (simulated run)
  --deploy_dir XXX   Directory used to deploy databases and conda environments used during pipeline execution. Default is in current working directory.
                     If the directory doesn't exist, it will be created during execution. If the directory already exists, it will be used as is, and 
                     the pipeline will check for the presence of the required databases and environments before execution.
  --config XXX       Configuration file location (Default ./config/config.yaml). If not specified, the config file will be overwritten during 
                     subsequent executions.
  --test             Perform test run of all modules. All conda environments and databases will be created in the deployment directory. 
                     Results will be stored in package folder Test/Results, and config file will be generated in config/Test.yaml
OUTPUT
  --outdir XXX       Directory used for storing analysis results. Default is in current working directory.
RESOURCES
  --threads XXX      Amount of threads (cores) to dedicate for executing the pipeline.
```


### Run the test set

Running the test set is a good way to verify that the pipeline is correctly installed and configured on your system. It will execute a predefined set of analyses using example data, and it will help you familiarize with the expected input and output formats. To run the test set, simply execute the following command:

```bash
mmaseq --test 
```



### Run the pipeline on your data
Run the pipeline with the following command:

```bash
mmaseq --samplesheet path/to/samplesheet.tsv 
       --deploy_dir path/to/deploy_dir 
       --outdir path/to/output_dir
```


## Final report

All the workflow results are stored in the output directory specified by `--outdir`. 
The final output is a file named "all_results.tsv". This file is a concatenated file of all the collected output that can be further processed by the user.
It is in the following format:

| sample_name | tool | filename | organism | row_index | row | value |
| --- | --- | --- | --- | --- | --- | --- |
| Sample1 | resfinder | ResFinder_results_tab.txt | xxx | 1 | Resistance gene | xxx |
| Sample1 | plasmidfinder | results_tab.tsv |  | 1 | Accession number | xxx |
| Sample1 | virulencefinder | results_tab.tsv |  | 1 | Database | xxx |
....

## Configure workflow

This section is dedicated on how to configure and set up the different workflow components such as the input manager, the samplesheet and the species-specific configuration files. 

### Input manager configuration

We refer the main configuration file as "input manager" and contains all the informations regarding input, metadata, output etc...
It is formatted as the following:


```yaml
deploy_dir: folder/where/the/deployed/conda/envs/and/databases/are
outdir: /folder/where/to/store/the/results
samplesheet: path/to/samplesheet.tsv
```

Generally, this file is created automatically during pipeline execution, based on the input parameters provided by the user. However, it can also be created manually before running the pipeline, in which case it will be used as-is. If the file already exists and is correctly specified in the launcher command, the pipeline will not overwrite it.  

This behavior allows the user to maintain full control over configuration settings and tailor them to their specific needs. If the file does not exist, it will be generated during execution using the parameters supplied by the user, ensuring that all required configurations are available for the pipeline to run smoothly.

### Samplesheet configuration

The pipeline employs a samplesheet in the .tsv format to run all the analysis  in batch for the samples.
The format is the following:


| sample_name   |      Illumina_mate1     |    Illumina_mate2       |        Assembly        |          config             |
|---------------|-------------------------|-------------------------|------------------------|-----------------------------|
| Sample1       | path/to/read2.fastq.gz  | path/to/read2.fastq.gz  | path/to/assembly.fasta | path/to/species_config.yaml |
| Sample2       | path/to/read2.fastq.gz  | path/to/read2.fastq.gz  | path/to/assembly.fasta | path/to/species_config.yaml |
---

We included a "samplesheet.tsv" in the "data/" folder. This can be modified by the user to suit its need. 


### Species-specific configuration

Certain species may require different tools than those currently defined or tool-specific parameters that differ from those used in other species. The options which distinguish them are defined in the species-specific configuration files.
The format is typical of a .yaml file, and is the following:

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


The analysis supported and their option are listed in the [supported analysis](supported_analysis.md) section of the documentation. Each species-specific configuration file can be linked to a sample in the samplesheet, allowing for flexible and tailored analysis across different organisms. This modular approach ensures that the pipeline can be easily adapted to a wide range of bacterial species and genomic contexts, making it a versatile tool for microbial genomics research and surveillance:


