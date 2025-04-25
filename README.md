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