#!/usr/bin/env python

argparser(
        deployment_dir=None # Will use to store databases AND conda environments. Default: workflowdir/deploy
        database_force_update=y/N, # Removes deployment_dir/Resources
        conda_force_reinstall=y/N, # Removes deployment_dir/conda (If deployment_dir is None, then deployment_dir/.snakemake/conda 
        samplesheet_file=examples/samplesheet.tsv, # Points to file with samplesheet, if not existing then creates it provided dataset_dir is provided [FUTURE ENHANCEMNET]
        input_dir=None, # Used to create samplesheet, in case samplesheet file doesn't exists. Will be used to scan the input dir recursively, and register paired end reads and assembly files automatically
        outdir=path, # Path to output
        results_everything=y/N, # Rule all --notemp
        results_cleaned=Y/n, # Rule clean -> copy resultcatalougue to new folder, delete old sample folders
        )



def create_config(x,y,z):
    """
    Creates configuration file for snakemake. (config/config.yaml)
    
    # Values and defaults
    samplesheet: data/samplesheet.tsv
    deployment_dir: Deploy
    output_folder: Test/Results
    """
    pass


def link_assemblies(x,y,z):
    """
    Checks config[samplesheet] assembly variable and config variable to determine preexisting assemblies.
    If the files exists, assemblies will be artificially created using symbolic links, to circumvent the pipeline to perform assembly itself.
    The samplesheet config variable will be used to extract the assembly keywords of the config/analysis/species_configs/*yaml files,
    and a symbolic link will be created to the rules assembly: output location for the given samples.
    """
    pass


def create_samplesheet(samplesheet_file, input_dir):
    """
    Will scan the input dir recursively, and register paired end reads and assembly files automatically.
    A samplesheet.tsv will be created at samplesheet_file location, it contains sample_name, read1, read2, assembly, and config variables
    sample_name is devised from either {sample_name}_1.fastq.gz or {sample_name}_R1_*.fastq.gz of read1 OR alternately, from {sample_name}.fasta
    """
    pass


def create_command(commands_list):
    """
    Unstacks the command list into a overall snakemake command.
    """
    # add --config-file flag to smk command to ensure latest config file is used
    # Add --conda-prefix deployment dir to ensure consistent deployment
    pass
