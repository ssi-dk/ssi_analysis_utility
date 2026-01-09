#!/usr/bin/env python

argparser(
        deployment_dir=None # Will use to store databases AND conda environments. Default: workflowdir/deploy
        database_force_update=y/N, # Removes deployment_dir/Resources
        conda_force_reinstall=y/N, # Removes deployment_dir/conda (If deployment_dir is None, then deployment_dir/.snakemake/conda 
        samplesheet_file=examples/samplesheet.tsv, # Points to file with samplesheet, if not existing then creates it provided dataset_dir is provided [FUTURE ENHANCEMNET]
        outdir=path, # Path to output
        results_everything=y/N, # Rule all --notemp
        results_cleaned=Y/n, # Rule clean -> copy resultcatalougue to new folder, delete old sample folders
        )



def create_config(x,y,z):
    """
    Create configuration file for snakemake.
    """
    pass


def link_assemblies(x,y,z):
    """
    Checks samplesheet for assembly files and artificially creates assembly output to enable use of custom assemblies and skip workflow assembly. Beware to use assembly keyword which corresponds with the one used in species configurations.
    """
    samplesheet = read(samplesheet_file)

    assemblies = {samplesheet.loc("sample_name"): samplesheet.loc("assembly")}

    for sample in assemblies.keys():
        assembly_file = assemblies.get(sample)

        if os.ispath(assembly_file):
            assembly_dir = f"{outdir}/{sample}/Assemblies/"

            try:
                os.createdir(assembly_dir)
            except PermissionErr as err:
                print(f"Unable to create sample folder {assembly_dir} due to permission issues.\n{err}"



            assembly_target = f"{assembly_dir}/{sample}_{type}.fasta"

            try:
                os.makelink(assembly_file, assembly_target)
            except PermissionErr as err:
                print(f"Unable to make symbolic link {assembly_target} due to permission issues.\n{err}"

    pass


def create_samplesheet(x,y,z):
    """
    Checks whether samplesheet file exists, if not it will scan the input dir and create a samplesheet with default settings from all samples
    """
    find_cmd = f"find {sampledir} -iname '.f*a*'"
    files = subprocess.run(find_cmd, check = True)
    print(f"I exepct a lot of files here:\n{files}")
    # TO BE CONTINUED
    # [] Regex patterns for files to ensure fastq and fasta files only
    # [] Distinguish assemblies and reads, use default scheme
    pass


def create_command(commands_list):
    """
    Unstacks the command list into a overall snakemake command.
    """
    # add --config-file flag to smk command to ensure latest config file is used
    # Add --conda-prefix deployment dir to ensure consistent deployment
    pass



