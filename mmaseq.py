#!/usr/bin/env python

import argparse
import logging
from pathlib import Path
import pandas
import re
import os
from datetime import datetime
import yaml
from workflow.scripts.helper_functions import determine_sample_configs 
import subprocess
import shutil
import sys


def parse_arguments():
    parser = argparse.ArgumentParser(
        description = "Configure and execute Seq And Destroy pipeline"
    )

    parser.add_argument(
        "--samplesheet",
        dest = "samplesheet_file",
        default = None,
        help = """
        Path to samplesheet TSV used by the pipeline. (Mandatory)
        If the samplesheet doesn't exist, `input_dir` must be specified to create one.
        """
    )

    parser.add_argument(
        "--input_dir",
        dest = "input_dir",
        default = None,
        help = """
        Input directory MUST be specified if the samplesheet does not yet exist.
        Input directory will be screened for `.fasta` and `fastq.gz` files,
        sample_names will be infered from the detected files, and used to populate a samplesheet.
        After samplesheet creation, the pipeline will be executed in dry-run mode (simulated run)
        """
    )

    parser.add_argument(
        "--deploy_dir",
        dest = "deploy_dir",
        default = None,
        help = """
            Directory used to deploy databases and conda environments used during peipeline execution. (Default: ./Deploy)
            To reinstall conda environments remove the folder: {deployment_dir}/conda
            To fetch the latest databases remove the folder: {deployment_dir}/Databases
        """
    )

    parser.add_argument(
        "--outdir",
        dest = "outdir",
        default = None,
        help = """
            Directory used for storing analysis results. (Mandatory)
        """
    )

    parser.add_argument(
        "--threads",
        dest = "threads",
        default = 4,
        help = """
            Amount of threads (cores) to dedicate for executing the pipeline. (Default 4)
        """
    )

    parser.add_argument(
        "--config",
        dest = "config",
        default = None,
        help = """
            Configuration file location. (Default ./config/config.yaml)
            If not specified, the config file will be overwritten during subsequent executions.
        """
    )

    parser.add_argument(
        "--test",
        dest = "test",
        action = "store_true",
        help = "Perform test run of all modules.\n" +
        "As a side effect, all conda environments and databases will be created in the deployment directory.\n" +
        "Results will be stored in Test/Results, and config file will be generated in config/Test.yaml "
    )

    parser.add_argument(
        "--debug",
        dest = "debug",
        action = "store_true",
        help = """
            Add debug messages during execution. Mostly used for development and debugging purposes
        """
    )

    return parser.parse_args()


def create_samplesheet(samplesheet_file, input_dir):
    """
    Will scan the input dir recursively, and register paired end reads and assembly files automatically.
    A samplesheet.tsv will be created at samplesheet_file location, it contains sample_name, read1, read2, assembly, and config variables
    sample_name is devised from either {sample_name}_1.fastq.gz or {sample_name}_R1_*.fastq.gz of read1 OR alternately, from {sample_name}.fasta
    """
    if not os.path.isdir(input_dir):
        logger.error(f"Input directory doesn't exist. Aborting!\n - {input_dir}")
        sys.exit(1)

    root = Path(input_dir).expanduser().resolve()

    records = {}

    for path in root.rglob("*"):
        if not path.is_file():
            continue

        fname = path.name

        if fname.endswith(".fasta"):
            logger.debug(f"Assembly file detected: {fname}")

            # Extract sample name
            sample = fname.replace(".fasta", "")
            records.setdefault(sample, {"read1": "NA", "read2": "NA", "assembly": "NA", "config": "default.yaml"})
            records[sample]["assembly"] = str(path)

        elif fname.endswith(".fastq.gz"):
            logger.debug(f"Read file detected: {fname}")

            # remove read indicators to get sample name
            sample = re.sub(r'(_R?[12].*)\.fastq\.gz$', '', fname)
            records.setdefault(sample, {"read1": "NA", "read2": "NA", "assembly": "NA", "config": "default.yaml"})

            if re.search(r'(_R?1|_1)', fname):
                records[sample]["read1"] = str(path)
            elif re.search(r'(_R?2|_2)', fname):
                records[sample]["read2"] = str(path)

    samplesheet = (
        pandas.DataFrame.from_dict(records, orient="index")
        .reset_index()
        .rename(columns={"index": "sample_name"})
    )

    samplesheet.to_csv(samplesheet_file, sep = "\t", index = False)

    return None


def create_config(samplesheet_file, outdir, deploy_dir, config_file, force):
    """
    Creates configuration file for snakemake.
    """    

    # Check for existing config file
    if not os.path.isfile(config_file):
        config_dir = os.path.dirname(config_file)

        # Check for existing config dir
        if not os.path.isdir(config_dir):
            logging.info("Config directory does not exist, creating directory")
            os.makedirs(config_dir)

    elif not force:
        timestamp = datetime.now().strftime("%y_%m_%d-%H_%M")
        config_file = f"{os.path.dirname(config_file)}/{timestamp}_{os.path.basename(config_file)}"
        logging.warning(f"Configuration file already exists, will change current config_file to:\n - {config_file}")

    # Cleanup outdir
    outdir = os.path.abspath(outdir).rstrip("/")
    
    config = {"samplesheet": samplesheet_file, "deploy_dir": deploy_dir, "outdir": outdir}

    with open(config_file, 'w') as config_yaml:
        yaml.dump(config, config_yaml)

    return config_file


def link_assemblies(samplesheet_file, config_dir, outdir):
    logger.debug("Reading samplesheet")
    samplesheet = pandas.read_csv(samplesheet_file, sep='\t').set_index("sample_name")

    logger.debug("Importing sample configs")
    sample_configs = determine_sample_configs(samplesheet, config_dir)

    logger.debug("Initiating symlink generation")
    # Iterate through sample specific configs
    for sample, configs in sample_configs.items():
        
        # Extract assembly column from samplesheet
        assembly_sheet = samplesheet.at[sample, "assembly"]
        assembly_source = os.path.abspath(assembly_sheet)

        # Ensure that assembly file exists
        if os.path.isfile(assembly_source):
            
            # Iterate through each module
            for module, options in configs.items():

                # Ignore modules without extra settings
                if type(options) == dict:

                    # Ensure that assembly type is specified in settings
                    if "assemblers" in options.keys():
                        assemblers = options.get("assemblers")

                        # Handle multiple assembly types
                        if isinstance(assemblers, str):
                            assemblers = [assemblers]

                        # Iterate through assembly types
                        for assembly_type in assemblers:

                            # Handle symlink directory
                            assembly_dir = f"{outdir}/{sample}/{assembly_type}"

                            # Ensure target directory exists
                            if not os.path.isdir(assembly_dir):
                                os.makedirs(assembly_dir)

                            # Handle assembly symlink
                            assembly_destination = f"{assembly_dir}/{sample}.fasta"    

                            # Remove preexisting links
                            if os.path.islink(assembly_destination):
                                logger.debug("Symlink already exists. Removing old link!")
                                os.unlink(assembly_destination)
                            
                            # Generate links if nothing exists at destination
                            if not os.path.isfile(assembly_destination):
                                logger.debug(f"Creating symlink:\n - {assembly_source} -> {assembly_destination}")
                                os.symlink(assembly_source, assembly_destination)
                            else:
                                logger.warning(f"File exists and is not a symbolic link:\n - {assembly_destination}")

        elif assembly_sheet.lower() not in ["na", ""]:
            logger.warning(f"Assembly file specified in samplesheet could not be found:\n - {assembly_source}")
    return None


def create_command(threads, config, conda_dir, arguments = None, rules = None):
    """
    Unstacks the command list into a overall snakemake command.
    """
    additionals = ""
    if arguments is not None:
        additionals = " ".join(arguments)

    target_rules = ""
    if rules is not None:
        target_rules = " ".join(rules)

    logger.debug(
        f"Creating command from arguments\n - cores: {threads}\n - configfile: {config} \n - conda-prefix: {conda_dir}\n - "+
        f"snakemake_arguments: {additionals}\n - snakemake_rules: {target_rules}"
    )
    command = f"snakemake --use-conda --cores {threads} --keep-going --configfile {config} --conda-prefix {conda_dir} {additionals} {target_rules}"


    return command




def execute_snakemake(command):
    logger.debug(f"Executing command\n - {command}")
    status = subprocess.Popen(command, shell = True).wait()

    return status


    

def mmaseq(args):

    samplesheet_file = args.samplesheet_file
    input_dir = args.input_dir
    deploy_dir = args.deploy_dir
    outdir = args.outdir
    threads = args.threads
    config = args.config
    # update = args.update
    test = args.test
    debug = args.debug

    logger.debug(
        f"User variables:\n - samplesheet: {samplesheet_file}\n - input_dir: {input_dir}\n - deploy_dir: {deploy_dir}\n - "+
        f"outdir: {outdir}\n - threads: {threads}\n - config: {config}\n - test_active: {test}\n - debug_active: {debug}"
    )

    force = False
    if config is None:
        config = "config/config.yaml"
        force = True
        logger.debug(f"Config file not specified")
    if deploy_dir is None:
        deploy_dir = "Deploy"
        logger.debug(f"Deployment directory not specified.")
    conda_dir = f"{deploy_dir}/conda"

    arguments = []
    rules = []

    if test:
        logger.info("Test run initiated. Will ignore irrelevant user arguments!")
        config = "config/Test.yaml"

        samplesheet_file = "data/samplesheet.tsv"
        outdir = "Test/Results"

        rules.append("all")

        logger.debug(f"Variables before pipeline:\n - config: {config}\n - deploy_dir: {deploy_dir}\n - conda_dir: {conda_dir}\n - force: {force}")

        config = create_config(samplesheet_file = samplesheet_file, outdir = outdir, deploy_dir = deploy_dir, config_file = config, force = True)

        logger.info("Investigating preexisting assemblies")
        config_dir = os.path.dirname(config)
        species_configs_path = f"{config_dir}/species_configs"
        link_assemblies(samplesheet_file = samplesheet_file, config_dir = species_configs_path, outdir = outdir)

        
        command = create_command(threads, config, conda_dir, arguments = arguments, rules = rules)

    else:
        if samplesheet_file is None:
            logger.error("Samplesheet file user argument missing. Aborting!")
            sys.exit(1)

        if outdir is None:
            logger.error("Output directory must be specified. Aborting!")
            sys.exit(1)

        if not os.path.isfile(samplesheet_file) and input_dir is not None:
            logger.info(f"Samplesheet not found, attempting to create at: {samplesheet_file}")
            create_samplesheet(samplesheet_file = samplesheet_file, input_dir = input_dir)
            arguments.append("--dry-run")

        elif not os.path.isfile(samplesheet_file):
            logger.error(f"Samplesheet not found, and input_dir not provided, can't continue. Aborting!")
            sys.exit(1)

        else:
            logger.warning("Samplesheet file exists, but input_dir is also specified. Will ignore input_dir and continue!")

        logger.info("Processing pipeline configurations")
        config = create_config(samplesheet_file = samplesheet_file, outdir = outdir, deploy_dir = deploy_dir, config_file = config, force = force)

        logger.info("Investigating preexisting assemblies")
        config_dir = os.path.dirname(config)
        species_configs_path = f"{config_dir}/species_configs"

        link_assemblies(outdir = outdir, samplesheet_file = samplesheet_file, config_dir = species_configs_path)

        # Communicate defined variables
        logger.debug(f"Variables before pipeline:\n - config: {config}\n - deploy_dir: {deploy_dir}\n - conda_dir: {conda_dir}\n - force: {force}")

        command = create_command(threads, config, conda_dir, arguments = arguments, rules = rules)

    logger.info("Executing pipeline: Mixed Microbial Analysis on Sequencing data")
    status = execute_snakemake(command)
    if status != 0:
        logger.error("Something went wrong while executing snakemake")
    else:
        logger.info("Pipeline successful")


if __name__ == "__main__":
    args = parse_arguments()

    level = logging.INFO
    if args.debug:
        level = logging.DEBUG

    # Configuring logging
    logging.basicConfig(
        level=level,
        format="%(asctime)s-%(levelname)s [%(name)s]: %(message)s"
    )

    # 2) Create a module‑level logger
    logger = logging.getLogger("MMAseq")

    mmaseq(args)
