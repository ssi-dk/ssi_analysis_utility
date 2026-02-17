#!/usr/bin/env python

import argparse
import logging
import os
import re
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd
import yaml

from .helper_functions import determine_sample_configs


# Package location discovery
# ---------------------------------------------------------------------
SRC = Path(__file__).resolve()
ROOT_LIB = SRC.parents[2]
SNAKEFILE = ROOT_LIB / "workflow" / "Snakefile"


# Logging configuration
# ---------------------------------------------------------------------
logger = logging.getLogger("MMAseq")
logger.setLevel(logging.INFO)

_handler = logging.StreamHandler()
_formatter = logging.Formatter("%(levelname)s:%(name)s:%(message)s")
_handler.setFormatter(_formatter)

logger.addHandler(_handler)
logger.propagate = False


# Argument parsing
# ---------------------------------------------------------------------
def parse_arguments():
    parser = argparse.ArgumentParser(
        description = """Configure and execute Seq And Destroy pipeline"""
    )

    parser.add_argument(
        "--samplesheet",
        dest = "samplesheet_file",
        default = None,
        help = """
               Path to samplesheet TSV used by the pipeline (Mandatory). If the 
               samplesheet doesn't exist, `input_dir` must be specified to create one.
            """
    )

    parser.add_argument(
        "--input_dir",
        dest = "input_dir",
        default = None,
        help = """
            Input directory MUST be specified if the samplesheet does not yet exist.
            Input directory will be screened for `.fasta` and `fastq.gz` files, 
            sample_names will be infered from the detected files, and used to populate 
            a samplesheet.After samplesheet creation, the pipeline will be executed in 
            dry-run mode (simulated run)
            """
    )

    parser.add_argument(
        "--deploy_dir",
        dest = "deploy_dir",
        default = ROOT_LIB / "Deploy",
        help = """
            Directory used to deploy databases and conda environments used during 
            pipeline execution. Default is in current working directory.
            If the directory doesn't exist, it will be created during execution.
            To reinstall conda environments remove the folder: {deployment_dir}/conda. 
            To fetch the latest databases remove the folder: {deployment_dir}/Databases
            """
    )

    parser.add_argument(
        "--outdir",
        dest = "outdir",
        default = Path.cwd() / "Results",
        help = """
            Directory used for storing analysis results.
            Default is in current working directory. 
            """
    )

    parser.add_argument(
        "--threads",
        dest = "threads",
        default = 4,
        help = """
               Amount of threads (cores) to dedicate for executing the pipeline.
            """
    )
    
    ## TOCHECK
    parser.add_argument(
        "--config",
        dest = "config",
        default = None,
        help = """
               Configuration file location. (Default ./config/config.yaml)
               If not specified, the config file will be overwritten during subsequent 
               executions.
            """
    )

    parser.add_argument(
        "--test",
        dest = "test",
        action = "store_true",
        help = """
               Perform test run of all modules.\n 
               As a side effect, all conda environments and databases will be created in
               the deployment directory. Results will be stored in Test/Results, and
               config file will be generated in config/Test.yaml 
            """
    )

    parser.add_argument(
        "--debug",
        dest = "debug",
        action = "store_true",
        help = """
               Add debug messages during execution. Mostly used for development and 
               debugging purposes
            """
    )

    return parser.parse_args()



# Path resolution from samplesheet
# ---------------------------------------------------------------------
def resolve_sample_path(path_from_sheet, 
                        samplesheet_file,
                        test=False):
    
    p = Path(path_from_sheet)
    
    # If already absolute → keep it
    if p.is_absolute():
        return p
    
    if test:
        # If test run, resolve relative to root lib
        return (ROOT_LIB / p).resolve()

    # Otherwise resolve relative to samplesheet location
    samplesheet_dir = Path(samplesheet_file).resolve().parent
    return (samplesheet_dir / p).resolve()



# Samplesheet creation
# ---------------------------------------------------------------------
def create_samplesheet(samplesheet_file, 
                       input_dir):
    """
    Will scan the input dir recursively, and register paired end reads and assembly
    files automatically. A samplesheet.tsv will be created at samplesheet_file location,
    it contains sample_name, read1, read2, assembly, and config variables sample_name is
    devised from either {sample_name}_1.fastq.gz or {sample_name}_R1_*.fastq.gz of read1 7
    OR alternately, from {sample_name}.fasta
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
            records.setdefault(sample, {"read1": "NA", 
                                        "read2": "NA", 
                                        "assembly": "NA", 
                                        "config": "default.yaml"})
            records[sample]["assembly"] = str(path)

        elif fname.endswith(".fastq.gz"):
            logger.debug(f"Read file detected: {fname}")

            # remove read indicators to get sample name
            sample = re.sub(r'(_R?[12].*)\.fastq\.gz$', '', fname)
            records.setdefault(sample, {"read1": "NA", 
                                        "read2": "NA", 
                                        "assembly": "NA", 
                                        "config": "default.yaml"})

            if re.search(r'(_R?1|_1)', fname):
                records[sample]["read1"] = str(path)
            elif re.search(r'(_R?2|_2)', fname):
                records[sample]["read2"] = str(path)

    samplesheet = (
        pd.DataFrame.from_dict(records, 
                               orient="index")
        .reset_index()
        .rename(columns={"index": "sample_name"})
    )

    samplesheet.to_csv(samplesheet_file, 
                       sep = "\t", 
                       index = False)

    return None



# Configuration file creation
# ---------------------------------------------------------------------
def create_config(samplesheet_file, 
                  outdir, 
                  deploy_dir, 
                  config_file, 
                  root, 
                  force):
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
        logger.warning(
        (
            "Configuration file already exists, will change current config_file to:\n"
            f"- {config_file}"
        )
    )
    
    # Cleanup outdir
    outdir = os.path.abspath(outdir).rstrip("/")
    
    config = {
        "samplesheet": str(samplesheet_file),
        "deploy_dir": str(deploy_dir),
        "outdir": str(outdir),
        "root": str(root)
    }

    with open(config_file, 'w') as config_yaml:
        yaml.safe_dump(config, 
                       config_yaml)

    return config_file



# Create assembly symlinks
#----------------------------------------------------------------------
def link_assemblies(samplesheet_file,
                    config_dir,
                    outdir,
                    test):
    """
    Generate assembly symlinks for each sample and assembler type.
    """

    logger.debug("Reading samplesheet")
    samplesheet = pd.read_csv(samplesheet_file, 
                              sep="\t").set_index("sample_name")
    

    logger.debug("Importing sample configs")
    sample_configs = determine_sample_configs(samplesheet=samplesheet,
                                              config_dir=config_dir)

    outdir = Path(outdir)  ### TO RECHECK

    logger.debug("Initiating symlink generation")

    for sample, configs in sample_configs.items():

        assembly_sheet = samplesheet.at[sample, "assembly"]
        assembly_source = resolve_sample_path(assembly_sheet,
                                              samplesheet_file,
                                              test=test)

        if not assembly_source.is_file():
            if assembly_sheet.lower() not in ["na", ""]:
                logger.warning(
                    "Assembly file not found: %s", assembly_source
                )
            continue


        # Extract ALL assembler types once per sample
        assemblers = {
            assembler
            for options in configs.values()
            if isinstance(options, dict)
            and "assemblers" in options
            for assembler in (
                options["assemblers"]
                if isinstance(options["assemblers"], list)
                else [options["assemblers"]]
            )
        }

        # Create symlinks (single loop now)
        for assembler in assemblers:

            assembly_dir = outdir / sample / assembler
            assembly_dir.mkdir(parents=True, 
                               exist_ok=True)

            destination = assembly_dir / f"{sample}.fasta"

            if destination.is_symlink():
                logger.debug("Removing existing symlink: %s", destination)
                destination.unlink()

            if not destination.exists():
                logger.info(
                    "Creating symlink: %s -> %s",
                    assembly_source,
                    destination,
                )
                destination.symlink_to(assembly_source)
            else:
                logger.warning(
                    "File exists and is not a symlink: %s",
                    destination,
                )



# Command creation and execution
# ---------------------------------------------------------------------
def create_command(threads, 
                   config, 
                   conda_dir, 
                   arguments = None, 
                   rules = None):
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
        (
            "Creating command from arguments\n"
            f" - cores: {threads}\n"
            f" - configfile: {config}\n"
            f" - conda-prefix: {conda_dir}\n"
            f" - snakemake_arguments: {additionals}\n"
            f" - snakemake_rules: {target_rules}"
        )
    )
    
    command = (
        "snakemake --use-conda "
        f"--cores {threads} "
        "--keep-going "
        f"--configfile {config} "
        f" --snakefile {SNAKEFILE} "
        f"--conda-prefix {conda_dir} "
                       f"{additionals} "
                       f"{target_rules}"
        )

    return command



# Execute snakemake command
#---------------------------------------------------------------------
def execute_snakemake(command):
    
    logger.debug(f"Executing command\n - {command}")
    status = subprocess.Popen(command, shell = True).wait()

    return status
    


# Main function
#----------------------------------------------------------------------
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
        (
            "User variables:\n"
            f" - samplesheet: {samplesheet_file}\n"
            f" - input_dir: {input_dir}\n"
            f" - threads: {threads}\n "
            f" - config: {config}\n"
            f" - test_active: {test}\n"
            f"- debug_active: {debug}"
        )
    )

    force = False
    if config is None:
        config = ROOT_LIB / "config/config.yaml" # Revert to default config location if not specified
        force = True
        logger.debug(f"Config file not specified")
    
    # if deploy_dir is None:
    #     deploy_dir = ROOT_LIB / "Deploy" # Default deploy location if not specified 
    #     logger.debug(f"Deployment directory not specified.")
    conda_dir = f"{deploy_dir}/conda"

    arguments = []
    rules = []

    if test:
        logger.info("Test run initiated. Will ignore irrelevant user arguments!")
        config = f"{ROOT_LIB}/config/Test.yaml" # Default test config location
        samplesheet_file = f"{ROOT_LIB}/data/samplesheet.tsv" # Default test samplesheet location
        outdir = f"{ROOT_LIB}/Test/Results" # Default test output location
        rules.append("all")


        logger.info(
            (
                "Variables before pipeline:\n"
                f" - config: {config}\n"
                f" - deploy_dir: {deploy_dir}\n"
                f" - conda_dir: {conda_dir}\n"
                f" - force: {force}"
                f" --snakefile {SNAKEFILE}"
                
            )
        )
        config = create_config(samplesheet_file = samplesheet_file, 
                               outdir = outdir, 
                               deploy_dir = deploy_dir, 
                               config_file = config, 
                               force = True, 
                               root = ROOT_LIB)

        logger.info("Investigating preexisting assemblies")
        species_configs_path = f"{ROOT_LIB}/config/species_configs"
        link_assemblies(samplesheet_file = samplesheet_file, 
                        config_dir = species_configs_path, 
                        outdir = outdir,
                        test=True)

        
        command = create_command(threads, 
                                 config, 
                                 conda_dir, 
                                 arguments = arguments, 
                                 rules = rules)

    else:
        if samplesheet_file is None:
            logger.error("Samplesheet file user argument missing. Aborting!")
            sys.exit(1)

        # if outdir is None:
        #     logger.error("Output directory must be specified. Aborting!")
        #     sys.exit(1)

        if not os.path.isfile(samplesheet_file):# and input_dir is not None:
            logger.info(f"Samplesheet not found, attempting to create at: {samplesheet_file}")
            create_samplesheet(samplesheet_file = samplesheet_file, 
                               input_dir = input_dir)
            arguments.append("--dry-run")

        # elif not os.path.isfile(samplesheet_file):
        #     logger.error(f"Samplesheet not found, and input_dir not provided, can't continue. Aborting!")
        #     sys.exit(1)

        else:
            logger.warning("Samplesheet file exists and input_dir is also specified. Will ignore input_dir and continue!")

        logger.info("Processing pipeline configurations")
        config = create_config(samplesheet_file = samplesheet_file, 
                               outdir = outdir, 
                               deploy_dir = deploy_dir, 
                               config_file = config, 
                               force = force, 
                               root = ROOT_LIB)

        logger.info("Investigating preexisting assemblies")
        species_configs_path =  f"{ROOT_LIB}/config/species_configs"
        link_assemblies(samplesheet_file = samplesheet_file, 
                        outdir = outdir, 
                        config_dir = species_configs_path,
                        test=False)

        # Communicate defined variables        
        logger.debug(
            (
                "Variables before pipeline:\n"
                f" - config: {config}\n"
                f" - deploy_dir: {deploy_dir}\n"
                f" - conda_dir: {conda_dir}\n"
                f" - force: {force}"
                f"--snakefile {SNAKEFILE}"
                
            )
        )

        command = create_command(threads, 
                                 config, 
                                 conda_dir, 
                                 arguments = arguments, 
                                 rules = rules)

    logger.info("Executing pipeline: Mixed Microbial Analysis on Sequencing data")
    status = execute_snakemake(command)
    if status != 0:
        logger.error("Something went wrong while executing snakemake")
    else:
        logger.info("Pipeline successful")



# Launcher
#---------------------------------------------------------------------
def launcher() -> None:
    args = parse_arguments()

    # Determine level of logging
    if args.debug:
        level = logging.DEBUG
        logger.setLevel(logging.INFO)

    mmaseq(args)
