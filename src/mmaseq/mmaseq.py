#!/usr/bin/env python
from .utils import pkg_logging, determine_sample_configs, parse_mmaseq
from .utils import *
from pathlib import Path
import pandas as pd
import os
import re
import subprocess
import sys
from datetime import datetime
import yaml

print(SPE_CONFIGS)

# Initiate logging
logger = pkg_logging.initiate_log("MMAseq")

def resolve_path(path, cwd, samplesheet_file):

    # Determine possible file paths
    path = Path(path)
    path_absolute = absolute.is_absolute()
    path_from_cwd = cwd / path
    samplesheet_dir = Path(samplesheet_file).resolve().parent
    path_from_samplesheet_dir = samplesheet_dir / path

    # Hierichically look through potential file paths
    if path_absolute:
        absolute = path
    elif not path_absolute & path_from_cwd.exists():
        absolute = path_from_cwd.resolve()
    elif not path_absolute & path_from_samplesheet_dir.exists(): # I suggest removing this, as default behavoir is files relative to cwd NOT input file
        absolute = path_from_samplesheet_dir.resolve()
    else:
        # Abort if file path was unsolvable
        logger.error(f"""Unable to resolve sample path: {path}.
             - Resolve the paths in your samplesheet and try again. Aborting!""")
        sys.exit(1)

    return absolute



def resolve_samplesheet_paths(samplesheet_file):
    """
    Converts all read1/read2/assembly paths in the samplesheet
    into absolute paths using resolve_sample_path().
    This guarantees that Snakemake always receives valid paths.
    """
    def fix(path):
        # Ensure that path is specified
        if isinstance(path, str) and path.lower() not in ["na", ""]:

            # Resolve path to absolute
            path = resolve_path(path, CWD, samplesheet_file)

        return path

    samplesheet = pd.read_csv(samplesheet_file, sep="\t").copy()

    # Resolve read and assembly file paths to absolute
    samplesheet["read1"] = samplesheet["read1"].apply(fix)
    samplesheet["read2"] = samplesheet["read2"].apply(fix)
    samplesheet["assembly"] = samplesheet["assembly"].apply(fix)

    samplesheet.to_csv(samplesheet_file,
        sep = "\t", 
        index = False)

    return True


def create_config(samplesheet_file, 
                  outdir, 
                  deploy_dir, 
                  config_file):

    # Check whether config file exists
    if not os.path.isfile(config_file):
        config_dir = os.path.dirname(config_file)
        if not os.path.isdir(config_dir):
            logging.info("Config directory does not exist, creating directory")
            os.makedirs(config_dir)
    else:
        # Ensure config is not overwritten
        timestamp = datetime.now().strftime("%y_%m_%d-%H_%M")
        config_file = f"{os.path.dirname(config_file)}/{timestamp}_{os.path.basename(config_file)}"
        logger.warning(
            (
                "Configuration file already exists, will change current config_file to:\n"
                f"- {config_file}"
            )
        )
    
    # Record config contents
    outdir = os.path.abspath(outdir).rstrip("/")
    config = {
        "samplesheet": str(samplesheet_file),
        "deploy_dir": str(deploy_dir),
        "outdir": str(outdir)
    }

    logger.debug(f"Creating config file at {config_file}")
    with open(config_file, "w") as config_yaml:
        yaml.safe_dump(config, config_yaml)

    # Return the path of newly create config file
    return config_file


def link_assemblies(samplesheet_file, 
                    config_dir, 
                    outdir):

    logger.debug(f"Reading samplesheet from {samplesheet_file}")
    samplesheet = pd.read_csv(samplesheet_file, 
                              sep = "\t").set_index("sample_name")


    logger.debug("Importing sample configs")
    sample_configs = determine_sample_configs(samplesheet = samplesheet, 
                                              config_dir = config_dir)

    # Iterating over sample configurations
    for sample, configs in sample_configs.items():

        # Reading assembly entry from samplesheet
        assembly_sheet = samplesheet.at[sample, "assembly"]
        assembly_source = Path(assembly_sheet)

        logger.debug(f"Assembly for {sample} in samplehseet is {assembly_sheet}")

        if pd.isna(assembly_sheet):
            # Handle no determined assembly sheet
            logger.debug(f"No assembly provided for {sample} — skipping!")

        elif not assembly_source.exists():
            # Handle assembly path is invalid or not existing
            logger.warning(f"""Assembly path does not exist for {sample} 
                at location {assembly_sheet}\nSkipping!""")
        else:
            logger.debug(f"Assembly found at {assembly_source}")

            # Determine assemblers specified in sample configs
            assemblers = {
                assembler ### I DON'T GET THIS LOOP...
                for options in configs.values()
                if isinstance(options, dict)
                and "assemblers" in options
                for assembler in (
                    options["assemblers"]
                    if isinstance(options["assemblers"], list)
                    else [options["assemblers"]]
                )
            }

            # Define assembly type from sample configurations
            for assembler in assemblers:
                assembly_dir = outdir / sample / assembler

                # Ensure output assembly results directory exists
                if not assembly_dir.exists():
                    assembly_dir.mkdir(parents = True)

                destination = assembly_dir / f"{sample}.fasta"

                # Handle destination when being broken links
                if destination.is_symlink() and not destination.exists(follow_symlinks = True):
                    logger.warning(f"""Assembly results directory is a 
                        broken link -> unlinking: {destination}""")
                    destination.unlink()

                # Ignore prexisting functional symbolic links at destination
                if destination.is_symlink():
                    logger.info(f"""Assembly allready linked to results
                     directory at {destination}.\nSkipping!""")

                # Initiate symlink creation if destination is empty
                if not destination.exists(follow_symlinks = False):
                    logger.info(f"""Creating symlink: 
                        {assembly_source} -> {destination}""")
                    destination.symlink_to(assembly_source)

                # Warn if assembly file allready exists, but not as a link
                else:
                    logger.warning(f"""File exists and is not a symlink {destination}.
                        Skipping!""")


def create_command(threads, 
                   config, 
                   conda_dir, 
                   arguments = None, 
                   rules = None):

    # Define arguments and rules as a single string
    additionals = " ".join(arguments) if arguments else ""
    target_rules = " ".join(rules) if rules else ""

    # Determine command
    command = (
        "snakemake --use-conda "
        f"--cores {threads} "
        "--keep-going "
        f"--configfile {config} "
        f"--snakefile {SNAKEFILE} "
        f"--conda-prefix {conda_dir} "
        f"{additionals} "
        f"{target_rules}"
    )

    logger.debug(f"Pipeline command created: {command}")

    return command


def execute_snakemake(command):
    status = subprocess.Popen(command, 
                              shell = True).wait()
    return status


def mmaseq(args):

    # Define user input
    samplesheet_file = Path(args.samplesheet_file)
    deploy_dir = Path(args.deploy_dir)
    outdir = Path(args.outdir)
    threads = args.threads
    config = args.config
    resolve = args.resolve
    ignore_assemblies = args.ignore_assemblies

    # Resolve other objects
    conda_dir = (deploy_dir / "conda").resolve()
    arguments = []
    rules = []

    # Handle missing samplesheet
    if not os.path.isfile(samplesheet_file):
        logger.error("""Samplesheet not found.\n - Execute the create module 
            first.\nAborting!""")
        sys.exit(1)

    # Enable path normalization
    if resolve:
        logger.info("Normalizing paths")
        resolve_samplesheet_paths(samplesheet_file)

    logger.info("Creating pipeline configurations")
    config = create_config(samplesheet_file, 
                           outdir, 
                           deploy_dir, 
                           config)

    if not ignore_assemblies:
        logger.info("Creating symbolic links for assemblies")
        link_assemblies(samplesheet_file, SPE_CONFIGS, outdir)

    logger.info("Creating pipeline command")
    command = create_command(threads, 
                             config, 
                             conda_dir, 
                             arguments, 
                             rules)

    logger.info("Executing pipeline: Mixed Microbial Analysis on Sequencing data")
    status = execute_snakemake(command)

    if status != 0:
        logger.error("Something went wrong while executing snakemake.")
    else:
        logger.info("Success!")


def launcher() -> None:
    # Parse user input
    args = parse_mmaseq()

    # Generate logger
    pkg_logging.adjust_log_level(logger, args.debug)

    logger.info("Initiating MMAseq")
    mmaseq(args)

    logger.info("Exitting")