#!/usr/bin/env python
from .utils import pkg_logging
from .utils.paths import *
import argparse
from pathlib import Path
import sys
import re
import pandas as pd
from datetime import datetime
import yaml
from .utils import helper_functions
import subprocess

# Initiate logging
logger = pkg_logging.initiate_log("MMAseq")

def parse_mmaseq():
    parser = argparse.ArgumentParser(
        description = """
        Execute Mixed Microbial Analysis on Sequencing data (MMAseq)
        """
    )

    parser.add_argument(
        "--samplesheet",
        dest = "samplesheet_file",
        required = True,
        help = """
            Path to samplesheet TSV used by the pipeline. 
            If the samplesheet doesn't exist, `indir` must be 
            specified to create one.
        """
    )
    

    parser.add_argument(
        "--deploy_dir",
        dest = "deploy_dir",
        default = PKG_DIR / "Deploy",
        help = f"""
            Directory used to deploy virtual environment and databases 
            used during pipeline execution. To reinstall environments 
            and/or databases, remove the `conda/` and/or the `Databases/` 
            folders in the deployment directory. (Default {PKG_DIR}/Deploy)
        """
    )

    parser.add_argument(
        "--outdir",
        dest = "outdir",
        default = CWD / "MMAseq_Results",
        help = """
            Directory used for storing analysis results. (Default ./MMAseq_Results)
        """
    )

    parser.add_argument(
        "--config",
        dest = "config",
        default = PKG_CONFIGS / "config.yaml",
        help = f"""
            Configuration file location. (Default {PKG_CONFIGS}/config.yaml)
        """
    )

    parser.add_argument(
        "--threads",
        dest = "threads",
        default = 4,
        help = """
            Amount of threads (cores) to dedicate for executing the pipeline. 
            (Default 4)
        """
    )

    parser.add_argument(
        "--resolve",
        dest = "resolve",
        action = "store_true",
        help = """
            Resolves absolute paths from samplesheet columns, will 
            overwrite samplesheet. (Default False)
        """
    )

    parser.add_argument(
        "--ignore_assemblies",
        dest = "ignore_assemblies",
        action = "store_true",
        help = """
            Avoid creating symbolic links of the assemblies into the 
            pipeline output directory. (Default False) 
            If this is not specified, symbolic links will be created 
            to the output directory, hence avoidubg assembly steps in the 
            pipeline. Use this option to enforce the pipeline to create 
            assemblies (Will take extra time) rather than relying on those 
            specified in the samplesheet.
        """
    )

    parser.add_argument(
        "--debug",
        dest = "debug",
        action = "store_true",
        help = """
            Add debug messages during execution. (Default False) 
            Mostly used for development and debugging purposes.
        """
    )

    return parser.parse_args()



def resolve_path(path, samplesheet_file):

    # Determine possible file paths
    path_absolute = path.is_absolute()
    path_from_cwd = CWD / path
    samplesheet_dir = samplesheet_file.resolve().parent
    path_from_samplesheet_dir = samplesheet_dir / path

    # Hierichically look through potential file paths
    if path_absolute & path.exists():
        absolute = path
        logger.debug(f"Path is absolute, no resolving needed {absolute}")
    elif (not path_absolute) & path_from_cwd.exists():
        absolute = path_from_cwd.resolve()
        logger.debug(f"Path not absolute, but found relative to current location {absolute}")
    elif (not path_absolute) & path_from_samplesheet_dir.exists():
        absolute = path_from_samplesheet_dir.resolve()
        logger.debug(f"Path not absolute, but found relative to samplesheet directory {absolute}")

    else:
        # Abort if file path was unsolvable
        logger.error((
            f"Unable to resolve sample path: {path}.\n"
            "- Resolve the paths in your samplesheet and try again. Aborting!"
        ))
        sys.exit(1)

    return absolute


def resolve_samplesheet_paths(samplesheet_file, outdir):
    """
    Converts all read1/read2/assembly paths in the samplesheet
    into absolute paths using resolve_sample_path().
    This guarantees that Snakemake always receives valid paths.
    """
    def fix(path):
        # Ensure that path is specified
        if isinstance(path, str) and path.lower() not in ["na", ""]:

            # Resolve path to absolute
            path = resolve_path(Path(path), samplesheet_file)

        return path

    samplesheet = pd.read_csv(samplesheet_file, sep="\t").copy()


    # Resolve read and assembly file paths to absolute
    samplesheet["read1"] = samplesheet["read1"].apply(fix)
    samplesheet["read2"] = samplesheet["read2"].apply(fix)
    samplesheet["assembly"] = samplesheet["assembly"].apply(fix)

    # Write samplesheet with resolved paths
    if not outdir.exists():
        outdir.mkdir(parents = True)

    samplesheet_resolved_file = outdir / re.sub(
        ".tsv", "_resolved.tsv", str(samplesheet_file.name)
    )
    logger.debug((
        "Writting samplesheet with resolved paths to "
        f"{samplesheet_resolved_file}"
    ))
    samplesheet.to_csv(samplesheet_resolved_file,
        sep = "\t", 
        index = False)

    return samplesheet_resolved_file


def create_config(samplesheet_file, 
                  outdir, 
                  deploy_dir
                  ):

    # Determine config file
    outdir = outdir.resolve()
    config_file = outdir / "config.yaml"

    # Check whether config file exists
    if not config_file.exists():
        logger.debug(f"Creating new config file {config_file}")
        if not outdir.exists():
            logger.debug((
                "Output directory does not exists, "
                f"creating directory {outdir}"
            ))

            outdir.mkdir(parents = True)
    else:
        # Ensure config is not overwritten
        timestamp = datetime.now().strftime("%y_%m_%d-%H_%M")
        config_file = config_file.parent / f"{timestamp}_{config_file.name}"
        logger.warning((
            "Old configuration file detected, "
            f"creating new file {config_file}"
        ))
    
    # Record config contents
    config = {
        "samplesheet": str(samplesheet_file),
        "deploy_dir": str(deploy_dir),
        "outdir": str(outdir)
    }

    logger.debug(f"Creating config file {config_file}")
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


    logger.debug(f"Importing sample configs from {config_dir}")
    sample_configs = helper_functions.determine_sample_configs(samplesheet = samplesheet, 
                                              config_dir = config_dir)

    # Iterating over sample configurations
    for sample, configs in sample_configs.items():

        # Reading assembly entry from samplesheet
        assembly_from_sheet = samplesheet.at[sample, "assembly"]
        
        # Handle if assembly is determined as NA
        if pd.isna(assembly_from_sheet):
            logger.debug(f"No assembly provided for {sample} — skipping!")
            continue

        assembly_source = Path(assembly_from_sheet)
        logger.debug(f"Assembly for {sample} in samplehseet is {assembly_source}")

        # Attempt to locate relative paths from assembly listed in sheet
        if assembly_source.exists(follow_symlinks = True):
            # Handle if assembly file exists with a valid path
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



        else:

            # Handle if assembly file path is invalid or not existing
            logger.warning(f"""Assembly path does not exist for {sample} 
                at location {assembly_source}\nSkipping!""")


def create_command(threads, 
                   config_file, 
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
        f"--configfile {config_file} "
        f"--snakefile {SNAKEFILE} "
        f"--conda-prefix {conda_dir} "
        f"{additionals} "
        f"{target_rules}"
    )

    logger.debug(f"Pipeline command created: {command}")

    return command


def execute_snakemake(command):
    status = subprocess.Popen(command, shell = True).wait()
    return status


def mmaseq(args):

    # Define user input
    samplesheet_file = Path(args.samplesheet_file)
    deploy_dir = Path(args.deploy_dir)
    outdir = Path(args.outdir)
    threads = args.threads
    resolve = args.resolve
    ignore_assemblies = args.ignore_assemblies

    # Resolve other objects
    conda_dir = (deploy_dir / "conda").resolve()
    arguments = []
    rules = []

    # Handle missing samplesheet
    if not samplesheet_file.exists():
        logger.error("""Samplesheet not found.\n - Execute the create module 
            first.\nAborting!""")
        sys.exit(1)

    # Enable path normalization
    if resolve:
        logger.info("Creating new samplesheet with resolved paths")
        samplesheet_file = resolve_samplesheet_paths(samplesheet_file, outdir)

    logger.info("Creating configurations")
    config_file = create_config(samplesheet_file, 
                           outdir, 
                           deploy_dir
                           )

    if ignore_assemblies:
        logger.info("Will not skip assembly parts with assemblies from samplehseet")
    else:
        logger.info("Creating symbolic links for assemblies")
        link_assemblies(samplesheet_file, str(SPE_CONFIGS), outdir)


    logger.info("Creating pipeline command")
    command = create_command(threads, 
                             config_file, 
                             conda_dir, 
                             arguments, 
                             rules)

    logger.info("Executing pipeline: Mixed Microbial Analysis on Sequencing data")
    status = execute_snakemake(command)

    if status != 0:
        logger.error("Something went wrong while executing snakemake.")
        sys.exit(1)


def launcher() -> None:
    # Parse user input
    args = parse_mmaseq()

    # Generate logger
    pkg_logging.adjust_log_level(logger, args.debug)

    logger.info("Initiating MMAseq")
    mmaseq(args)

    logger.info("MMAseq successful!")