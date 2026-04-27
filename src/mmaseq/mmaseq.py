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
    # Define command line arguments
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=(
        f"Mixed Microbial Analysis on Sequencing data\n"
        f"This is the main module used for executing the MMAseq pipeline.\n"
        f"MMAseq is split into the following three executable commands:\n"
        f"    * mmacreate - Create input samplesheets for MMAseq\n"
        f"    * mmadeploy - Install environment and create databases\n"
        f"    * mmaseq - Execute the pipeline.\n"
        )
        
    )

    parser.add_argument(
        "--samplesheet",
        dest = "samplesheet_file",
        required = True,
        help = 
           """
            Path to samplesheet TSV used by the pipeline. 
            If the samplesheet doesn't exist, `indir` must be 
            specified to create one.
            """
    )
    

    parser.add_argument(
        "--deploy_dir",
        dest = "deploy_dir",
        default = PKG_DIR / "Deploy",
        help = 
            f"""
            Directory used to deploy virtual environment and databases 
            used during pipeline execution. To reinstall environments 
            and/or databases, remove the `conda/` and/or the `Databases/` 
            folders in the deployment directory. (Default: %(default)s)
            """
    )

    parser.add_argument(
        "--outdir",
        dest = "outdir",
        required = True,
        help = 
            """
            Directory used for storing analysis results.
            """
    )

    parser.add_argument(
        "--config",
        dest = "config",
        default = PKG_CONFIGS / "config.yaml",
        help = 
            f"""
            Configuration file location. (Default: %(default)s)
            """
    )

    parser.add_argument(
        "--threads",
        dest = "threads",
        default = 4,
        help = 
            f"""
            Amount of threads (cores) to dedicate for executing the pipeline. 
            (Default: %(default)s)
            """
    )

    parser.add_argument(
        "--resolve",
        dest = "resolve",
        action = "store_true",
        help = 
            f"""
            Resolves absolute paths from samplesheet columns, will 
            overwrite samplesheet. (Default: %(default)s)
            """
    )

    parser.add_argument(
        "--force",
        dest = "force",
        action = "store_true",
        help = 
            """
            Force running all rules. (Default: %(default)s) 
            This will cause all necessary rules for a given run to be executed, 
            which e.g. run assemblies despite preexisting ones existing.
            Mostly useful for forcing deployment, or rerunning a suspicious batch.
            """
    )

    parser.add_argument(
        "--ignore_assemblies",
        dest = "ignore_assemblies",
        action = "store_true",
        help = 
            """
            Avoid creating symbolic links of the assemblies into the 
            pipeline output directory. (Default: %(default)s) 
            If this is not specified, symbolic links will be created 
            to the output directory, hence avoiding assembly steps in the 
            pipeline. Use this option to enforce the pipeline to create 
            assemblies (Will take extra time) rather than relying on those 
            specified in the samplesheet.
            """
    )

    parser.add_argument(
        "--verbosity",
        dest = "verbosity",
        type = int,
        choices = [0, 1, 2],
        default = 0,
        help = 
            """
            Adjust the verbosity (Default: %(default)s); 
            0: Minimal messages, 
            1: Debug messages, 
            2: Trace messages (development only)
            """
    )

    parser.add_argument(
        "--logfile",
        dest = "logfile",
        type = str,
        default = None,
        help = 
            f"""
            If provided, will redirect log messages from STDOUT to logfile. (Default: %(default)s) 
            Will be ignored if logfile parent folder doesn't exists.
            """
    )

    return parser.parse_args()


def resolve_path(path, samplesheet_file):
    logger.trace((
                  "resolve_path(\n - "
                  f"path: {path}\n - "
                  f"samplesheet_file: {samplesheet_file})"))

    # Determine possible file paths
    path_absolute = path.is_absolute()
    path_from_cwd = CWD / path
    path_from_samplesheet_dir = samplesheet_file.resolve().parent / path

    # Hierichically look through potential file paths
    if path_absolute & path.exists():
        absolute = path
        logger.trace(
            f"Path is absolute, no resolving needed:\n - {absolute}"
        )
    elif (not path_absolute) & path_from_cwd.exists():
        absolute = path_from_cwd.resolve()
        logger.trace((
            f"Path is not absolute, but found relative to current location:"
            f"\n - {absolute}"
        ))
    elif (not path_absolute) & path_from_samplesheet_dir.exists():
        absolute = path_from_samplesheet_dir.resolve()
        logger.trace((
            f"Path is not absolute, but found relative to samplesheet:"
            f"\n - {absolute}"
        ))

    else:
        # Abort if file path was unsolvable
        logger.error((
            f"Unable to resolve sample path: {path}.\n"
            f"- Resolve the paths in your samplesheet {samplesheet_file} "
            "and try again.\nAborting!"
        ))
        sys.exit(1)

    return absolute


def resolve_samplesheet_paths(samplesheet_file, outdir):
    """
    Converts all read1/read2/assembly paths in the samplesheet
    into absolute paths using resolve_path().
    This guarantees that Snakemake always receives valid paths.
    """
    logger.trace((
            "resolve_samplesheet_paths(\n - "
            f"samplesheet_file: {samplesheet_file}\n - "
            f"outdir: {outdir})"))

    def fix(path):
        # Ensure that path is specified
        if isinstance(path, str) and path.lower() not in ["na", ""]:

            # Resolve path to absolute
            path = resolve_path(Path(path), samplesheet_file)

        return path


    samplesheet = pd.read_csv(samplesheet_file, sep="\t").copy()


    # Resolve read and assembly file paths
    samplesheet["read1"] = samplesheet["read1"].apply(fix)
    samplesheet["read2"] = samplesheet["read2"].apply(fix)
    samplesheet["assembly"] = samplesheet["assembly"].apply(fix)

    # Write samplesheet with resolved paths"
    samplesheet_resolved_file = outdir / re.sub(
        ".tsv", "_resolved.tsv", str(samplesheet_file.name)
    )

    logger.debug(("Writing samplesheet with resolved paths "
                  f"{samplesheet_resolved_file}"))

    if not outdir.exists():
        outdir.mkdir(parents = True)
    samplesheet.to_csv(samplesheet_resolved_file,
        sep = "\t", 
        index = False)

    return samplesheet_resolved_file


def create_config(samplesheet_file, 
                  outdir, 
                  deploy_dir,
                  verbosity
                  ):

    logger.trace(("create_config(\n - "
                  f"samplesheet_file: {samplesheet_file}\n - "
                  f"outdir: {outdir}\n - "
                  f"deploy_dir: {deploy_dir}\n - "
                  f"verbosity: {verbosity}"))

    # Determine config file
    outdir = outdir.resolve()
    config_file = outdir / "config.yaml"

    # Check whether config file exists
    if not config_file.exists():
        logger.trace(f"Creating new config file {config_file}")
        if not outdir.exists():
            logger.trace((
                "Output directory does not exist, "
                f"creating directory {outdir}"
            ))

            outdir.mkdir(parents = True)
    else:
        # Ensure config is not overwritten
        timestamp = datetime.now().strftime("%y_%m_%d-%H_%M")
        config_file = config_file.parent / f"{timestamp}_{config_file.name}"
        logger.debug("Old configuration file detected. Creating new file!")
    
    # Record config contents
    config = {
        "samplesheet": str(samplesheet_file),
        "deploy_dir": str(deploy_dir),
        "outdir": str(outdir),
        "verbosity": int(verbosity)
    }

    logger.debug(f"Writing config file to {config_file}")
    with open(config_file, "w") as config_yaml:
        yaml.safe_dump(config, config_yaml)

    # Return the path of newly create config file
    return config_file


def link_assemblies(samplesheet_file, 
                    config_dir, 
                    outdir):
    logger.trace((
        "link_assemblies(\n - "
        f"samplesheet_file: {samplesheet_file}\n - "
        f"config_dir: {config_dir}\n - "
        f"outdir: {outdir})"))

    logger.debug("Initiating assembly symlinking")
    logger.trace(f"Reading samplesheet from {samplesheet_file}")
    samplesheet = pd.read_csv(samplesheet_file, 
                              sep = "\t").set_index("sample_name")


    logger.trace(f"Importing sample configs from {config_dir}")
    sample_configs = helper_functions.determine_sample_configs(samplesheet = samplesheet, 
                                                               config_dir = config_dir)

    # Iterating over sample configurations
    for sample, configs in sample_configs.items():

        # Reading assembly entry from samplesheet
        assembly_from_sheet = samplesheet.at[sample, "assembly"]
        
        # Handle if assembly is determined as NA
        if pd.isna(assembly_from_sheet):
            logger.trace(f"No assembly provided for {sample} — skipping!")
            continue

        assembly_source = Path(assembly_from_sheet)
        logger.trace(f"Assembly for {sample} in samplesheet is {assembly_source}")

        # Attempt to locate relative paths from assembly listed in sheet
        if assembly_source.exists(follow_symlinks = True):
            # Handle if assembly file exists with a valid path
            logger.trace(f"Assembly found at {assembly_source}")

            # Determine assemblers specified in sample configs
            assemblers = set() 
            for options in configs.values():
                if not isinstance(options, dict) or "assemblers" not in options:
                    continue
                raw = options["assemblers"]
                assembler_list = raw if isinstance(raw, list) else [raw]
                assemblers.update(assembler_list)

            # Define assembly type from sample configurations
            for assembler in assemblers:
                assembly_dir = outdir / sample / assembler
                destination = assembly_dir / f"{sample}.fasta"

                # Ensure output assembly results directory exists
                if not assembly_dir.exists():
                    assembly_dir.mkdir(parents = True)

                # Handle destination when being broken links
                if destination.is_symlink() and not destination.exists(follow_symlinks = True):
                    logger.warning(f"""Assembly results directory is a 
                        broken link -> unlinking: {destination}""")
                    destination.unlink()

                # Ignore pre-existing functional symbolic links at destination
                if destination.is_symlink():
                    logger.debug(f"""Assembly already linked to results
                     directory at {destination}.\nSkipping!""")

                # Initiate symlink creation if destination is empty
                if not destination.exists(follow_symlinks = False):
                    logger.debug(f"""Creating symlink: 
                        {assembly_source} -> {destination}""")
                    destination.symlink_to(assembly_source)

                # Warn if assembly file already exists, but not as a link
                else:
                    logger.warning(f"""File exists and is not a symlink {destination}.
                        Skipping!""")



        else:

            # Handle if assembly file path is invalid or not existing
            logger.warning(f"Assembly path does not exist for {sample}"
                           f"at location {assembly_source}\nSkipping!")


def create_command(threads, 
                   config_file, 
                   conda_dir,
                   force,
                   arguments = None, 
                   rules = None):
    logger.trace(("create_command(\n - "
                  f"config_file: {config_file}\n - "
                  f"conda_dir: {conda_dir}\n - "
                  f"force: {force}\n - "
                  f"arguments: {arguments}\n - "
                  f"rules: {rules})"))

    # Define arguments and rules as a single string
    additionals = " ".join(arguments) if arguments else ""
    target_rules = " ".join(rules) if rules else ""

    if force:
        additionals += " --forceall"

    # Determine command
    command = (
        "snakemake --use-conda "
        f"--cores {threads} "
        "--keep-going "
        "--rerun-incomplete "
        f"--configfile {config_file} "
        f"--snakefile {SNAKEFILE} "
        f"--conda-prefix {conda_dir} "
        f"{additionals} "
        f"{target_rules}"
    )

    return command


def execute_snakemake(command):
    logger.debug(f"Executing:\n{command}")
    status = subprocess.Popen(command, shell = True).wait()
    return status


def mmaseq(args):

    # Define user input
    samplesheet_file = Path(args.samplesheet_file)
    deploy_dir = Path(args.deploy_dir)
    outdir = Path(args.outdir)
    threads = args.threads
    resolve = args.resolve
    force = args.force
    ignore_assemblies = args.ignore_assemblies

    # Resolve other objects
    conda_dir = (deploy_dir / "conda").resolve()
    logger.trace((f"User arguments and deduced variables:\n - "
                  f"samplesheet_file: {samplesheet_file}\n - "
                  f"deploy_dir: {deploy_dir}\n - "
                  f"outdir: {outdir}\n - "
                  f"threads: {threads}\n - "
                  f"resolve: {resolve}\n - "
                  f"ignore_assemblies: {ignore_assemblies}\n - "
                  f"conda_dir: {conda_dir}"))

    # Handle missing samplesheet
    if not samplesheet_file.exists():
        logger.error(f"Samplesheet not found.\n "
                     f"Execute the create module first\n"
                     f"Aborting!")
        sys.exit(1)

    # Enable path normalization
    if resolve:
        samplesheet_file = resolve_samplesheet_paths(samplesheet_file, outdir)
        logger.info("Resolved the file paths stated in the samplesheet")

    logger.debug("Creating pipeline configuration file")
    config_file = create_config(samplesheet_file, 
                           outdir, 
                           deploy_dir,
                           args.verbosity
                           )

    if ignore_assemblies or force:
        logger.info("Assemblies in samplesheet will not replace "
                    "assembly steps in the pipeline. This might take some time!")
    else:
        logger.info("Assemblies in samplesheet will be used to skip "
                    "assembly steps in the pipeline, where applicable!")
        link_assemblies(samplesheet_file, str(SPE_CONFIGS), outdir)


    logger.debug("Creating pipeline command")
    command = create_command(threads,
                             config_file, 
                             conda_dir,
                             force
                             )

    logger.info(
        "Executing pipeline for Mixed Microbial Analysis on Sequencing data"
    )
    status = execute_snakemake(command)

    if status != 0:
        logger.error("Something went wrong while executing snakemake.")
        sys.exit(1)


def launcher() -> None:
    # Parse user input
    args = parse_mmaseq()

    # Initiate logging
    pkg_logging.adjust_log(logger, args.verbosity, args.logfile)

    mmaseq(args)

    logger.info("MMAseq successful!")