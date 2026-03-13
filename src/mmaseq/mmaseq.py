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

# Initiate logging
logger = pkg_logging.initiate_log("MMAseq")

def resolve_path(path, cwd, samplesheet_file):

    # Determine possible file paths
    path = Path(path)
    path_absolute = absolute.is_absolute()
    path_from_cwd = cwd / path
    samplesheet_dir = Path(samplesheet_file).resolve().parent
    path_from_samplesheet_dir = samplesheet_dir / path

    if path_absolute:
        absolute = path
    elif not path_absolute & path_from_cwd.exists():
        absolute = path_from_cwd.resolve()
    elif not path_absolute & path_from_samplesheet_dir.exists(): # I suggest removing this, as default behavoir is files relative to cwd NOT input file
        absolute = path_from_samplesheet_dir.resolve()
    else:
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
                  config_file, 
                  force):

    if not os.path.isfile(config_file):
        config_dir = os.path.dirname(config_file)
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
    
    outdir = os.path.abspath(outdir).rstrip("/")

    config = {
        "samplesheet": str(samplesheet_file),
        "deploy_dir": str(deploy_dir),
        "outdir": str(outdir)
    }

    with open(config_file, "w") as config_yaml:
        yaml.safe_dump(config, config_yaml)

    return config_file


def link_assemblies(samplesheet_file, 
                    config_dir, 
                    outdir):

    logger.debug("Reading samplesheet")
    samplesheet = pd.read_csv(samplesheet_file, 
                              sep = "\t").set_index("sample_name")


    logger.debug("Importing sample configs")
    sample_configs = determine_sample_configs(samplesheet = samplesheet, 
                                              config_dir = config_dir)

    logger.debug("Initiating symlink generation")

    for sample, configs in sample_configs.items():

        assembly_sheet = samplesheet.at[sample, "assembly"]

        # Case 1: NA or missing
        if pd.isna(assembly_sheet):
            logger.debug(f"No assembly provided for {sample} — skipping!")
            continue

        assembly_source = Path(str(assembly_sheet))

        # Case 2: Path is invalid or does not exist
        if not assembly_source.exists():
            logger.warning(f"""Assembly path does not exist for {sample} 
                at location {assembly_sheet}\nSkipping!""")
            continue

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

        for assembler in assemblers:

            assembly_dir = outdir / sample / assembler

            if not assembly_dir.exists():
                assembly_dir.mkdir(parents = True)

            destination = assembly_dir / f"{sample}.fasta"

            if destination.is_symlink() and not destination.exists(follow_symlinks = True):
                destination.unlink()

            if not destination.exists(follow_symlinks = False):
                logger.info(f"""Creating symlink: 
                    {assembly_source} -> {destination}""")
                destination.symlink_to(assembly_source)
            else:
                logger.warning(f"""File exists and is not a symlink {destination}.\n
                    Skipping!""")


def create_command(threads, 
                   config, 
                   conda_dir, 
                   arguments = None, 
                   rules = None):

    additionals = " ".join(arguments) if arguments else ""
    target_rules = " ".join(rules) if rules else ""

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

    return command


def execute_snakemake(command):
    status = subprocess.Popen(command, 
                              shell = True).wait()
    return status


def mmaseq(args):

    samplesheet_file = Path(args.samplesheet_file)
    input_dir = Path(args.input_dir)
    deploy_dir = Path(args.deploy_dir)
    outdir = Path(args.outdir)
    threads = args.threads
    config = args.config
    resolve = args.reslove
    ignore_assemblies = args.ignore_assemblies

    # Outdated
    # force = False
    # if config is None:
    #     config = PKG_CONFIGS / "config.yaml"
    #     force = True

    conda_dir = (deploy_dir / "conda").resolve()
    arguments = []
    rules = []

    # if test:
    #     logger.info("Test run initiated. Will ignore irrelevant user arguments!")
    #     config = f"{PKG_CONFIGS}/Test.yaml"
    #     samplesheet_file = f"{DATA_DIR}/samplesheet.tsv"
    #     outdir = outdir / "Test"
    #     rules.append("all")

    #     # Fix: normalize paths in test
    #     if resolve:
    #         resolve_samplesheet_paths(samplesheet_file)

    #     config = create_config(samplesheet_file, 
    #                            outdir, 
    #                            deploy_dir, 
    #                            config, 
    #                            force = True)

    #     link_assemblies(samplesheet_file, 
    #                     SPE_CONFIGS, 
    #                     outdir)

    #     command = create_command(threads, 
    #                              config, 
    #                              conda_dir, 
    #                              arguments, 
    #                              rules)


        # if samplesheet_file is None:
        #     logger.error("Samplesheet file user argument missing. Aborting!")
        #     sys.exit(1)

    if not os.path.isfile(samplesheet_file):
        logger.error("""Samplesheet not found.\n - Execute the create module 
            first.\nAborting!""")
        sys.exit(1)

    # Fix: normalize user paths
    if resolve:
        resolve_samplesheet_paths(samplesheet_file)

    config = create_config(samplesheet_file, 
                           outdir, 
                           deploy_dir, 
                           config, 
                           force)

    if not ignore_assemblies:
        logger.info("")
        link_assemblies(samplesheet_file, SPE_CONFIGS, outdir)

    command = create_command(threads, 
                             config, 
                             conda_dir, 
                             arguments, 
                             rules)

    logger.info("Executing pipeline: Mixed Microbial Analysis on Sequencing data")
    status = execute_snakemake(command)

    if status != 0:
        logger.error("Something went wrong while executing snakemake")
    else:
        logger.info("Pipeline successful")


def launcher() -> None:
    args = parse_mmaseq()

    adjust_log_level(logger, args.debug)

    mmaseq(args)
