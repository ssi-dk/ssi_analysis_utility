#!/usr/bin/env python

import argparse
import logging
import os
import re
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from importlib import resources

import pandas as pd
import yaml

from .utils.helper_functions import determine_sample_configs

# Determining system paths
PKG_DIR = resources.files("mmaseq")
WORKFLOW_DIR = PKG_DIR / "workflow"
SNAKEFILE =  WORKFLOW_DIR / "Snakefile"
CONFIG_DIR = PKG_DIR / "config"
SPE_CONFIG_DIR = CONFIG_DIR / "species_configs"
CWD = Path.cwd()

# Determining test data path
DATA_DIR = PKG_DIR / "data"    

# Logging configuration
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
        description="""
        Mixed Microbial Analysis on Sequencing data (MMAseq)
        """
    )

    parser.add_argument(
        "--samplesheet",
        dest="samplesheet_file",
        default=None,
        help="""
            Path to samplesheet TSV used by the pipeline. 
            If the samplesheet doesn't exist, `input_dir` must be 
            specified to create one.
        """
    )

    parser.add_argument(
        "--input_dir",
        dest="input_dir",
        default=None,
        help="""
            Input directory MUST be specified if the samplesheet does 
            not yet exist. (Default None)
            Input directory will be screened for `.fasta` and `fastq.gz` 
            files, sample_names will be infered from the detected files, 
            and used to populate a samplesheet. After samplesheet creation, 
            the pipeline will be executed in dry-run mode (simulated run)
        """
    )

    parser.add_argument(
        "--deploy_dir",
        dest="deploy_dir",
        default=PKG_DIR / "Deploy",
        help=f"""
            Directory used to deploy virtual environment and databases 
            used during pipeline execution. To reinstall environments 
            and/or databases, remove the `conda/` and/or the `Databases/` 
            folders in the deployment directory. (Default {PKG_DIR}/Deploy)
        """
    )

    parser.add_argument(
        "--outdir",
        dest="outdir",
        default=CWD / "Results",
        help="""
            Directory used for storing analysis results. (Default ./Results)
        """
    )

    parser.add_argument(
        "--threads",
        dest="threads",
        default=4,
        help="""
            Amount of threads (cores) to dedicate for executing the pipeline. 
            (Default 4)
        """
    )

    parser.add_argument(
        "--config",
        dest="config",
        default=CONFIG_DIR / "config.yaml",
        help=f"""
            Configuration file location. (Default {CONFIG_DIR}/config.yaml)
        """
    )

    parser.add_argument(
        "--test",
        dest="test",
        action="store_true",
        help="""
            Perform test run of all modules. Can be used to install all 
            conda environments and databases. The config file will be 
            generated as {CONFIG_DIR}/Test.yaml and output will be stored 
            in a Test/ folder of the specified output directory. 
            (Default ./Results/Test)
        """
    )

    parser.add_argument(
        "--debug",
        dest="debug",
        action="store_true",
        help="""
            Add debug messages during execution. Mostly used for development 
            and debugging purposes
        """
    )

    return parser.parse_args()



# Path resolution from samplesheet
# ---------------------------------------------------------------------
def resolve_sample_path(path_from_sheet, 
                        samplesheet_file):

    p = Path(path_from_sheet)

    if p.is_absolute():
        return p

    samplesheet_dir = Path(samplesheet_file).resolve().parent
    return (samplesheet_dir / p).resolve()



# NEW: Normalize samplesheet paths (fix)
# ---------------------------------------------------------------------
def normalize_samplesheet_paths(samplesheet_file):
    """
    Converts all read1/read2/assembly paths in the samplesheet
    into absolute paths using resolve_sample_path().
    This guarantees that Snakemake always receives valid paths.
    """
    df = pd.read_csv(samplesheet_file, sep="\t").copy()

    def fix(p):
        if isinstance(p, str) and p.lower() not in ["na", ""]:
            return str(resolve_sample_path(p, 
                                           samplesheet_file))
        return p

    df["read1"] = df["read1"].apply(fix)
    df["read2"] = df["read2"].apply(fix)
    df["assembly"] = df["assembly"].apply(fix)

    normalized = Path(samplesheet_file).with_suffix(".normalized.tsv")
    df.to_csv(normalized, 
              sep="\t", 
              index=False)

    return str(normalized)



# Samplesheet creation
# ---------------------------------------------------------------------
def create_samplesheet(samplesheet_file, 
                       input_dir):

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
            sample = fname.replace(".fasta", "")
            records.setdefault(sample, {
                "read1": "NA",
                "read2": "NA",
                "assembly": "NA",
                "config": "default.yaml"
            })
            records[sample]["assembly"] = str(path)

        elif fname.endswith(".fastq.gz"):
            sample = re.sub(r'(_R?[12].*)\.fastq\.gz$', '', fname)
            records.setdefault(sample, {
                "read1": "NA",
                "read2": "NA",
                "assembly": "NA",
                "config": "default.yaml"
            })

            if re.search(r'(_R?1|_1)', fname):
                records[sample]["read1"] = str(path)
            elif re.search(r'(_R?2|_2)', fname):
                records[sample]["read2"] = str(path)

    samplesheet = (
        pd.DataFrame.from_dict(records, orient="index")
        .reset_index()
        .rename(columns={"index": "sample_name"})
    )

    samplesheet.to_csv(samplesheet_file, 
                       sep="\t", 
                       index=False)
    return None



# Configuration file creation
# ---------------------------------------------------------------------
def create_config(samplesheet_file, 
                  outdir, 
                  deploy_dir, 
                  config_file, 
                  root, 
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
        "outdir": str(outdir),
        "root": str(root)
    }

    with open(config_file, "w") as config_yaml:
        yaml.safe_dump(config, config_yaml)

    return config_file



# Create assembly symlinks
# ---------------------------------------------------------------------
def link_assemblies(samplesheet_file, 
                    config_dir, 
                    outdir):

    logger.debug("Reading samplesheet")
    samplesheet = pd.read_csv(samplesheet_file, 
                              sep="\t").set_index("sample_name")

    logger.debug("Importing sample configs")
    sample_configs = determine_sample_configs(samplesheet=samplesheet, 
                                              config_dir=config_dir)

    outdir = Path(outdir)

    logger.debug("Initiating symlink generation")

    for sample, configs in sample_configs.items():

        assembly_sheet = samplesheet.at[sample, "assembly"]

        # Case 1: NA or missing
        if pd.isna(assembly_sheet):
            logger.warning("No assembly provided for sample %s — skipping.", sample)
            continue

        assembly_source = Path(str(assembly_sheet))

        # Case 2: Path is invalid or does not exist
        if not assembly_source.exists():
            logger.warning("Assembly path does not exist for sample %s: %s", sample, assembly_sheet)
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
            assembly_dir.mkdir(parents=True, 
                               exist_ok=True)

            destination = assembly_dir / f"{sample}.fasta"

            if destination.is_symlink():
                destination.unlink()

            if not destination.exists():
                logger.info("Creating symlink: %s -> %s", 
                            assembly_source, 
                            destination)
                destination.symlink_to(assembly_source)
            else:
                logger.warning("File exists and is not a symlink: %s", destination)



# Command creation
# ---------------------------------------------------------------------
def create_command(threads, 
                   config, 
                   conda_dir, 
                   arguments=None, 
                   rules=None):

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



# Execute snakemake
# ---------------------------------------------------------------------
def execute_snakemake(command):
    status = subprocess.Popen(command, 
                              shell=True).wait()
    return status



# Main
# ---------------------------------------------------------------------
def mmaseq(args):

    samplesheet_file = args.samplesheet_file
    input_dir = args.input_dir
    deploy_dir = args.deploy_dir
    outdir = args.outdir
    threads = args.threads
    config = args.config
    test = args.test


    force = False
    if config is None:
        config = CONFIG_DIR / "config.yaml"
        force = True

    conda_dir = (Path(deploy_dir) / "conda").resolve()
    arguments = []
    rules = []

    if test:
        logger.info("Test run initiated. Will ignore irrelevant user arguments!")
        config = f"{CONFIG_DIR}/Test.yaml"
        samplesheet_file = f"{DATA_DIR}/samplesheet.tsv"
        outdir = outdir / "Test"
        rules.append("all")

        # Fix: normalize paths in test
        samplesheet_file = normalize_samplesheet_paths(samplesheet_file)

        config = create_config(samplesheet_file, 
                               outdir, 
                               deploy_dir, 
                               config, 
                               PKG_DIR, 
                               force=True)

        link_assemblies(samplesheet_file, 
                        SPE_CONFIG_DIR, 
                        outdir)

        command = create_command(threads, 
                                 config, 
                                 conda_dir, 
                                 arguments, 
                                 rules)

    else:

        if samplesheet_file is None:
            logger.error("Samplesheet file user argument missing. Aborting!")
            sys.exit(1)

        if not os.path.isfile(samplesheet_file):
            logger.info("Samplesheet not found, attempting to create at: %s", samplesheet_file)
            create_samplesheet(samplesheet_file, 
                               input_dir)
            arguments.append("--dry-run")
        else:
            logger.warning("Samplesheet exists and input_dir is also specified. Ignoring input_dir.")

        # Fix: normalize user paths
        samplesheet_file = normalize_samplesheet_paths(samplesheet_file)

        config = create_config(samplesheet_file, 
                               outdir, 
                               deploy_dir, 
                               config, 
                               PKG_DIR, 
                               force)

        link_assemblies(samplesheet_file, 
                        SPE_CONFIG_DIR, 
                        outdir)

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



# Launcher
# ---------------------------------------------------------------------
def launcher() -> None:
    args = parse_arguments()

    if args.debug:
        logger.setLevel(logging.DEBUG)

    mmaseq(args)
