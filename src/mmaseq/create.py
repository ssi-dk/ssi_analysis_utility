from .utils import pkg_logging
from .utils.paths import *
import argparse
from pathlib import Path
import re
import pandas as pd
import sys

# Initiate logging
logger = pkg_logging.initiate_log("MMAcreate")


def parse_create():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description = (
            "MMAseq create\n"
            "Create a samplesheet by scanning the input directory for sample files."
        ),
        epilog = (
            "This is the MMAseq Create module.\n"
            "For details on the deployment module execute 'mmadeploy -h'"
            "For details on the main module execute 'mmaseq -h'"
        )
    )

    parser.add_argument(
        "--indir",
        dest = "indir",
        default = None,
        help = """
            Input directory MUST be specified if the samplesheet does 
            not yet exist. (Default None)
            Input directory will be screened for `.fasta` and `fastq.gz` 
            files, sample_names will be infered from the detected files, 
            and used to populate a samplesheet. After samplesheet creation, 
            the pipeline will be executed in dry-run mode (simulated run)
        """
    )

    parser.add_argument(
    "--outdir",
    dest = "outdir",
    default = CWD / "MMAseq_Results",
    help = """
        Directory used for storing the samplesheet. (Default ./MMAseq_Results)
        Samplesheet will be stored in ourdir as 'samplesheet.tsv'
    """
    )

    parser.add_argument(
        "--verbosity",
        dest = "verbosity",
        default = 0,
        help = """
            Adjust the verbosity level of running with integers between 0 and 2.
            0: Show standard messages
            1: Provide debug messages (Usefull for inspecting errors)
            2: Provide detailed trace messages (Usefull for development)

            d debug messages during execution. (Default False) 
            Mostly used for development and debugging purposes.
        """
    )

    return parser.parse_args()


def create_samplesheet(args):
    # Define user arguments
    indir = Path(args.indir)
    outdir = Path(args.outdir)

    # Ensure that input directory exists
    if not indir.is_dir():
        logger.error(f"Input directory doesn't exist. Aborting!\n - {indir}")
        sys.exit(1)


    logger.debug(f"Scanning for sample files in {indir}")
    root = indir.expanduser().resolve()
    records = {}
    # Iterate recursively over everything in input directory
    for path in root.rglob("*"):
        file = Path(path.name)
        
        # Only examine files
        if path.is_file():

            # Investegate assembly files
            if path.suffix in [".fasta", ".fa"]:
                sample = str(file.with_suffix(""))

                logger.trace(f"Found {sample} with assembly file: {path}")
                records.setdefault(sample, {
                    "read1": "NA",
                    "read2": "NA",
                    "assembly": "NA",
                    "config": "default.yaml"
                })

                # Record assembly file location
                records[sample]["assembly"] = str(path)

            # Investegate read files
            if path.suffixes in [[".fastq", ".gz"], [".fq", ".gz"]]:
                sample = re.sub(r'(_R?[12])\D*\.fastq\.gz$', '', str(file))

                logger.trace(f"Found {sample} with read file: {path}")
                records.setdefault(sample, {
                    "read1": "NA",
                    "read2": "NA",
                    "assembly": "NA",
                    "config": "default.yaml"
                })

                # Record read mates based on specified suffices
                if re.search(r"(_R?1\D?|_1\D?)", str(file)):
                    records[sample]["read1"] = str(path)
                elif re.search(r"(_R?2\D?|_2\D?)", str(file)):
                    records[sample]["read2"] = str(path)
                else:
                    # Provide warnings for unknown naming convensions
                    logger.warning(f"""Read mate {path} was not recognized. 
                        Inspect samplesheet manually! Ignoring file for {sample}...""")

    # Collect records into a data frame
    samplesheet = (
        pd.DataFrame.from_dict(records, orient = "index")
        .reset_index()
        .rename(columns = {"index": "sample_name"})
    )
    logger.info("Creating samplesheet")

    samplesheet_file = outdir / "samplesheet.tsv"
    try:
        if not outdir.exists():
            logger.debug(f"Creating output directory {outdir}")
            outdir.mkdir(parents = True)

        logger.debug(f"Writing sample sheet to {samplesheet_file}")
        samplesheet.to_csv(samplesheet_file, 
            sep = "\t",
            index = False)
    # Handle permission errors
    except PermissionError as e:
        logger.error(f"""You don't have permission to write to: 
            {samplesheet_file.parrent}\n - Aborting!\n{e}""")
        sys.exit(1)

    return samplesheet_file


def launcher():
    # Read user arguments
    args = parse_create()

    # Generate logger
    pkg_logging.adjust_log_level(logger, args.verbosity)

    logger.info("Initiating samplesheet creation")
    samplesheet_file = create_samplesheet(args)

    logger.info(f"Samplesheet creation successful!")
