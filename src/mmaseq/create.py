from .utils import pkg_logging, parse_create, adjust_log_level 
from .utils import *
from pathlib import Path
import re
import pandas as pd

# Initiate logging
logger = pkg_logging.initiate_log("MMAseq Create")

def create_samplesheet(args):

    indir = Path(args.indir)
    outdir = Path(args.outdir)

    if not indir.is_dir():
        logger.error(f"Input directory doesn't exist. Aborting!\n - {indir}")
        sys.exit(1)

    root = indir.expanduser().resolve()

    records = {}

    logger.info(f"Scanning for sample files in {indir}")

    for path in root.rglob("*"):

        file = Path(path.name)
        
        # Only examine files
        if path.is_file():

            # Examining assembly file
            if path.suffix in [".fasta", ".fa"]: # Was fname
                logger.debug(f"Found assembly: {path}")
                sample = str(file.with_suffix(""))

                # Ensure sample is created with reasonable defaults
                records.setdefault(sample, {
                    "read1": "NA",
                    "read2": "NA",
                    "assembly": "NA",
                    "config": "default.yaml"
                })

                # Record assembly file location
                records[sample]["assembly"] = str(path)

            # Examining read file
            if path.suffixes in [[".fastq", ".gz"], [".fq", ".gz"]]:
                logger.debug(f"Found read: {path}")
                sample = re.sub(r'(_R?[12])\D*\.fastq\.gz$', '', str(file))

                # Ensure sample is created with reasonable defaults
                records.setdefault(sample, {
                    "read1": "NA",
                    "read2": "NA",
                    "assembly": "NA",
                    "config": "default.yaml"
                })

                # Search for read mate annotations and record file location
                if re.search(r"(_R?1\D?|_1\D?)", str(file)):
                    records[sample]["read1"] = str(path)
                elif re.search(r"(_R?2\D?|_2\D?)", str(file)):
                    records[sample]["read2"] = str(path)
                else:
                    logger.warning(f"""Read mate {path} was not recognized. 
                        Inspect samplesheet manually! Ignoring file for {sample}...""")

    samplesheet = (
        pd.DataFrame.from_dict(records, orient = "index")
        .reset_index()
        .rename(columns = {"index": "sample_name"})
    )

    samplesheet_file = outdir / "samplesheet.tsv"
    try:
        if not outdir.exists():
            outdir.mkdir(parents = True)

        samplesheet.to_csv(samplesheet_file, 
            sep = "\t",
            index = False)
    except PermissionError as e:
        logger.error(f"""You don't have permission to write to: 
            {os.path.dirname(samplesheet_file)}\n - Aborting!\n{e}""")
        sys.exit(1)

    return samplesheet_file


def launcher():

    args = parse_create()

    adjust_log_level(logger, args.debug)

    samplesheet_file = create_samplesheet(args)

    logger.info(f"Samplesheet successfully created: {samplesheet_file}")
