from .utils import system_paths, pkg_logging, parse_create, adjust_log_level 
from pathlib import Path
import re
import pandas as pd

# Initiate logging
logger = pkg_logging.initiate_log("MMAseq Create")

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

                logger.debug(f"Found {sample} with assembly file: {path}")
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

                logger.debug(f"Found {sample} with read file: {path}")
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
            {os.path.dirname(samplesheet_file)}\n - Aborting!\n{e}""")
        sys.exit(1)

    return samplesheet_file


def launcher():
    # Read user arguments
    args = parse_create()

    # Generate logger
    adjust_log_level(logger, args.debug)

    logger.info("Initiating samplesheet creation")
    samplesheet_file = create_samplesheet(args)

    logger.info(f"Samplesheet successfully created: {samplesheet_file}")
