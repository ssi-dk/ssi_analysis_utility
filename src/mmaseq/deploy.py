from .utils import pkg_logging
from .utils.paths import *
import argparse
import subprocess
import sys
import ftplib

logger = pkg_logging.initiate_log("MMAdeploy")


def parse_deploy():
    parser = argparse.ArgumentParser(
        description = """
        Run pipeline on inbuild test set and setup module environments and databases for deployment of all modules and databases.
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
        "--small",
        dest = "small",
        action = "store_true",
        help = f"""
            Use a small dataset rather than running on  full deployment dataset. (Default False)
            The small dataset consists of a single isolate, executed on ALL modules.
            If enabled, Read data of a single isolate will be downloaded rather a selection of species.
            Read data will be downloaded to {READ_DIR}
        """
    )

    parser.add_argument(
        "--retries",
        dest = "retries",
        default = 3,
        help = f"""
            Amount of attempts allwoed for each file, when downloading the dataset. (Default 3) 
            Setting this to 0 will lead to failure if any short instance of disconnect occurs. 
            Contrarily setting this to a too high value could lead to long run time, 
            if any continous connection issue ARE occuring. 
            It's recommended to allow for a handfull of attempts.
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
        "--debug",
        dest = "debug",
        action = "store_true",
        help = """
            Add debug messages during execution. Mostly used for development 
            and debugging purposes
        """
    )

    return parser.parse_args()


def download_ftp_file(url, destination, max_retries):

    status = False
    url = url.strip()
    retries = 0
    while retries <= max_retries:
        retries += 1
        try:
            # Parse the URL
            ftp_url = url.replace('ftp://', '')
            host, *path = ftp_url.split('/')        

            # Define the target file path
            target_file = destination / path[-1]
            
            # Abort if the file exists 
            if target_file.exists():
                logger.debug((f"File {target_file.name} is skipped "
                    f"because it allready exists: {target_file}"))
                return None
            
            logger.debug(f"Connecting to FTP server at {host}")
            ftp = ftplib.FTP(host)

            # Anonymous login
            ftp.login()

            logger.info(f"Downloading {url} -> {target_file}")
            with open(target_file, 'wb') as local_file:
                ftp.retrbinary(f'RETR {"/".join(path)}', local_file.write)
            
            logger.debug(f"Downloaded successfully to {target_file}")

            status = True

        except Exception as e: # Look for Login exceptions, connection not found exceptions etc.
            logger.error((f"Failed to download {url} on attempt #{retries}\n"
                f"Error: {str(e)}"))
        finally:
            logger.debug("Closing ftp connection.")
            try:
                ftp.quit()
            except UnboundLocalError:
                logger.debug(("Closing FTP failed because it was never "
                    "established. This was expected behavior!"))
            except AttributeError as e:
                logger.error(f"Yeeeerp.. I have no idea yet what goes wrong, ignoring...\n{e}")


    return status


def deploy_dataset(small, max_retries):

    with open(URL_FILE, "r") as url_file:
        urls = url_file.readlines()

    # Reduce dataset size if small is selected
    size = "the full"
    if small:
        urls = urls[0:2]
        size = "a subselection of the"


    logger.debug((
        f"Attempting to download {size} deployment dataset "
        f"allowing up to {max_retries} retries on each file."))

    for url in urls:
        downloaded = download_ftp_file(url, READ_DIR, max_retries)

        # Want to introduce status messages here.


def deploy(args):

    deploy_dir = Path(args.deploy_dir)
    small = args.small
    retries = args.retries
    threads = args.threads

    logger.info("Downloading deployment dataset.")
    deploy_dataset(small, retries)

    config = f"{PKG_CONFIGS}/Test.yaml"

    samplesheet_file = f"{DATA_DIR}/samplesheet.tsv"

    # Create arguments for command
    dataset = "full"
    additional_cmds = ""
    if small:
        dataset = "small"
        samplesheet_file = f"{DATA_DIR}/samplesheet_small.tsv"
        additional_cmds = "--ignore_assemblies "


    outdir = deploy_dir / "MMAseq_Test"
    if args.debug:
        additional_cmds += "--debug "

    # Create command
    command = (
        f"mmaseq --samplesheet {samplesheet_file} "
        f"--deploy_dir {deploy_dir} "
        f"--outdir {outdir} "
        f"--threads {threads} "
        f"--config {config} "
        "--resolve "
        f"{additional_cmds}"
    )

    logger.debug(f"Executing MMAseq with following command: {command}")
    status = subprocess.Popen(command, shell = True).wait()

    if status != 0:
        logger.error((
            "Something went wrong during deployment. "
            "Rerun command with '--debug' for more details."
        ))
        sys.exit(1)

    else:
        logger.info((
            f"Deployment complete on {dataset} dataset. "
            f"Environments isntalled and databases downloaded to {deploy_dir}.\n"
            f"Results from Test dataset stored in {outdir}"
        ))


def launcher() -> None:
    # Parse user input
    args = parse_deploy()

    # Generate logger
    pkg_logging.adjust_log_level(logger, args.debug)

    logger.info("Running MMAseq Deployment")
    deploy(args)

    logger.info("Deployment successful!")