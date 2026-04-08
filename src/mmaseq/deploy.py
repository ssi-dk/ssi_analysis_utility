from .utils import pkg_logging
from .utils.paths import *
import argparse
import subprocess
import sys
import collections
import ftplib

logger = pkg_logging.initiate_log("MMAdeploy")


def parse_deploy():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description = (
            "MMAseq deploy\n"
            "Install environments and creates databases by "
            "executing MMAseq on an inbuilt test dataset.\n"
            "Pipeline results files are written to the deployment directory."
        ),
        epilog = (
            "This is the MMAseq Deploy module.\n"
            "For details on samplesheet creation execute 'mmacreate -h'\n"
            "For details on the main module execute 'mmaseq -h'"
        )
    )

    parser.add_argument(
        "--deploy_dir",
        dest = "deploy_dir",
        default = PKG_DIR / "Deploy",
        help = f"""
            Directory used to deploy virtual environment and databases 
            used during pipeline execution. To reinstall environments 
            and/or databases, remove the `conda/` and/or the `Databases/` 
            folders in the deployment directory. (Default: %(default)s)
        """
    )

    parser.add_argument(
        "--small",
        dest = "small",
        action = "store_true",
        help = f"""
            Use a small dataset rather than running on  full deployment dataset. (Default: %(default)s)
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
            Amount of attempts allwoed for each file, when downloading the dataset. (Default: %(default)s) 
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
            (Default: %(default)s)
        """
    )

    parser.add_argument(
        "--verbosity",
        dest = "verbosity",
        type = int,
        choices = [0, 1, 2],
        default = 0,
        help = (
            "Adjust the verbosity (Default: %(default)s); "
            "0: Minimal messages, "
            "1: Debug messages, "
            "2: Trace messages (deveolpment only)"
        )
    )

    parser.add_argument(
        "--logfile",
        dest = "logfile",
        type = str,
        default = None,
        help = (
            "If provided, will redirect log messages from STDOUT to logfile. (Default: %(default)s) "
            "Will be ignored if logfile parent folder deosn't exists."
        )
    )


    return parser.parse_args()


def extract_hosts(urls):
    logger.trace(f"extract_hosts(urls = {urls})")
    hosts = collections.defaultdict(list)
    for url in urls:
        host, *path_list = url.replace('ftp://', '').split('/')
        path = f"/{'/'.join(path_list)}"
        hosts[host].append(path)

    return hosts


def connect_ftp(host):
    logger.trace(f"Connecting to {host}")

    ftp = ftplib.FTP(host)

    # Anonymous login
    ftp.login()

    return(ftp)


def download_ftp_file(ftp, paths, destination, max_retries):

    logger.trace(("download_ftp_file(\n - "
        f"ftp: {ftp}\n - "
        f"paths: {paths}\n - "
        f"destination: {destination}\n - "
        f"max_retries: {max_retries})"))

    for path in paths:

        # Define target file and download chunk file
        target_file = destination / path.split('/')[-1]
        target_chnk = target_file.with_suffix(f"{target_file.suffix}.chunk")
                
        # Remove old chunks if allready exists
        if target_chnk.exists():
            logger.warning(
                f"Old chunk file detected: {target_chnk} "
                "this is not intended. Removing!"
            )
            target_chnk.unlink()

        # Abort if the file exists 
        if target_file.exists():
            logger.debug((
                f"File allready downloaded. Skipping {target_file.name}"
            ))
            continue

        logger.info(f"Test sample missing. Downloading {target_file.name}")

        success = False    
        try:

            logger.trace(
                f"Downloading {target_file.name} as {target_chnk}"
            )

            with open(target_chnk, 'wb') as local_file:
                ftp.retrbinary(f'RETR {path}', local_file.write)

            logger.trace(
                f"Renaming {target_chnk.name} to {target_file.name}"
            )
            target_chnk.replace(target_file)
            
            success = True

        except Exception as e:
            logger.error((
                f"Something went wrong with download.\n"
                f"{e}"
            ))
        finally:
            if target_file.exists() and not success:
                logger.warning("Download was unsuccessful, "
                    f"but target does exist: {target_file}. "
                    "Something is wrong - Deleting!")
                target_file.unlink()

        # Want to introduce status messages here.
        if success:
            logger.trace(f"{target_file.name} was successfully downloaded into {READ_DIR}")
        else:
            logger.warning(f"{target_file.name} failed to download!")

    return None


def deploy_dataset(small, max_retries):

    logger.trace(("deploy_dataset(\n - "
        f"small: {small}\n - "
        f"max_retries: {max_retries})"))

    with open(URL_FILE, "r") as url_file:
        urls = url_file.read().splitlines()

    # Reduce dataset size if small is selected
    size = "the full"
    if small:
        urls = urls[0:2]
        size = "a subselection of the"

    hosts = extract_hosts(urls)

    for host in hosts.keys():

        paths = hosts.get(host)

        logger.debug(f"Examining {host} for test dataset")
        
        retries = 0
        while retries <= max_retries:
            retries += 1
            try:        

                ftp = connect_ftp(host)

                download_ftp_file(ftp, paths, READ_DIR, max_retries)

                retries = max_retries + 1

            except Exception as e:
                logger.error((f"Failed to download {paths} on attempt #{retries}\n"
                    f"Error: {str(e)}"))
            finally:
                logger.trace(f"Attempting to close the FTP connection to {host}.")
                try:
                    ftp.quit()
                    logger.trace("FTP connection soft close successful!")
                except UnboundLocalError:
                    logger.trace(("Closing FTP failed because it was never "
                        "established in the first place. "
                        "This was expected behavior!"))
                except AttributeError as e:
                    logger.warning((
                        "FTP connection can't be terminated softly. "
                        "Attempting aggressive termination!"))
                    try:
                        ftp.close()
                        logger.trace("FTP connection aggressive close successful!")
                    except Exception as e:
                        logger.error(f"FTP disconnect failed. ")


def deploy(args):

    deploy_dir = Path(args.deploy_dir)
    small = args.small
    retries = args.retries
    threads = args.threads
    verbosity = args.verbosity

    logger.info(f"Inspecting the deployment dataset")
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
    additional_cmds += f"--verbosity {verbosity}"

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
    logger.debug(f"Created command for MMAseq:\n{command}")

    logger.info(f"Executing MMAseq")
    status = subprocess.Popen(command, shell = True).wait()

    if status != 0:
        logger.error((
            "Something went wrong during deployment. "
            "Rerun command with '--verbosity 1' for more details."
        ))
        sys.exit(1)

    else:
        logger.info((
            f"Deployment complete on {dataset} dataset. "
            f"Environments installed and databases downloaded to {deploy_dir}.\n"
            f"Results from Test dataset stored in {outdir}"
        ))


def launcher() -> None:

    # Parse user input
    args = parse_deploy()

    # Adjust logger
    pkg_logging.adjust_log_level(logger, args.verbosity)

    deploy(args)

    logger.info("Deployment successful!")
