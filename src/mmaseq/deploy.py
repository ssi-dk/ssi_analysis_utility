from .utils import logging_setup
from .utils.PATH import *
import argparse
import subprocess
import sys
import collections
import ftplib
import shutil

logger = logging_setup.initiate_log("MMAdeploy")


def parse_deploy():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=(
            "MMAseq deploy\n"
            "Install environments and creates databases by "
            "executing MMAseq on an inbuilt test dataset.\n"
            "Pipeline results files are written to the deployment directory."
        ),
        epilog=(
            "This is the MMAseq Deploy module.\n"
            "For details on samplesheet creation execute 'mmacreate -h'\n"
            "For details on the main module execute 'mmaseq -h'"
        )
    )

    parser.add_argument(
        "--deploy_dir",
        dest="deploy_dir",
        default=PKG_DIR / "Deploy",
        help=(
            "Directory used to deploy virtual environment and databases "
            "used during pipeline execution. To reinstall environments "
            "and/or databases, remove the `conda/` and/or the `Databases/` "
            "folders in the deployment directory. (Default: %(default)s)"
        )
            
    )

    parser.add_argument(
        "--update",
        dest="update",
        action="store_true",
        help=(
            "Will force running all rules to ensure issuing database updates. "
            "(Default: %(default)s) The small dataset consists of a single "
            "isolate, executed on ALL modules, thus all results should be "
            f"considered wrong. Read data will be downloaded to {READ_DIR}"
        )
    )

    parser.add_argument(
        "--test",
        dest="test",
        action="store_true",
        help=(
            "Will run rules to complete a quick pipeline test. "
            "(Default: %(default)s) The test dataset consist of exactly "
            "400001 paired end reads created synthetically from AI. "
            "Certain modules will fail on these reads and are "
            " excluded from the test. "
            "Excluded; resfinder, pointfinder, kleborate, shovill"
        )
    )

    parser.add_argument(
        "--retries",
        dest="retries",
        default=3,
        type=int,
        help=(
            "Amount of attempts allowed for each file, when downloading the "
            "dataset. (Default: %(default)s) Setting this to 0 will lead to "
            "failure if any short instance of disconnect occurs. Contrarily, "
            "setting this to a too high value could lead to long run time if "
            "any continuous connection issue occurs. It's recommended to allow "
            "for a handful of attempts."
        )
    )

    parser.add_argument(
        "--threads",
        dest="threads",
        default=4,
        help=(
            "Amount of threads (cores) to dedicate for executing the pipeline. "
            "(Default: %(default)s)"
        )
    )

    parser.add_argument(
        "--verbosity",
        dest="verbosity",
        type=int,
        choices=[0, 1, 2],
        default=0,
        help = (
            "Adjust the verbosity (Default: %(default)s); "
            "0: Minimal messages, "
            "1: Debug messages, "
            "2: Trace messages (development only)"
        )
    )

    parser.add_argument(
        "--logfile",
        dest="logfile",
        type=str,
        default=None,
        help=(
            "If provided, will redirect log messages from STDOUT to logfile. "
            "(Default: %(default)s) Will be ignored if logfile parent folder "
            "doesn't exist."
        )
    )


    return parser.parse_args()


def deploy_spe_configs(deploy_dir):
    logger.trace(("deploy_spe_configs(\n "
        f"deploy_dir = {deploy_dir}"
        ")"))

    spe_configs_dir = deploy_dir / "spe_configs"
    
    logger.trace("Checking whether config dir allready exists")
    if not spe_configs_dir.exists():
        logger.trace(f"Couldn't locate {spe_configs_dir}, clonig from {SPE_CONFIGS}")
        shutil.copytree(SPE_CONFIGS, spe_configs_dir)
        #SPE_CONFIGS.copy(spe_configs_dir) # Use from python 3.14+ and remove import shutils
        logger.info("Copied species configs directory from installation folder into deployment dir.")
    else:
        logger.info("Species configuration directory allready exists. Skipping!")

    return None


def extract_hosts(urls):
    logger.trace(f"extract_hosts(urls = {urls})")
    hosts = collections.defaultdict(list)
    for url in urls:
        host, *path_list = url.replace('ftp://', '').split('/')
        path = f"/{'/'.join(path_list)}"
        hosts[host].append(path)

    return hosts


def connect_ftp(host, timeout = 15):
    logger.trace(f"Connecting to {host}")

    ftp = ftplib.FTP(host, timeout = timeout)

    # Anonymous login
    ftp.login()

    return ftp


def disconnect_ftp(ftp):
    logger.trace("Attempting to close the FTP connection.")
    try:
        ftp.quit()
        logger.trace("FTP connection soft close successful!")
    except UnboundLocalError:
        logger.trace(
            "Closing FTP failed because it was never established in the "
            "first place. This was expected behavior!"
        )
    except AttributeError:
        logger.warning(
            "FTP connection can't be terminated softly. "
            "Attempting aggressive termination!"
        )
        try:
            ftp.close()
            logger.trace("FTP connection aggressive close successful!")
        except Exception:
            logger.error("FTP disconnect failed.")


def download_ftp_file(ftp, paths, destination, max_retries):

    logger.trace(
        f"download_ftp_file(\n - ftp: {ftp}\n - paths: {paths}\n - "
        f"destination: {destination}\n - max_retries: {max_retries})"
    )

    for path in paths:

        # Define target file and download chunk file
        target_file = destination / path.split('/')[-1]
        target_chnk = target_file.with_suffix(f"{target_file.suffix}.chunk")
                
        # Remove old chunks if already exists
        if target_chnk.exists():
            logger.warning(
                f"Old chunk file detected: {target_chnk} "
                "this is not intended. Removing!"
            )
            target_chnk.unlink()

        # Abort if the file exists
        if target_file.exists():
            logger.debug(f"File already downloaded. Skipping {target_file.name}")
            continue

        logger.info(f"Test sample missing. Downloading {target_file.name}")

        success = False    
        retries = 0
        while retries <= max_retries:
            retries += 1

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

                retries = max_retries + 1

            except Exception as e:
                logger.error(
                    f"Failed to download {path} on attempt #{retries}\n{e}"
                )
            finally:
                if target_file.exists() and not success:
                    logger.warning(
                        f"Download was unsuccessful, but target does exist: "
                        f"{target_file}. Something is wrong - Deleting!"
                    )
                    target_file.unlink()

        # Want to introduce status messages here.
        if success:
            logger.trace(
                f"{target_file.name} was successfully downloaded into {READ_DIR}"
            )
        else:
            logger.warning(f"{target_file.name} failed to download!")


def deploy_dataset(update, max_retries):

    logger.trace(
        f"deploy_dataset(\n - update: {update}\n - "
        f"max_retries: {max_retries})"
    )



    with open(URL_FILE, "r") as url_file:
        urls = url_file.read().splitlines()

    # Reduce dataset size if small is selected
    size = "the full"
    if update:
        urls = urls[0:2]
        size = "a subselection of the"

    hosts = extract_hosts(urls)

    for host in hosts.keys():

        paths = hosts.get(host)

        logger.debug(f"Examining {host} for test dataset")

        try:        
            ftp = connect_ftp(host)
        except TimeoutError as e:
            logger.error((
                f"ftp connection to {host} could not be established. Are you firewalled?\n"
                "Check whether ftp ports are openned (default is often 20, 21 or 990). "
                "Skipping host!"
            ))
            continue
        except Exception as e:
            logger.error(
                f"What? Something bad is going on... Skipping!!!\n{e}"
            )
            continue

        download_ftp_file(ftp, paths, READ_DIR, max_retries)
        
        disconnect_ftp(ftp)

    return None


def deploy(args):

    deploy_dir = Path(args.deploy_dir)
    update = args.update
    test = args.test
    retries = args.retries
    threads = args.threads
    verbosity = args.verbosity

    logger.info("Inspecting species configuration directory")
    deploy_spe_configs(deploy_dir)

    if not test:
        logger.info(f"Inspecting the deployment dataset")
        deploy_dataset(update, retries)

    samplesheet_file = f"{DATA_DIR}/samplesheet.tsv"

    # Create arguments for command
    additional_cmds = "--resolve "
    if update:
        dataset = "small"
        samplesheet_file = f"{DATA_DIR}/samplesheet_small.tsv"
        additional_cmds += "--ignore_assemblies --force --clean "
    elif test:
        dataset = "test"
        samplesheet_file = f"{DATA_DIR}/samplesheet_test.tsv"
        additional_cmds += "--clean "
    else:
        dataset = "full"



    outdir = deploy_dir / "MMAseq_Test"
    additional_cmds += f"--clean --verbosity {verbosity}"

    # Create command
    command = (
        f"mmaseq --samplesheet {samplesheet_file} "
        f"--deploy_dir {deploy_dir} "
        f"--outdir {outdir} "
        f"--threads {threads} "
        "--resolve "
        f"{additional_cmds}"
    )
    logger.debug(f"Created command for MMAseq:\n{command}")

    logger.info("Executing MMAseq")
    status = subprocess.Popen(command, shell=True).wait()

    if status != 0:
        verb = ""
        if verbosity < 1:
            verb = "Rerun command with '--verbosity 1' for more details.\n"
        logger.error(
            "Something went wrong during deployment. ",
            verb,
            "Either way you are welcome to post an issue on our Github."
            
        )
        sys.exit(1)

    else:
        logger.info(
            f"Deployment complete on {dataset} dataset. "
            f"Environments installed and databases downloaded to {deploy_dir}.\n"
            f"Results from Test dataset stored in {outdir}"
        )


def launcher() -> None:

    # Parse user input
    args = parse_deploy()

    # Adjust logger
    logging_setup.adjust_log(logger, args.verbosity)

    deploy(args)

    logger.info("Deployment successful!")

