import os
import sys
import logging
from datetime import datetime

def setup_logging(log_file: str) -> None:
    """
    Sets up logging to a file and console for a given sample.

    Args:
        log_file (str): Path to the log file, missing directoreis will be created.

    Returns:
        None
    """

    log_dir = os.path.dirname(log_file)
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    
    # Ensure log directory exists
    if not os.path.isdir(log_dir):
        print("Creating directory for Log file", flush = True)
        os.makedirs(log_dir, exist_ok = True)

    # Get root logger
    logger = logging.getLogger()

    # Remove all handlers to reset logging to a new file
    while logger.hasHandlers():
        logger.removeHandler(logger.handlers[0])

    # Set up logging with a new file for each run
    logging.basicConfig(
        filename=log_file,
        filemode="w",
        format="%(asctime)s - %(levelname)s - %(message)s",
        level=logging.DEBUG
    )

    # Add console logging
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(logging.Formatter("%(message)s"))
    logger.addHandler(console_handler)

    logging.info(f"===============================================")
    logging.debug(f"Logging to {log_file}")
