import os
import sys
import logging
from datetime import datetime

def setup_logging(log_dir: str, sample_name: str, notes: str) -> None:
    """
    Sets up logging to a file and console for a given sample.

    Args:
        log_dir (str): Path to the directory where logs should be saved.
        sample_name (str): Identifier for the sample (used in log file name).

    Returns:
        None
    """
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    
    # Ensure log directory exists
    os.makedirs(log_dir, exist_ok=True)

    # Construct log file path correctly
    log_file = os.path.join(log_dir, f"{sample_name}_{notes}.log")

    # Get root logger
    logger = logging.getLogger()

    # Remove all handlers to reset logging to a new file
    while logger.hasHandlers():
        logger.removeHandler(logger.handlers[0])

    # Set up logging with a new file for each run
    logging.basicConfig(
        filename=log_file,
        filemode="a",  # Append mode
        format="%(asctime)s - %(levelname)s - %(message)s",
        level=logging.INFO
    )

    # Add console logging
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(logging.Formatter("%(message)s"))
    logger.addHandler(console_handler)

    logging.info(f"Logging started for {log_file}")
