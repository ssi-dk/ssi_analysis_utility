import logging


def initiate_log(name):
    """
    Initiates a logger with the given name.

    Args:
        name (str): Name of the logger.

    Returns:
        logging.Logger: The logger instance.
    """
    logger = logging.getLogger(name)

    return logger


def adjust_log(logger, verbosity, logfile=None):
    """
    Adjusts the logger with verbosity and optional logfile.

    Args:
        logger (logging.Logger): The logger to adjust.
        verbosity (int): Verbosity level (0, 1, 2).
        logfile (str, optional): Path to logfile.
    """
    def trace(self, msg, *args, **kwargs):
        if self.isEnabledFor(TRACE_LEVEL):
            self._log(TRACE_LEVEL, msg, args, **kwargs)

    # Adding TRACE level
    TRACE_LEVEL = 5
    logging.addLevelName(TRACE_LEVEL, "TRACE")
    logging.TRACE = TRACE_LEVEL
    logging.Logger.trace = trace

    # Define text formatting options
    logformat_simple = "[%(name)s] %(levelname)s: %(message)s"
    logformat_advanced = "[%(name)s] %(levelname)s %(asctime)s: %(message)s"

    # Determine verbosity settings
    level = logging.INFO
    logformat = logformat_simple
    if int(verbosity) == 1:
        level = logging.DEBUG
        logformat = logformat_advanced
    elif int(verbosity) == 2:
        level = logging.TRACE
        logformat = logformat_advanced
    elif int(verbosity) != 0:
        print(
            "Warning: "
            f"Couldn't interpret verbosity level, it was set to {verbosity}.\n"
            "Ignoring user input, verbosity set to 0!"
        )

    # Generate log formatter
    _formatter = logging.Formatter(logformat, "%H:%M:%S")

    _handler = logging.StreamHandler()
    _handler.setFormatter(_formatter)

    # Setup logger object
    logger.addHandler(_handler)
    logger.propagate = False
    logger.setLevel(level)

    # Define log file
    if logfile is not None:
        print(f"Writing messages to logfile: {logfile}")

        # Define where to write logfile
        logging.basicConfig(
            filename=logfile,
            filemode="w",
            format="%(asctime)s - %(levelname)s: %(message)s",
            level=level
        )
