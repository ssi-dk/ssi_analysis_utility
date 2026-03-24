import logging

def initiate_log(name):
    def trace(self, msg, *args, **kwargs):
        if self.isEnabledFor(TRACE_LEVEL):
            self._log(TRACE_LEVEL, msg, args, **kwargs)


    # Adding TRACE level
    TRACE_LEVEL = 5
    logging.addLevelName(TRACE_LEVEL, "TRACE")
    logging.TRACE = TRACE_LEVEL
    logging.Logger.trace = trace

    # Logging configuration and setting default level
    logger = logging.getLogger(name)

    _handler = logging.StreamHandler()
    _formatter = logging.Formatter(
        "[%(name)s] %(levelname)s %(asctime)s: %(message)s", "%H:%M:%S"
    )
    _handler.setFormatter(_formatter)

    logger.addHandler(_handler)
    logger.propagate = False

    return logger


def adjust_log_level(logger, verbosity, logfile = None):

    # Determine log level
    if verbosity == 0:
        _formatter = logging.Formatter(
            "[%(name)s] %(asctime)s: %(message)s", "%H:%M:%S"
        )
        level = logging.INFO
    elif int(verbosity) == 1:
        level = logging.DEBUG
    elif int(verbosity) == 2:
        level = logging.TRACE
    elif int(verbosity) != 0:
        logger.warning((
            f"Couldn't interpret verbosity level, it was set to {verbosity}.\n"
            "Ignoring user input, verbosity set to 0!"
            ))

    # Set log level
    logger.setLevel(level)

    if logfile is not None:

        # Reset handlers (idk if this is nescesary)
        while logger.hasHandlers():
            logger.removeHandler(logger.handlers[0])

        # Define where to write logfile
        logging.basicConfig(
            filename=logfile,
            filemode="w",
            format="%(asctime)s - %(levelname)s: %(message)s",
            level=level
        )