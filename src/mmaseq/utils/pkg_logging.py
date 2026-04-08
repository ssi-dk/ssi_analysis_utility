import logging

def initiate_log(name):
    logger = logging.getLogger(name)

    return logger


def adjust_log(logger, verbosity, logfile = None):
    def trace(self, msg, *args, **kwargs):
        if self.isEnabledFor(TRACE_LEVEL):
            self._log(TRACE_LEVEL, msg, args, **kwargs)


    # Adding TRACE level
    TRACE_LEVEL = 5
    logging.addLevelName(TRACE_LEVEL, "TRACE")
    logging.TRACE = TRACE_LEVEL
    logging.Logger.trace = trace

    # Define text formatting options
    logformat_simple = "[%(name)s] %(asctime)s: %(message)s"
    logformat_advced = "[%(name)s] %(levelname)s %(asctime)s: %(message)s"

    # Determine verbosity settings
    level = logging.INFO
    logformat = logformat_simple
    if int(verbosity) == 1:
        level = logging.DEBUG
        logformat = logformat_advced
    elif int(verbosity) == 2:
        level = logging.TRACE
        logformat = logformat_advced
    elif int(verbosity) != 0:
        print((
            "Warning: "
            f"Couldn't interpret verbosity level, it was set to {verbosity}.\n"
            "Ignoring user input, verbosity set to 0!"
            ))

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

        print(f"Writting messages to logfile: {logfile}")

        # Define where to write logfile
        logging.basicConfig(
            filename=logfile,
            filemode="w",
            format="%(asctime)s - %(levelname)s: %(message)s",
            level=level
        )

    return None
