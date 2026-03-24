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
    logger.setLevel(logging.INFO)

    _handler = logging.StreamHandler()
    _formatter = logging.Formatter("[%(levelname)s] %(name)s: %(message)s")
    _handler.setFormatter(_formatter)

    logger.addHandler(_handler)
    logger.propagate = False

    return logger


def adjust_log_level(logger, verbosity):
    if int(verbosity) == 0:
        return None
    elif int(verbosity) == 1:
        logger.setLevel(logging.DEBUG)
    elif int(verbosity) == 2:
        logger.setLevel(5)
    else:
        logger.warning((
            f"Couldn't interpret verbosity level, it was set to {verbosity}.\n"
            "Ignoring user input, verbosity set to 0!"
            ))

    return None