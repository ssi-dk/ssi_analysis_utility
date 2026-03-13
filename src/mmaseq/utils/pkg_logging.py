import logging

def initiate_log(name):
    # Logging configuration
    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)

    _handler = logging.StreamHandler()
    _formatter = logging.Formatter("[%(levelname)s] %(name)s: %(message)s")
    _handler.setFormatter(_formatter)

    logger.addHandler(_handler)
    logger.propagate = False

    return logger


def adjust_log_level(logger, debug):
    if debug:
        logger.setLevel(logging.DEBUG)

    return None