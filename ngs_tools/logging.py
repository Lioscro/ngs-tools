import functools
import logging
from contextlib import contextmanager
from typing import Optional


def silence_logger(name: str):
    """Given a logger name, silence it completely.

    Args:
        name: Logger name to silence
    """
    package_logger = logging.getLogger(name)
    package_logger.setLevel(logging.CRITICAL + 10)
    package_logger.propagate = False


class Logger:
    """Custom logger that provides namespacing functionality.

    The following loggers are silenced by default:
        anndata, h5py, numba, pysam, pystan

    Attributes:
        FORMAT: Static attribute that encodes the logging format string
        logger: The logger object
        ch: A :class:`logging.Streamhandler` object that sets the format
    """
    FORMAT = '[%(asctime)s] %(levelname)7s %(message)s'

    def __init__(self, name: Optional[str] = None):
        self.namespace = 'main'

        self.logger = logging.getLogger(name or __name__)
        self.ch = logging.StreamHandler()
        self.ch.setFormatter(logging.Formatter(self.FORMAT))
        self.logger.addHandler(self.ch)
        self.logger.propagate = False

        # Other global initialization
        silence_logger('anndata')
        silence_logger('h5py')
        silence_logger('numba')
        silence_logger('pysam')
        silence_logger('pystan')

    def namespaced(self, namespace: str):
        """Function decorator to set the logging namespace for the duration of
        the function.

        Args:
            namespace: The namespace
        """

        def wrapper(func):

            @functools.wraps(func)
            def inner(*args, **kwargs):
                previous = self.namespace
                try:
                    self.namespace = namespace
                    return func(*args, **kwargs)
                finally:
                    self.namespace = previous

            return inner

        return wrapper

    @contextmanager
    def namespaced_context(self, namespace: str):
        """Context manager to set the logging namespace.

        Args:
            namespace: The namespace
        """
        previous = self.namespace
        self.namespace = namespace
        yield
        self.namespace = previous

    def namespace_message(self, message: str) -> str:
        """Add namespace information at the beginning of the logging message.

        Args:
            message: The logging message

        Returns:
            The namespaced message
        """
        return f'[{self.namespace}] {message}'

    def addHandler(self, hdlr: logging.Handler, format: bool = True):
        if format:
            hdlr.setFormatter(logging.Formatter(self.FORMAT))
        return self.logger.addHandler(hdlr)

    def removeHandler(self, *args, **kwargs):
        return self.logger.removeHandler(*args, **kwargs)

    def setLevel(self, *args, **kwargs):
        return self.logger.setLevel(*args, **kwargs)

    def debug(self, message, *args, **kwargs):
        return self.logger.debug(
            self.namespace_message(message), *args, **kwargs
        )

    def info(self, message, *args, **kwargs):
        return self.logger.info(
            self.namespace_message(message), *args, **kwargs
        )

    def warning(self, message, *args, **kwargs):
        return self.logger.warning(
            self.namespace_message(message), *args, **kwargs
        )

    def exception(self, message, *args, **kwargs):
        return self.logger.exception(
            self.namespace_message(message), *args, **kwargs
        )

    def critical(self, message, *args, **kwargs):
        return self.logger.critical(
            self.namespace_message(message), *args, **kwargs
        )

    def error(self, message, *args, **kwargs):
        return self.logger.error(
            self.namespace_message(message), *args, **kwargs
        )


logger = Logger()


def set_logger(log: Logger):
    """Set the logger to the provided :class:`Logger` instance. Use this function
    to override the default logger for this (ngs-tools) library from libraries that
    use this library (ngs-tools) as a dependency.

    Args:
        log: The logger
    """
    global logger
    logger = log
