"""
Logging configuration for the maxiprot package.

This module provides functionality to initialize and configure logging with rich
formatting for better readability in terminal output. It uses the 'rich' library
to create visually enhanced log messages.
"""

import logging
from pathlib import Path
from typing import Optional, Union

from rich.console import Console
from rich.logging import RichHandler

#: Valid values for CLI ``--log-level`` flags.
LOG_LEVELS = ('DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL')


def init_logging(
    loglevel: str = 'INFO', logfile: Optional[Union[str, Path]] = None
) -> None:
    """
    Initialize root logger with specified log level and rich formatting.

    Configures the global logging system with rich formatting for console output
    and optionally writes logs to a file. Log output goes to stderr; data output
    streams (stdout, or files) are never written to by logging.

    Parameters
    ----------
    loglevel : str, optional
        The log level to use (e.g., "DEBUG", "INFO", "WARNING", "ERROR",
        "CRITICAL"), by default "INFO".
    logfile : str or Path, optional
        Path to log file. If provided, logs will be written to both console and file.

    Returns
    -------
    None
        This function configures the global logging system and doesn't return a value.

    Raises
    ------
    ValueError
        If the provided log level is invalid.
    """
    numeric_level = getattr(logging, loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError(f'Invalid log level: {loglevel}')

    root_logger = logging.getLogger()

    # Remove only handlers this function previously installed, so repeated CLI
    # invocations don't stack handlers while host-installed handlers (e.g.
    # pytest's caplog) are left alone.
    for handler in list(root_logger.handlers):
        if getattr(handler, '_maxiprot_handler', False):
            root_logger.removeHandler(handler)

    root_logger.setLevel(numeric_level)

    console_handler = RichHandler(console=Console(stderr=True), show_path=False)
    console_handler.setLevel(numeric_level)
    console_handler._maxiprot_handler = True  # type: ignore[attr-defined]
    root_logger.addHandler(console_handler)

    if logfile:
        try:
            logfile_path = Path(logfile)
            logfile_path.parent.mkdir(parents=True, exist_ok=True)

            file_handler = logging.FileHandler(logfile_path, mode='w', encoding='utf-8')
            file_handler.setLevel(numeric_level)
            file_handler._maxiprot_handler = True  # type: ignore[attr-defined]
            file_handler.setFormatter(
                logging.Formatter(
                    fmt='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                )
            )
            root_logger.addHandler(file_handler)
            logging.debug('Log file: %s', logfile_path.absolute())
        except OSError as e:
            logging.warning('Failed to create log file %s: %s', logfile, e)
            logging.warning('Continuing with console-only logging')

    logging.debug('Logging initialized with level: %s', loglevel)
