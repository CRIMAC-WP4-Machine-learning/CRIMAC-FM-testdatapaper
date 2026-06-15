import logging
import logging.handlers
import sys
from pathlib import Path


class ColorFormatter(logging.Formatter):
    COLORS = {
        logging.DEBUG: "\033[36m",  # cyan
        logging.INFO: "\033[32m",  # green
        logging.WARNING: "\033[33m",  # yellow
        logging.ERROR: "\033[31m",  # red
        logging.CRITICAL: "\033[1;31m",  # bold red
    }
    RESET = "\033[0m"

    def format(self, record):
        color = self.COLORS.get(record.levelno, self.RESET)
        message = super().format(record)
        return f"{color}{message}{self.RESET}"


def setup_logging(
    name: str | None = None,
    log_file: Path | str = "app.log",
    level: int = logging.DEBUG,
):
    """
    Configure logging with:
      - stdout handler (<= INFO, colored)
      - stderr handler (>= WARNING, colored)
      - rotating file handler (DEBUG+, no color)

    DEBUG < INFO < WARNING < ERROR < CRITICAL
    """

    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.propagate = False

    if logger.handlers:
        return logger  # already configured

    log_file = Path(log_file)

    # ---- Base format ----
    fmt = "%(asctime)s | %(levelname)-8s | %(name)-15s | %(message)s"
    datefmt = "%Y-%m-%d %H:%M:%S"

    plain_formatter = logging.Formatter(fmt=fmt, datefmt=datefmt)
    color_formatter = ColorFormatter(fmt=fmt, datefmt=datefmt)

    # ---- stdout handler (INFO and below) ----
    stdout_handler = logging.StreamHandler(sys.stdout)
    stdout_handler.setLevel(logging.INFO)
    stdout_handler.addFilter(lambda r: r.levelno <= logging.INFO)
    stdout_handler.setFormatter(color_formatter)

    # ---- stderr handler (WARNING and above) ----
    stderr_handler = logging.StreamHandler(sys.stderr)
    stderr_handler.setLevel(logging.WARNING)
    stderr_handler.setFormatter(color_formatter)

    # ---- rotating file handler (everything, no color) ----
    file_handler = logging.handlers.RotatingFileHandler(
        log_file,
        maxBytes=10 * 1024 * 1024,  # 10 MB
        backupCount=5,
    )
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(plain_formatter)

    # ---- attach handlers ----
    logger.addHandler(stdout_handler)
    logger.addHandler(stderr_handler)
    logger.addHandler(file_handler)

    return logger
