import os
import logging
from datetime import datetime

DATE_FORMAT = "%m_%d_%y"
DATETIME_FORMAT = "%Y-%m-%d--%H:%M:%S"

def append_date(base: str, time: bool = False) -> str:
    """
    Append today's date to base string.
    """
    today = datetime.now()
    if time:
        date_str = today.strftime(DATETIME_FORMAT)
    else:
        date_str = today.strftime(DATE_FORMAT)
    name = base + "-" + date_str
    return name

def setup_logging():
    logs_dir = "logs"
    if not os.path.exists(logs_dir):
        os.makedirs(logs_dir)
    basename = "onlinekcenter"
    filename = append_date(basename, time=True) + ".log"
    filepath = os.path.join(logs_dir, filename)
    terminal_format = "[%(levelname)s]: %(message)s"
    file_format = "[%(asctime)s - %(name)s - %(levelname)s]: %(message)s"
    logging.basicConfig(
        filename=filepath,
        level=logging.INFO,
        format=file_format,
        datefmt=DATETIME_FORMAT,
    )

    # Create a console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)

    # Create a formatter for the console handler
    console_formatter = logging.Formatter(terminal_format)
    console_handler.setFormatter(console_formatter)

    # Get the root logger and add the console handler to it
    logging.getLogger().addHandler(console_handler)


def gurobi_log_file():
    logs_dir = "logs"
    if not os.path.exists(logs_dir):
        os.makedirs(logs_dir)
    basename = "onlinekcenter_gurobi"
    filename = append_date(basename, time=True) + ".log"
    filepath = os.path.join(logs_dir, filename)
    return filepath
