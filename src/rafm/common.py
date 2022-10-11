# -*- coding: utf-8 -*-
"""Global constants and CLI functions in common."""
# standard library imports
from __future__ import annotations

import sys
from typing import TypedDict

import loguru
import typer
from loguru import logger
from statsdict import StatsDict

from . import __doc__ as docstring
from . import NAME

# third-party imports
# module imports

# global constants
DEFAULT_STDERR_LOG_LEVEL = "INFO"
NO_LEVEL_BELOW = 30  # Don't print level for messages below this level


class GlobalState(TypedDict):
    """Dictionary of global state variables."""

    verbose: bool
    log_level: str


STATE: GlobalState = {"verbose": False, "log_level": DEFAULT_STDERR_LOG_LEVEL}


def _stderr_format_func(record: loguru.Record) -> str:
    """Do level-sensitive formatting."""
    if record["level"].no < NO_LEVEL_BELOW:
        return "<level>{message}</level>\n"
    return "<level>{level}</level>: <level>{message}</level>\n"


logger.remove()
logger.add(sys.stderr, level=STATE["log_level"], format=_stderr_format_func)
APP = typer.Typer(help=docstring, name=NAME)
STATS = StatsDict(logger=logger, app=APP, module_name=NAME)
