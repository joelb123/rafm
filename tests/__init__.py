"""Base for pytest testing."""
# standard library imports
import contextlib
import functools
import os
from pathlib import Path
from typing import Callable

import pytest
import sh
from sh import ErrorReturnCode

# global constants
PROGRAM_NAME = "rafm"
ALPHAFOLD_SITE = "https://alphafold.ebi.ac.uk/files/"
MODEL_FILE_1_NAME = "AF-A0A075B6Y3-F1-model_v1.pdb"
MODEL_FILE_2_NAME = "AF-Q9NRI7-F1-model_v1.pdb"
MODEL_FILE_3_NAME = "AF-P0DP23-F1-model_v1.cif"
STATS_INPUTS = [MODEL_FILE_1_NAME, MODEL_FILE_2_NAME, MODEL_FILE_3_NAME]
STEM = "mytest"
STATS_FILE = f"{STEM}_plddt_stats.tsv"
GLOBAL_STATS_FILE = f"{PROGRAM_NAME}_stats.json"
SELECT_INPUTS = [MODEL_FILE_2_NAME, STATS_FILE]
STATS_OUTPUTS = [STATS_FILE, GLOBAL_STATS_FILE]
RESIDUE_FILE = f"{STEM}_plddt80_91.2.tsv"
SELECT_OUTPUTS = [RESIDUE_FILE, GLOBAL_STATS_FILE]
PLOT_INPUTS = [STATS_FILE, RESIDUE_FILE]
IMAGE_FILE = f"{STEM}_dists.png"
COMMAND = sh.Command(PROGRAM_NAME)
TOLERANCE = 0.01


@contextlib.contextmanager
def working_directory(path: str) -> None:
    """Change working directory in context."""
    prev_cwd = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(prev_cwd)


def help_check(subcommand: str) -> None:
    """Test help function for subcommand."""
    print(f"Test {subcommand} help.")
    if subcommand == "global":
        help_command = ["--help"]
    else:
        help_command = [subcommand, "--help"]
    try:
        output = sh.rafm(help_command)
    except ErrorReturnCode as errors:
        print(errors)
        pytest.fail(f"{subcommand} help test failed")
    print(output)
    assert "Usage:" in output
    assert "Options:" in output


def print_docstring() -> Callable:
    """Decorator to print a docstring."""

    def decorator(func: Callable) -> Callable:
        """Define decorator."""

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            """Print docstring and call function."""
            print(func.__doc__)
            return func(*args, **kwargs)

        return wrapper

    return decorator


def run_command(args, component):
    """Run command with args."""
    command_string = " ".join(args)
    print(f"Testing {component} with" + f'"{PROGRAM_NAME} {command_string}"')
    try:
        COMMAND(args)
    except ErrorReturnCode as errors:
        print(errors)
        pytest.fail(f"{component} failed")
