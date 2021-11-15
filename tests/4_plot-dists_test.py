# -*- coding: utf-8 -*-
"""Tests for plotting distributions."""
import sys
from pathlib import Path

import pytest
import sh

from . import COMMAND
from . import help_check
from . import IMAGE_FILE
from . import PLOT_INPUTS
from . import print_docstring
from . import STEM

# global constants
SUBCOMMAND = "plddt-plot-dists"


def test_subcommand_help():
    """Test subcommand help message."""
    help_check(SUBCOMMAND)


@print_docstring()
def test_plot_dists(datadir_mgr):
    """Test plddt-plot-dists command."""
    with datadir_mgr.in_tmp_dir(
        inpathlist=PLOT_INPUTS,
        save_outputs=True,
        outscope="module",
    ):
        args = [
            "--verbose",
            SUBCOMMAND,
            "--file-stem",
            STEM,
        ]
        try:
            COMMAND(
                args,
                _out=sys.stderr,
            )
        except sh.ErrorReturnCode as errors:
            print(errors)
            pytest.fail(f"{SUBCOMMAND} failed")
        assert Path(IMAGE_FILE).exists()
