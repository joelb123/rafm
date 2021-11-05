# -*- coding: utf-8 -*-
"""Tests for data ingestion."""
# standard library imports
import sys
from pathlib import Path

import pytest
import sh

from . import help_check
from . import COMMAND
from . import INPUTS
from . import PDB_1
from . import PDB_2
from . import STATS_OUTPUTS
from . import STEM
from . import print_docstring

# global constants
SUBCOMMAND = "plddt-stats"


def test_subcommand_help():
    """Test subcommand help message."""
    help_check(SUBCOMMAND)


@print_docstring()
def test_plddt_stats(datadir_mgr):
    """Test plddt-stats command."""
    with datadir_mgr.in_tmp_dir(
        inpathlist=INPUTS,
        save_outputs=True,
        outscope="module",
    ):
        args = ["--verbose", SUBCOMMAND, "--file-stem", STEM, PDB_1, PDB_2]
        try:
            COMMAND(
                args,
                _out=sys.stderr,
            )
        except sh.ErrorReturnCode as errors:
            print(errors)
            pytest.fail(f"{SUBCOMMAND} failed")
        for filestring in STATS_OUTPUTS:
            assert Path(filestring).exists()
