# -*- coding: utf-8 -*-
"""Tests for selecting per-residue pLDDTs."""
import json
import sys
from pathlib import Path

import pandas as pd
import pytest
import sh

from . import COMMAND
from . import GLOBAL_STATS_FILE
from . import help_check
from . import print_docstring
from . import RESIDUE_FILE
from . import SELECT_INPUTS
from . import SELECT_OUTPUTS
from . import STEM
from . import TOLERANCE

# global constants
SUBCOMMAND = "plddt-select-residues"


def test_subcommand_help():
    """Test subcommand help message."""
    help_check(SUBCOMMAND)


@print_docstring()
def test_select_residues(datadir_mgr):
    """Test select-residues command."""
    with datadir_mgr.in_tmp_dir(
        inpathlist=SELECT_INPUTS,
        save_outputs=True,
        outscope="global",
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
        for filestring in SELECT_OUTPUTS:
            assert Path(filestring).exists()
        df = pd.read_csv(RESIDUE_FILE, sep="\t")
        assert abs(df["pLDDT"][20] - 86.5) <= TOLERANCE
        with Path(GLOBAL_STATS_FILE).open("r") as f:
            json_data = json.loads(f.read())
        assert json_data["usable_residues_pct"]["val"] == 52
