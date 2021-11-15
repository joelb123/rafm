# -*- coding: utf-8 -*-
"""Tests for calculating stats."""
import json
import sys
from pathlib import Path

import pandas as pd
import pytest
import sh

from . import COMMAND
from . import GLOBAL_STATS_FILE
from . import help_check
from . import MODEL_FILE_1_NAME
from . import MODEL_FILE_2_NAME
from . import print_docstring
from . import STATS_FILE
from . import STATS_INPUTS
from . import STATS_OUTPUTS
from . import STEM
from . import TOLERANCE

# global constants
SUBCOMMAND = "plddt-stats"


def test_subcommand_help():
    """Test subcommand help message."""
    help_check(SUBCOMMAND)


@print_docstring()
def test_plddt_stats(datadir_mgr):
    """Test plddt-stats command."""
    with datadir_mgr.in_tmp_dir(
        inpathlist=STATS_INPUTS,
        save_outputs=True,
        outscope="global",
    ):
        args = [
            "--verbose",
            SUBCOMMAND,
            "--file-stem",
            STEM,
            MODEL_FILE_1_NAME,
            MODEL_FILE_2_NAME,
        ]
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
        df = pd.read_csv(STATS_FILE, sep="\t")
        df = df.set_index("file")
        stat_vals_2 = {
            "residues_in_pLDDT": 21,
            "pLDDT_mean": 90.55,
            "pLDDT_median": 91.82,
            "pLDDT80_count": 21,
            "pLDDT80_frac": 1.0,
            "pLDDT80_mean": 90.55,
            "pLDDT80_median": 91.82,
            "LDDT_expect": 0.814,
        }
        for key in stat_vals_2:
            stat_val = stat_vals_2[key]
            if isinstance(stat_val, float):
                assert abs(df[key][MODEL_FILE_2_NAME] - stat_val) <= TOLERANCE
            else:
                assert df[key][MODEL_FILE_2_NAME] == stat_val
        assert df["passing"][MODEL_FILE_2_NAME]
        global_stats = {
            "models_in": 2,
            "min_length": 20,
            "min_count": 20,
            "plddt_lower_bound": 80,
            "plddt_upper_bound": 100,
            "plddt_criterion": 91.2,
            "total_residues": 41,
            "models_selected": 1,
            "model_selection_pct": 50,
            "selected_residues": 21,
        }
        with Path(GLOBAL_STATS_FILE).open("r") as f:
            json_data = json.loads(f.read())
        for key in global_stats:
            stat_val = global_stats[key]
            if isinstance(stat_val, float):
                assert abs(json_data[key]["val"] - stat_val) <= TOLERANCE
            else:
                assert json_data[key]["val"] == stat_val
