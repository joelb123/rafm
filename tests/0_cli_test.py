# -*- coding: utf-8 -*-
"""Tests for basic CLI function."""
# third-party imports
import pytest
import sh

from . import COMMAND
from . import help_check
from . import print_docstring
from . import working_directory


def test_cli():
    """Test global help function."""
    help_check("global")


@print_docstring()
def test_version(tmp_path):
    """Test version command."""
    with working_directory(tmp_path):
        try:
            output = COMMAND(["--version"])
        except sh.ErrorReturnCode as errors:
            print(errors)
            pytest.fail(errors)
        assert "version" in output
