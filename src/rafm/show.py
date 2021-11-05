# -*- coding: utf-8 -*-

from pathlib import Path

import streamlit as st
import typer
from stmol import component_3dmol

# from .common import APP

# @APP.command()
def show(pdb_path: Path = typer.Argument(...)):
    """Show AlphaFold structure with reliability measures."""
    st.title(f"{pdb_path}")
    component_3dmol()


if __name__ == "__main__":
    test_file = "tests/testdata/AF-A0A075B6Y3-F1-model_v1.pdb"
    show(test_file)
