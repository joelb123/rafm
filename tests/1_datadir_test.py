# -*- coding: utf-8 -*-
"""Create a clean directory for test data."""
# standard library imports
import shutil
from pathlib import Path

from . import ALPHAFOLD_SITE
from . import MODEL_FILE_3_NAME
from . import print_docstring


@print_docstring()
def test_clean_datadir(request):
    """Clean up datadir."""
    testdir = Path(request.fspath.dirpath())
    datadir = testdir / "data"
    if datadir.exists():
        shutil.rmtree(datadir)  # remove anything left in data directory


@print_docstring()
def test_setup_datadir(request, datadir_mgr, capsys):
    """Copy in and download static data."""
    testdir = Path(request.fspath.dirpath())
    datadir = testdir / "data"
    filesdir = testdir / "testdata"
    shutil.copytree(filesdir, datadir)
    with capsys.disabled():
        datadir_mgr.download(
            download_url=ALPHAFOLD_SITE,
            files=[MODEL_FILE_3_NAME],
            scope="global",
            md5_check=False,
            progressbar=True,
        )
