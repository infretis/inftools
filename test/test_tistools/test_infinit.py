"""Test methods for doing TIS."""
import difflib
import filecmp
import os
import shutil
from pathlib import PosixPath
from subprocess import STDOUT, check_output

import pytest
import tomli
import tomli_w

from inftools.misc.infinit_helper import read_toml


@pytest.mark.heavy
def test_infinit_1(tmp_path: PosixPath) -> None:
    """Test infinit from phase point"""
    folder = tmp_path / "temp"
    folder.mkdir()
    basepath = PosixPath(__file__).parent
    load_dir = (
        basepath / "../../examples/turtlemd/double_well/load_copy"
    ).resolve()
    toml_dir = basepath / "data/infretis.toml"
    conf_dir = basepath / "data/initial.xyz"
    shutil.copy(str(toml_dir), str(folder))
    shutil.copy(str(conf_dir), str(folder))
    os.chdir(folder)

    # Run infinit
    os.system("inft infinit")

    # Check if files exist and that p is same or growing
    tomls = [f"infretis_{i}.toml" for i in range(2, 6)] + ["infretis.toml"]
    assert os.path.isfile(tomls[0])

    config0 = read_toml(tomls[0])
    assert "infinit" in config0
    for i,toml in enumerate(tomls):
        assert os.path.isfile(toml)
        assert os.path.isfile(f"combo_{i}.toml")
        assert os.path.isfile(f"combo_{i}.txt")
        # number of lines in combined infretis_data
        num_lines = num_lines = sum(1 for _ in open(f"combo_{i}.txt"))
        # check that we actually have data
        assert num_lines > 0
        # check that we have more data than previous run
        if i > 1:
            assert num_lines > prev_run_num_lines
        prev_run_num_lines = num_lines
