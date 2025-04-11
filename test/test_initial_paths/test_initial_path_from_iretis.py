import os
import pathlib

import numpy as np
import pytest
from inftools.exercises.puckering import initial_path_from_iretis
import tomli

HERE = pathlib.Path(__file__).resolve().parent

def test_similar_interfaces(tmp_path: pathlib.PosixPath):
    """Test that we get pack same paths if using interfaces that
    are similar to the previous ones.
    """
    out0 = tmp_path / "test0"
    initial_path_from_iretis(
            traj = "data/run0",
            toml = "data/infretis0.toml",
            restart = "data/restart0.toml",
            keep_all_active = True,
            out_dir = out0,
            )
    for orderf in out0.glob("*/order.txt"):
        x0 = np.loadtxt(HERE/"data"/"run0"/str(orderf.parts[-2])/"order.txt")
        x = np.loadtxt(orderf)
        assert np.allclose(x0,x)


def test_add_active_paths(tmp_path: pathlib.PosixPath):
    """Test that we can add interfaces if active paths are not valid."""
    out0 = tmp_path / "test1"
    initial_path_from_iretis(
            traj = "data/run0",
            toml = "data/infretis1.toml",
            restart = "data/restart0.toml",
            keep_all_active = True,
            out_dir = out0,
            out_toml = tmp_path / "new_intf.toml"
            )
    with open(tmp_path / "new_intf.toml", "rb") as rfile:
        config = tomli.load(rfile)
    intfs = config["simulation"]["interfaces"]
    sh_m = config["simulation"]["shooting_moves"]
    assert len(intfs) == len(sh_m)

    # check that plus paths are valid above intfs
    for i,intf in enumerate(intfs[:-1]):
        orderf = out0 / str(i+1) / "order.txt"
        x = np.loadtxt(orderf)
        assert np.max(x[:,1]) > intf


