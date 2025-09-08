"""Test methods for combine_data."""
import os
from pathlib import PosixPath

import numpy as np
import pytest
import tomli

from inftools.analysis.wham import wham
from inftools.tistools.combine_results import combine_data
from inftools.misc.data_helper import data_reader

HERE = PosixPath(__file__).parent

DATD = HERE/"data"
TOMLS=[str(i) for i in [DATD/"sim1.toml", DATD/"sim2.toml", DATD/"sim3.toml"]]
DATAS = [str(i) for i in [DATD/"sim1.txt.gz", DATD/"sim2.txt.gz", DATD/"sim3.txt.gz"]]

def wham_and_assert_results(name, wham_dir):
    # Assert that files exist
    assert os.path.isfile(f"{name}.toml")
    assert os.path.isfile(f"{name}.txt")

    # Get the combo interfaces.
    with open(f"{name}.toml", "rb") as rfile:
        conf = tomli.load(rfile)
        combo_intf = conf["simulation"]["interfaces"]

    intfs = []
    for toml in TOMLS:
        with open(toml, "rb") as rfile:
            conf = tomli.load(rfile)
            intf = conf["simulation"]["interfaces"]
            assert len(combo_intf) >= len(intf)
            intfs += intf
    assert set(intfs) == set(combo_intf)

    # Run wham
    wham(toml=f"{name}.toml", data=f"{name}.txt", nskip = 0, folder=f"{wham_dir}")

    # Assert that important files are created
    assert os.path.isfile(f"{wham_dir}/runav_rate.txt")
    assert os.path.isfile(f"{wham_dir}/runav_flux.txt")
    assert os.path.isfile(f"{wham_dir}/runav_Pcross.txt")

    # Test that we get good results
    data = np.loadtxt(f"{wham_dir}/runav_rate.txt")
    calc_rate = data[-1, -1]/0.025
    kramer = 2.58 * 10 **(-7) # from wf paper
    assert 100*abs(kramer-calc_rate)/kramer < 7 # %

@pytest.mark.heavy
def test_infinit_1(tmp_path: PosixPath) -> None:
    """Test infinit from phase point"""
    data_dir = (HERE / "data").resolve()
    os.chdir(tmp_path)

    combine_data(
        tomls=TOMLS,
        datas=DATAS,
        skip=[100],
        out="combo"
    )

    wham_and_assert_results("combo", "wham")


@pytest.mark.heavy
def test_infinit_2(tmp_path: PosixPath) -> None:
    """Test that combining datas seperately gives the same result as at once.

    With same results we mean the same WHAM results.
    """
    os.chdir(tmp_path)
    # combine all at once
    combine_data(tomls=TOMLS, datas=DATAS, skip=[0], out="combo_all")
    # combine in 2 steps
    combine_data(tomls=TOMLS[:2], datas=DATAS[:2], skip=[0], out="combo0")
    combine_data(
            tomls = ["combo0.toml"] + TOMLS[2:],
            datas = ["combo0.txt"] + DATAS[2:],
            out="combo_sep",
            skip=[0],
    )

    combo_all_txt = data_reader("combo_all.txt")
    combo_sep_txt = data_reader("combo_sep.txt")
    # first check that we have same length in data files
    assert len(combo_all_txt) == len(combo_sep_txt)
    wham_and_assert_results("combo_all", "wham_all")
    wham_and_assert_results("combo_sep", "wham_sep")

    P_all = np.loadtxt(f"wham_all/Pcross.txt")
    P_sep = np.loadtxt(f"wham_sep/Pcross.txt")
    assert np.allclose(P_all,P_sep)
