import os
import shutil
from pathlib import PosixPath

import numpy as np

from inftools.tistools.shoot import shoot

HERE = PosixPath(__file__).parent

TRAJ = HERE/"data/path/traj.txt"
TOML = HERE/"data/infretis.toml"

ORD0 = np.loadtxt(HERE/"DATA/order.txt")
ENE0 = np.loadtxt(HERE/"DATA/energy.txt")
TRJ0 = np.loadtxt(HERE/"DATA/traj.txt", dtype=str)

def test_shoot_1(tmp_path: PosixPath) -> None:
    """Tests the inftool shooting function with the turtlemd engine.

    This should yield a vaild [0^-] path.
    """
    shutil.copy(TOML, tmp_path)
    os.chdir(tmp_path)

    sht_dir = "shooting_dir0"
    shoot(traj=TRAJ, toml=TOML, name=sht_dir, seed=0)

    ord1 = np.loadtxt(sht_dir + "/0/order.txt")
    ene1 = np.loadtxt(sht_dir + "/0/energy.txt")
    trj1 = np.loadtxt(sht_dir + "/0/traj.txt", dtype=str)

    assert np.sum(ord1[:, 0] - ORD0[:, 0]) < 10e-6
    assert np.sum(ord1[:, 1] - ORD0[:, 1]) < 10e-6
    assert np.sum(ene1[:, 0] - ENE0[:, 0]) < 10e-6
    assert np.sum(ene1[:, 1] - ENE0[:, 1]) < 10e-6
    assert np.sum(ene1[:, 2] - ENE0[:, 2]) < 10e-6
    assert np.sum(ene1[:, 0] - ENE0[:, 0]) < 10e-6
    assert np.all(trj1[:, 2] ==TRJ0[:, 2])
    assert np.all(trj1[:, 3] ==TRJ0[:, 3])
