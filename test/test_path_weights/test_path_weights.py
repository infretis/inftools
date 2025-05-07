import pytest
import numpy as np
import pathlib

from inftools.analysis.wham import wham
from inftools.tistools.path_weights import get_path_weights
from inftools.misc.data_helper import data_reader

HERE = pathlib.Path(__file__).parent

def test_get_path_weights_sh(tmp_path: pathlib.PosixPath) -> None:
    """Test that infretis_data.txt with sh returns correct weights and pcross.

    We have 10 entries in the [0+] ensembles and 1 in the [0-]. The [0+] paths
    have order parameter 1.25 x 2 and 1.3 x 8, so we should get 100% cross 1.25
    and 80% *cross* 1.3 (paths count as crossing if orderp >= x). The interfaces
    are [1.0, 2.0], so the pcross should read

        1.00    1.0
        1.25    1.0
        1.30    0.8

    and the expected path weights should all be 0.1 in magnitude.
    """
    get_path_weights(
            toml = HERE / "data/infretis.toml",
            data = HERE / "data/infretis_data.txt",
            out = tmp_path / "path_weights.txt",
            outP = tmp_path / "pcross.txt",
            )
    pcross = np.loadtxt(tmp_path / "pcross.txt")
    path_weights = np.loadtxt(tmp_path / "path_weights.txt")

    expected_pcross = [[1.0, 1.0],[1.25, 1.0],[1.30, 0.8]]
    expected_path_weights = [0.1 for i in range(10)]

    assert np.allclose(pcross, expected_pcross)
    assert np.allclose(path_weights[:,2], expected_path_weights)
