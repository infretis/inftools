import pytest
import numpy as np
import pathlib

from inftools.misc.infinit_helper import estimate_interface_positions

HERE = pathlib.Path(__file__).parent


# nice orderp values for testing
X = np.array([i for i in range(10)])
# bit harder with closely spaced x (but not duplicate values)
# inftools.tistools.estimate_interface_positions does not return duplicates
X2 = np.array([0] + [1 + 1e-4*i for i in range(4)] +  [i+2 for i in range(5)])

# local crossing probability of 0.5
PL0 = 0.5
# pcross for first test
P0 = np.array([PL0**i for i in range(10)])
# pcross for second test with a large drop
P1 = np.array([1, 1, 1] + [PL0**(i+3) for i in range(7)])

@pytest.mark.parametrize("x,p,pL", [
    (X, P0, 0.50),
    (X, P0, 0.49),
    (X, P0, 0.51),
    (X, P1, 0.55),
    (X, P1, 0.25),
    (X2,P1, 0.79),
    (X2,P1, 0.19),
    ])
def test_estimate_binless_pcross(x, p, pL):
    """Test that we can place interfaces exactly for some test cases.

    NOTE:
    It is assumed that x does not contain duplicates, as returned by
    `inft get_path_weights -outP pcross.txt`.
    """
    # last interface is not added
    interfaces, pL_used = estimate_interface_positions(x, p, pL)
    expected_ptot_at_interfaces = [pL_used**(i) for i in range(len(interfaces))]
    # get actual Ptot at interfaces by interpolating X vs log(P)
    ptot_at_interfaces = np.exp(np.interp(interfaces, x, np.log(p)))
    assert np.allclose(ptot_at_interfaces, expected_ptot_at_interfaces)
    # uncomment to see plots of pcross and estimated interfacess
    #import matplotlib.pyplot as plt
    #plt.plot(x, p, marker="x",lw=2.5)
    #for i, (intf,ptot) in enumerate(zip(interfaces, expected_ptot_at_interfaces)):
    #    plt.axvline(intf,c="k")
    #    plt.axhline(ptot,c="k")
    #plt.yscale("log")
    #plt.plot(interfaces, ptot_at_interfaces)
    #plt.plot(interfaces, expected_ptot_at_interfaces)
    #print("expected Ptot at interfaces", expected_ptot_at_interfaces)
    #plt.show()
