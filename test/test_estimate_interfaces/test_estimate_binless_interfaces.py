import numpy as np
import pathlib

from inftools.misc.infinit_helper import estimate_binless_interface_positions

HERE = pathlib.Path(__file__).parent


# orderp values for testing
X = np.array([i for i in range(10)])
X1 = np.array([0] + [1 for i in range(4)] +  [i+2 for i in range(5)])

# local crossing probability of 0.5
PL0 = 0.5
# pcross for first test
P0 = np.array([PL0**i for i in range(10)])

# pcross for second test with large drop
P1 = np.array([1, 1, 1] + [PL0**(i+3) for i in range(7)])

def test_estimate_binless_pcross():
    """Test that we can place interfaces in an exact test case"""

    # contains orderp and pcross with large drop
    interfaces = estimate_binless_interface_positions(X, P0, 0.5)
    # uncomment to see plots of pcross and estimated interfacess
    import matplotlib.pyplot as plt
    plt.plot(X, P0, marker="x",lw=2.5)
    for i, intf in enumerate(interfaces):
        plt.axvline(intf,c="k")
        plt.axhline(PL0**(i+1),c="k")
    plt.yscale("log")
    plt.show()
    assert np.allclose(X, interfaces)


def test_large_drop_binless_pcross():
    """Test that we place interfaces correctly with a large drop in pcross."""
    # contains orderp and pcross with large drop
    p_here=0.50
    interfaces = estimate_binless_interface_positions(X, P1, p_here)
    print(interfaces)
    import matplotlib.pyplot as plt
    plt.plot(X, P1, marker="x",lw=2.5)
    for i, intf in enumerate(interfaces):
        plt.axvline(intf,c="k")
        plt.axhline(p_here**(i+1),c="k")
    plt.yscale("log")
    plt.show()
    #assert np.allclose(X, interfaces)

#def test_multiple_points_binless_pcross():


test_large_drop_binless_pcross()
#test_estimate_binless_pcross()
