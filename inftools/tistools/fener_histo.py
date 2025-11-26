from typing import Annotated as Atd
from typer import Option as Opt

def plot_hist(
    xval: Atd[str, Opt(help="xval histogram file")] = "histo_xval.txt",
    yval: Atd[str, Opt(help="yval histogram file")] = "histo_yval.txt",
    hist: Atd[str, Opt(help = "the histogram or free-energy file")] = "histo_free_energy.txt",
    ):
    """
    Plot the output from the WHAM free-energy analysis.
    """
    import numpy as np
    import matplotlib.pyplot as plt

    x = np.loadtxt(xval)
    y = np.loadtxt(yval)
    h = np.loadtxt(hist)
    plt.pcolormesh(x, y, h)
    plt.show()
