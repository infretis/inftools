from typing import Annotated
import typer


def plot_ufree_2d(
    whams: Annotated[
        list[str],
        typer.Option("-whams", help="wham folder for the two simulations"),
    ],
    out: Annotated[str, typer.Option("-out", help="name for output .txt/toml file.")] = "",
):
    """
    Given wham folders containing A -> B and B -> A
    and their 2D calculated conditional free energies,
    calculate the merged 2D unconditional free energy.

    The wham folders must be in the order of 1. A -> B and 2. B -> A.
    The histo_xval and histo_yval must be the same between the folders.

    * P. G. Bolhuis, J. Chem. Phys. 129, 114108 (2008)
    """
    import matplotlib.pyplot as plt
    import numpy as np
    from scipy.integrate import simpson

    hists = []
    hints = []
    rates = []
    txts = ["histo_probability.txt", "runav_rate.txt"]

    for wham in whams:
        xval = np.loadtxt(f"{wham}/histo_xval.txt")
        yval = np.loadtxt(f"{wham}/histo_yval.txt")
        hists.append(np.loadtxt(f"{wham}/histo_probability.txt"))
        hints.append(simpson(simpson(hists[-1], xval), yval))  # 2d integral
        rates.append(np.loadtxt(f"{wham}/runav_rate.txt")[-1, 3])

    # From detailed balance relation:
    # p_a x k_ab = p_b x k_ba -> p_a x k_ab = c x p_b* x k_ba
    # we have p_a, p_b*, k_ab, k_ba. calculate c:
    k_ratio = rates[0] / rates[1]
    c = (hints[0] / hints[1]) * k_ratio

    # calcuate reweighted hist_rw
    frees = [-np.log(hists[0]), -np.log(hists[1] * c)]
    hist_rw = np.zeros(hists[0].shape)
    for idx_x, x in enumerate(xval):
        for idx_y, y in enumerate(xval):
            # get the x,y values of the frees in a list
            free_vals = np.array([free[idx_x, idx_y] for free in frees])

            # prune infinities
            inf = free_vals == np.inf
            not_inf = np.logical_not(inf)
            if np.all(inf):
                hist_rw[idx_x, idx_y] = np.inf
            else:
                hist_list = free_vals[not_inf == 1]
                exp_list = -np.log(np.sum([np.exp(-i) for i in hist_list]))
                hist_rw[idx_x, idx_y] = exp_list

    bar = plt.pcolormesh(yval, xval, hist_rw, rasterized=True) #, cmap="magma"
    if len(out) > 0:
        plt.savefig(out)

    plt.show()
