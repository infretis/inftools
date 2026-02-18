from typer import Option as Opt
from typing import Annotated as And


def plot_error(
    wham: And[str,  Opt("-wham", help="wham folder")] = "wham",
    frac: And[float,  Opt("-frac", help="From what frac of the data to take the average from")] = 0.5,

    ):
    """Plot the wham folder's running error with an estimated error based 
    on the last fraction of the running error.
    """
    import numpy as np
    import matplotlib.pyplot as plt
    import os

    from inftools.misc.calc_help import estimate_error

    fig, axs = plt.subplots(1, 2, figsize=(10,4))

    err_files0 = ["errFLUX.txt", "errPtot.txt", "errRATE.txt"]
    err_files = []
    for file in err_files0:
        path = os.path.join(wham,file)
        err_files.append(path)

    for idx, err_file in enumerate(err_files):
        edata = np.loadtxt(err_file)
        err = estimate_error(edata[:, -1], frac)
        print(f"{err_file:<20}", f"rel-err: {err:.04f}")

        axs[0].plot(edata[:, 0], edata[:, -1], color=f"C{idx}", label=f"{err_files0[idx][:-4]} | {err:.04f}")
        axs[0].axhline(err, color=f"C{idx}")

    edata = np.loadtxt(os.path.join(wham, "errploc.txt"))
    for idx, ens in enumerate(edata.T[1:]):
        err = estimate_error(ens, frac)
        print(f"ploc{idx:<20}", f"rel-err: {err:.04f}")
        axs[1].plot(edata[:, 0], ens, color=f"C{idx}", label=f"{ens} | {err:.04f}")
        axs[1].axhline(err, color=f"C{idx}")

    axs[0].legend()
    plt.show()
