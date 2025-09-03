from typing import Annotated
import typer


def calc_flow_op(
    rep: Annotated[int, typer.Option("-rep", help="the replica to plot")]= 0,
    toml: Annotated[str, typer.Option("-toml", help="The .toml input file defining the orderparameter")] = "infretis.toml",
    log: Annotated[str, typer.Option("-log", help="The .log file to read path numbers")] = "sim.log",
    load: Annotated[str, typer.Option("-load", help="The trajectory load folder containing order.txt")] = "load",
    ):
    """
    Plots the order.txt for a single replica as it visits different ensembles.

    The x axis becomes the total number of subcycles propagated for the replica.
    (if we ignore the rejected moves.)
    """
    import matplotlib.pyplot as plt
    import numpy as np

    from inftools.tistools.flow2 import calc_flow2
    from inftools.misc.infinit_helper import read_toml

    intfs = read_toml(toml)["simulation"]["interfaces"]
    # fig, axs = plt.subplots(1, 2, figsize=(10.5, 4.0))
    fig, axs = plt.subplots(1, 2, figsize=(10.5, 4.0))
    for intf in intfs:
        axs[1].axhline(intf, color="k", alpha=0.2)

    flow_map = calc_flow2(toml=toml, log=log)
    reps = np.arange(len(intfs))
    xacc = 0
    axs[0].plot(flow_map[reps[rep]]["step"], flow_map[reps[rep]]["ens"], marker = "o", markersize = 5, color="C0", label=f"Replica {reps[rep]}")
    axs[0].set_ylim([0, len(intfs)])
    axs[0].set_xlabel("Cycle")
    axs[0].set_ylabel("Ensemble")
    for idx, path in enumerate(flow_map[reps[rep]]["path"][1:]):
        ens = flow_map[reps[rep]]["ens"][1:][idx]
        try:
            data = np.loadtxt(load + f"/{path}/order.txt")
            axs[1].plot(data[:, 0] + xacc, data[:, 1], color=f"C{ens%8}")
            xacc += len(data[:, 0])
        except:
            continue
    axs[1].set_xlabel("Subcycles")
    axs[1].set_ylabel("Order Parameter")
    plt.show()

