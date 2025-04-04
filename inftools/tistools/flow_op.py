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
    for intf in intfs:
        plt.axhline(intf, color="k", alpha=0.2)

    flow_map = calc_flow2(toml=toml, log=log)
    xacc = 0
    for idx, path in enumerate(flow_map[rep]["path"][1:]):
        ens = flow_map[rep]["ens"][1:][idx]
        print(load + f"/{path}/order.txt")
        data = np.loadtxt(load + f"/{path}/order.txt")
        plt.plot(data[:, 0] + xacc, data[:, 1], color=f"C{ens%8}")
        xacc += len(data[:, 0])
    plt.xlabel("Subcycles")
    plt.ylabel("Order Parameter")
    plt.show()

