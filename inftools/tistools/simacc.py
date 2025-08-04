from typing import Annotated
import typer


def sim_acc(
    log: Annotated[str, typer.Option("-log")] = "sim.log",
    toml: Annotated[str, typer.Option("-toml")] = "infretis.toml",
    ):
    """Calculates and plots the acceptance rate of MC shooting moves
    for individual ensembles.
    """
    import os
    import matplotlib.pyplot as plt
    from inftools.misc.infinit_helper import read_toml
    import numpy as np
    enss = {}
    with open(log, "r") as read:
        for idx, line in enumerate(read):
            if "shooted" in line and "sh sh" not in line:
                rip = line.rstrip().split()
                ens = int(rip[5])
                if ens not in enss:
                    enss[ens] = []
                status = read.readline().split()
                enss[ens].append(status[3] == "ACC")
                if status[3] != "ACC" and ens not in (0, 1):
                    print(f"Failed move in ensemble {ens} with path {rip[8]} and status {status[3]}")

    keys = list(sorted(enss.keys()))
    pacc = [sum(enss[key])/len(enss[key]) for key in keys]

    if os.path.exists(toml):
        intfs = read_toml(toml)["simulation"]["interfaces"]
        x = [intfs[0]] + intfs[:-1]
        label = "Order parameter"
    else:
        x = keys
        label = "Ensemble"

    plt.scatter(x[0], np.array(pacc[0])*100,   color="C0", label="0-")
    plt.scatter(x[1:], np.array(pacc[1:])*100, color="C1", label="i+")
    plt.xlabel(label)
    plt.ylabel("Acceptance probability [%]")
    plt.legend()
    plt.ylim([0, 105])
    plt.show()
