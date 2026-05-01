from typing import Annotated
import typer

# Disable automatic underscore -> hyphen in CLI names
# typer.main.get_command_name = lambda name: name

def calc_iflow(
    plot: Annotated[str, typer.Option("-plot", help="Plot the flow for those paths, string of spaced idxes")]="",
    log: Annotated[str, typer.Option("-log", help="The .log file to read path numbers")] = "sim.log",
    ):
    """
    Calculates and plots the flow of individual replica across ensembles.

    Returns a flow_map dictionary.
    """
    import numpy as np
    import tomli
    import matplotlib.pyplot as plt

    from inftools.misc.misc import read_log
    from inftools.tistools.max_op import COLS


    for idx0, rep in enumerate([int(i) for i in plot.split(" ")]):
    
        ens, pns, follow_size, shootings = read_log("sim.log", rep)
        print("ens:", " ".join(ens))
        print("pns:", " ".join(pns))
    
        plt.plot(shootings, [int(i) for i in ens], color=COLS[idx0])
        plt.scatter(shootings, [int(i) for i in ens], color=COLS[idx0])
        if idx0 == 0:
            acc = 50
            for idx, fs in enumerate(follow_size):
                plt.plot([acc*idx, acc*(idx+1)], [fs[0]-1]*2, color="k")
            # acc += fs[1]
            # acc += 50
    
    
    plt.xlabel("MC Moves")
    plt.ylabel("Ensemble")
    plt.ylim([0, None])
    plt.savefig("iflow.png")
    plt.show()
