from typing import Annotated
import typer

# Disable automatic underscore -> hyphen in CLI names
typer.main.get_command_name = lambda name: name

def calc_flow2(
    plot: Annotated[str, typer.Option("-plot", help="Plot the flow for those paths, string of spaced idxes")]="",
    toml: Annotated[str, typer.Option("-toml", help="The .toml input file defining the orderparameter")] = "infretis.toml",
    log: Annotated[str, typer.Option("-log", help="The .log file to read path numbers")] = "sim.log",
    ):
    import numpy as np
    import tomli
    import matplotlib.pyplot as plt

    with open(toml, "rb") as toml_file:
        config = tomli.load(toml_file)
    n_ensembles = len(config["simulation"]["interfaces"])
    flow_map = {i: {"path": [i], "ens": [i], "step": [0], "cmd": ["init"]} for i in range(n_ensembles)}
    rep_hist = {i: i for i in range(n_ensembles)} # path : replica

    step = 0
    with open(log, "r") as rfile:
        for idx, line in enumerate(rfile):
            if "->" not in line:
                continue
            step += 1
            srip = [i for i in line.rstrip().split() if i.isnumeric()]
            rip = [int(i) for i in srip]

            # accepted 0 <->
            if len(set(srip)) == 6:
                po_m, po_p = rip[2], rip[3]
                pn_m, pn_p = rip[4], rip[5]
                ens_m, ens_p = rip[0], rip[1]

                replica_m = rep_hist[po_m]
                replica_p = rep_hist[po_p]
                # swap
                rep_hist[pn_p] = rep_hist[po_m]
                rep_hist[pn_m] = rep_hist[po_p]

                del rep_hist[po_m]
                del rep_hist[po_p]

                flow_map[replica_m]["path"].append(pn_p)
                flow_map[replica_m]["ens"].append(ens_p)
                flow_map[replica_m]["step"].append(step)
                flow_map[replica_p]["path"].append(pn_m)
                flow_map[replica_p]["ens"].append(ens_m)
                flow_map[replica_p]["step"].append(step)

            # accepted sh/wf
            elif len(set(srip)) == 3:
                ens = rip[0]
                p_o = rip[1]
                p_n = rip[2]

                replica = rep_hist[p_o]
                rep_hist[p_n] = rep_hist[p_o]

                del rep_hist[p_o]

                flow_map[replica]["path"].append(p_n)
                flow_map[replica]["ens"].append(ens)
                flow_map[replica]["step"].append(step)


    if len(plot) > 0:
        pn_numbs = []
        for idx, ens in enumerate([int(i) for i in plot.split(" ")]):
            # for ens in range(n_ensembles):
            if ens == -1:
                ens = n_ensembles - 1
            s1 = len(set(pn_numbs))
            s2 = len(set(flow_map[ens]["path"]))
            s3 = len(set(pn_numbs + flow_map[ens]["path"]))
            assert s1 + s2 == s3

            plt.plot(flow_map[ens]["step"], flow_map[ens]["ens"], marker = "o", markersize = 5, color=f"C{idx%8}", label=f"Replica {ens}")
            plt.xlabel(f"MC steps")
            plt.ylabel(f"Ensemble")
            plt.legend()
            plt.title(f"replica {ens:03.0f}")
        plt.show()

    return flow_map
