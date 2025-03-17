from typing import Annotated
import typer

# Disable automatic underscore -> hyphen in CLI names
typer.main.get_command_name = lambda name: name

def calc_flow2(
    toml: Annotated[str, typer.Option("-toml", help="The .toml input file defining the orderparameter")] = "infretis.toml",
    log: Annotated[str, typer.Option("-log", help="The .log file to read path numbers")] = "sim.log",
    out: Annotated[bool, typer.Option("-out", help="The output of the analysis")] = False,
    ):
    import numpy as np
    import tomli
    import matplotlib.pyplot as plt

    with open(toml, "rb") as toml_file:
        config = tomli.load(toml_file)
    n_ensembles = len(config["simulation"]["interfaces"])
    flow_map = {i: {"path": [i], "ens": [i], "step": [0]} for i in range(n_ensembles)}
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
                print(
                po_m, po_p = rip[10], rip[11]
                pn_m, pn_p = rip[13], rip[14]

                assert flow_map[0]["path"][-1] == po_m
                assert flow_map[1]["path"][-1] == po_m
                assert flow_map[0]["ens"][-1] == 1
                assert flow_map[1]["ens"][-1] == 0

                flow_map[0]["path"].append(rip[13])
                flow_map[0]["ens"].append(0)
                flow_map[1]["path"].append(rip[14])
                flow_map[1]["ens"].append(1)

            # accepted sh/wf
            elif len(set(srip)) == 3:
                replica = rep_hist[rip[1]]
                rep_hist[rip[2]] = rep_hist[rip[1]]
                if rip[1] in rep_hist:
                    del rep_hist[rip[1]]

                print("hehe", srip)
                print(replica, rep_hist)
                flow_map[replica]["path"].append(rip[2])
                flow_map[replica]["ens"].append(rip[0])
                flow_map[replica]["step"].append(step)

                # path = rep_hist[
                print("b", int(rip[0]), int(rip[1]), int(rip[2]))
                print(line)
                print(flow_map[replica])
                print("")
            # if step == 7:

            #     print("bear")
            #     exit()

