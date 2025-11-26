from typing import Annotated
import typer

# Disable automatic underscore -> hyphen in CLI names
# typer.main.get_command_name = lambda name: name

def calc_flow(
    plot: Annotated[str, typer.Option("-plot", help="Plot the flow for those paths, string of spaced idxes")],
    toml: Annotated[str, typer.Option("-toml", help="The .toml input file defining the orderparameter")] = "infretis.toml",
    log: Annotated[str, typer.Option("-log", help="The .log file to read path numbers")] = "sim.log",
    out: Annotated[str, typer.Option("-out", help="The output of the analysis")] = "",
    ):
    """
    Keep track of which parent paths are in which ensemble after each MC move.

    We start the simulation with a valid path in each ensemble, and we
    would like to see how the initial paths move through the ensembles.

    For example, when starting out, path nr. 4 is initially in ensemble 004.
    However, it can be swapped with other ensembles, and as such move upward
    all the way til it become a reactive path, or it may move down all the way
    til it becomes a path in [0-].

    The output of this function returns a text file where the column i specifies
    which ensemble the parent path i is currently in.
    """
    import numpy as np
    import tomli
    import matplotlib.pyplot as plt

    plot = [int(i) for i in plot.split(' ')]

    with open(toml, "rb") as toml_file:
        config = tomli.load(toml_file)

    n_ensembles = len(config["simulation"]["interfaces"])

    # map a path number from sim.log to its initial parent path
    path_map = {}
    # keep track of which parent path is in which ensemble
    flow_map = [[] for i in range(n_ensembles)]
    # if true we initialize by reading the first matrix in the sim.log
    # to get the parent path numbers
    initialize = True
    # we found the beginning of the first matrix in sim.log
    found_start_init = False
    # number of paths read from the first matrix in sim.log when initializing
    n_paths_read_init = 0

    with open(log, "r") as rfile:
        for i,line in enumerate(rfile):
            #print(i)
            if "->" in line and not initialize:
                line = line.split()
                # this is a zero swap
                #print(line)
                if len(line)==15:
                    ens0, ens1 = 0, 1
                    pathnr0_old, pathnr1_old = line[-5:-3]
                    pathnr1_new, pathnr0_new = line[-2:]

                    if pathnr1_new != pathnr0_old:
                        path_map[pathnr0_new] = path_map[pathnr0_old]
                        del path_map[pathnr0_old]

                        path_map[pathnr1_new] = path_map[pathnr1_old]
                        del path_map[pathnr1_old]

                    pathnr_ens0 = path_map[pathnr0_new]
                    pathnr_ens1 = path_map[pathnr1_new]
                    flow_map[pathnr_ens0].append(1)
                    flow_map[pathnr_ens1].append(0)
                    for j in range(n_ensembles):
                        if j != pathnr_ens0 and j!= pathnr_ens1:
                            flow_map[j].append(flow_map[j][-1])

                elif len(line)==11:
                    pathnri_old = line[-3]
                    pathnri_new = line[-1]
                    ensi = int(line[-6])

                    if pathnri_new != pathnri_old:
                        path_map[pathnri_new] = path_map[pathnri_old]
                        del path_map[pathnri_old]

                    pathnr_ensi = path_map[pathnri_new]

                    for j in range(n_ensembles):
                        if flow_map[j][-1] == ensi and j != pathnr_ensi:
                            flow_map[j].append(flow_map[pathnr_ensi][-1])
                        elif j != pathnr_ensi:
                            flow_map[j].append(flow_map[j][-1])

                    flow_map[pathnr_ensi].append(ensi)

                else:
                    ValueError("Wrong length in line.split()")
                #print("=="*10)
                #print(line)
                #print(path_map)
                #print(flow_map)
                #print("=="*10)

            elif initialize and "  -- |" in line:
                found_start_init = True

            elif initialize and found_start_init:
                line = line.split()
                pathnri = str(int(line[1][1:]))
                path_map[pathnri] = n_paths_read_init
                flow_map[n_paths_read_init].append(n_paths_read_init)
                n_paths_read_init += 1
                if n_paths_read_init == n_ensembles:
                    initialize = False
                    print(f"[INFO] Initialization done {path_map}")

            else:
                ValueError("Something wrong, but not sure what")

        # write to file
        if out:
            np.savetxt(out, np.array(flow_map).T, fmt="%4d")


        # plotting stuff
        if plot is not None:
            if plot == ['all'] or len(plot) == 0:
                plot_range = range(n_ensembles)
            else:
                plot_range = []
                for i in plot:
                    try:
                        plot_range.append(int(i))
                    except ValueError:
                        print(f"Invalid axis value: {plot}")
            for j in plot_range:
                plt.plot(flow_map[j], label = f"path {j}", marker = "o", markersize = 5)
            plt.axhline(0, c = "k")
            plt.axhline(n_ensembles - 1, c = "k")
            plt.legend()
            plt.show()
