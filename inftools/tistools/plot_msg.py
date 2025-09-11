from typing import Annotated

import typer


def plot_msg(
    toml: Annotated[str, typer.Option("-toml")] = "infretis.toml",
    add_markers: Annotated[bool, typer.Option("--add_markers",
        help = "Also plot markers on discrete phasepoints")] = False,
):
    """Plot the order printed in the worker*/msg* files, so we
    can visualize the progress of the shooting moves.

    Require that the user be in the root sim folder.

        TODO: color code last accepted WF subpath and the extension
        phase."""
    import os

    import matplotlib.pyplot as plt
    import numpy as np

    from inftools.misc.infinit_helper import read_toml

    intfs = []
    # if toml, plot interfaces
    if os.path.isfile(toml):
        config = read_toml(toml)
        intfs = config["simulation"]["interfaces"]
        cap = config["simulation"]["tis_set"].get("interface_cap", None)

    # iterate through worker folders
    workers = 0
    for i in range(9999):
        wfolder = f"worker{i}"
        if os.path.isdir(wfolder):
            workers += 1
        else:
            break

    # Set grid size
    x_grid, y_grid = 1, 1
    while x_grid * y_grid < workers:
        x_grid += 1
        if x_grid * y_grid >= workers:
            break
        y_grid += 1
        if x_grid * y_grid >= workers:
            break

    fig, axes = plt.subplots(x_grid, y_grid)
    # Try flatten if more than 1 worker
    try:
        axes = axes.flatten()
    except:
        axes = [axes]

    # iterate through worker folder
    for i in range(workers):
        wfolder = f"worker{i}"
        if not os.path.isdir(wfolder):
            break
        idx = 0
        ens = -1

        # get op from msg files in sorted order
        files = [j for j in os.listdir(wfolder) if "msg" in j]
        files_idx = [int(j.split("_")[-2]) for j in files]
        files = np.array(files)[np.argsort(files_idx)]

        for file in files:
            ens = int(file[4:7])
            try:
                msg = np.loadtxt(wfolder + "/" + file)
                if "trajB" in file:
                    axes[i].plot(
                        msg[:, 0] + idx, msg[::-1, 1], ls="--", color=f"C{i%8}",
                        marker="^" if add_markers else None,
                    )
                    axes[i].scatter(
                        [msg[0, 0] + idx],
                        [msg[-1, 1]],
                        color=f"C{i%8}",
                        marker="<",
                        edgecolor="k",
                        zorder=100,
                    )
                    axes[i].scatter(
                        [msg[-1, 0] + idx],
                        [msg[0, 1]],
                        color="k",
                        marker="x",
                        zorder=100,
                    )
                else:
                    axes[i].plot(
                        msg[:, 0] + idx - 1, msg[:, 1], color=f"C{i%8}",
                        marker="o" if add_markers else None,
                    )
                    axes[i].scatter(
                        [msg[-1, 0] + idx - 1],
                        [msg[-1, 1]],
                        color=f"C{i%8}",
                        marker=">",
                        edgecolor="k",
                    )
                idx += len(msg[:, 0])
            except:
                print(
                    "cannot plot",
                    wfolder + "/" + file,
                    "perhaps no op printed yet",
                )

        # plt ensemble intfs
        if intfs and ens > -1:
            intf_idx = ens if ens == 0 else ens - 1

            axes[i].axhline(intfs[intf_idx], color=f"C{i%8}", alpha=0.2)
            axes[i].set_title(f"{ens:03}-worker{i}")
            axes[i].axhline(intfs[0], color="k")
            axes[i].axhline(intfs[-1], color="k")
            if cap is not None:
                axes[i].axhline(cap, color="r", ls="--")

    if workers > 0:
        plt.show()
    else:
        print("Not in the root infretis simulation folder!")
