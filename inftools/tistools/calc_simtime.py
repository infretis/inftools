import typer
from typing import Annotated, Optional


def calc_simtime(
    log: Annotated[str, typer.Option("-log")] = "sim.log",
    skip: Annotated[int, typer.Option("-skip")] = 100,
    plot: Annotated[bool, typer.Option("-plot")] = True,
    ):
    """Calculate the total simulation wall time while
    considering restarts. Basically by calculating the delta time."""

    import matplotlib.pyplot as plt
    import numpy as np
    import time

    from datetime import datetime
    format_str = "%Y.%m.%d %H:%M:%S"

    paths = []
    starts = []
    # previous, current time
    ptime, ctime = None, None

    with open(log, "r") as read:
        for line in read:
            if "submit worker 0 START" in line:
                ptime = datetime.strptime(line[-20:-1], format_str)
                ctime = None
                starts.append(len(paths))
            if "shooting" in line:
                line = read.readline()
                if "date" in line:
                    if ctime is not None:
                        ptime = ctime
                    rip = " ".join(line.rstrip().split()[2:])
                    ctime = datetime.strptime(rip, format_str)
                    paths.append((ctime-ptime).total_seconds())

    # enforce skip
    paths = paths[skip:]

    if plot:
        plt.plot(np.arange(len(paths)), np.cumsum(paths)/3600/24)
        for start in starts:
            # plt.axvline(np.sum(paths[:start]))
            plt.axvline(start, color="k", ls="--")
        plt.xlabel("Sampled Paths")
        plt.ylabel("Time [Days]")
        plt.show()

    return np.sum(paths)/3600/24, len(starts)-1
