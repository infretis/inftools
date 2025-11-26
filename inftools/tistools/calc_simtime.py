import typer
from typing import Annotated, Optional


def calc_simtime(
    log: Annotated[str, typer.Option("-log")] = "sim.log",
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
    pstarts = []
    tstarts = []
    # previous, current time
    ptime, ctime = None, None

    with open(log, "r") as read:
        for line in read:
            if "submit worker 0 START" in line:
                ptime = datetime.strptime(line[-20:-1], format_str)
                ctime = None
                pstarts.append(len(paths))
                tstarts.append(np.sum(paths)/3600/24)
            if "shooting" in line:
                line = read.readline()
                if "date" in line:
                    if ctime is not None:
                        ptime = ctime
                    rip = " ".join(line.rstrip().split()[2:])
                    ctime = datetime.strptime(rip, format_str)
                    paths.append((ctime-ptime).total_seconds())

    # tstarts = [np.sum(paths[:i])/3600/24 for i in pstarts]

    if plot:
        # plt.plot(np.arange(len(paths)), np.cumsum(paths)/3600/24)
        plt.plot(np.cumsum(paths)/3600/24, np.arange(len(paths)))
        np.savetxt("simtime.txt", np.array([np.cumsum(paths)/3600/24, np.arange(len(paths))]).T)
        for pstart, tstart in zip(pstarts, tstarts):
            # plt.axvline(np.sum(paths[:start]))
            plt.axvline(tstart, color="k", ls="--")
            # plt.axhline(pstart, color="k", ls="--")
        plt.ylabel("Shooting Attempts")
        plt.xlabel("Time [Days]")
        plt.show()

    print(f"Total Wall Time: {np.sum(paths)/3600/24:.01f} Days")
    print(f"Total Restarts: {len(tstarts)-1}")
    print(f"Total Sampled Paths: {len(paths)}")

    return np.sum(paths)/3600/24, len(tstarts)-1
