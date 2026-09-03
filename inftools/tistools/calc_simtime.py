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

    mcmoves = []
    pstarts = []
    tstarts = []
    # previous, current time
    ptime, ctime = None, None

    with open(log, "r") as read:
        for line in read:
            if "submit worker 0 START" in line:
                ptime = datetime.strptime(line[-20:-1], format_str)
                ctime = None
                pstarts.append(len(mcmoves))
                tstarts.append(np.sum(mcmoves)/3600/24)
            if "[INFO]: date:" in line:
                rip = " ".join(line.rstrip().split()[2:])
                ctime = datetime.strptime(rip, format_str)
                end = read.readline()
                if "END" not in end:
                    continue
                if ptime is not None:
                    delta = (ctime - ptime).total_seconds()
                    mcmoves.append(delta)
                ptime = ctime

    if plot:
        plt.plot(np.cumsum(mcmoves)/3600/24, np.arange(len(mcmoves)))
        np.savetxt("simtime.txt", np.array([np.cumsum(mcmoves)/3600/24, np.arange(len(mcmoves))]).T)
        for pstart, tstart in zip(pstarts, tstarts):
            # plt.axvline(np.sum(paths[:start]))
            plt.axvline(tstart, color="k", ls="--")
            # plt.axhline(pstart, color="k", ls="--")
        plt.ylabel("Shooting Attempts")
        plt.xlabel("Time [Days]")
        plt.show()

    print(f"Total Wall Time: {np.sum(mcmoves)/3600/24:.02f} Days")
    print(f"Total Wall Time: {np.sum(mcmoves)/3600:.02f} Hours")
    print(f"Total Restarts: {len(tstarts)-1}")
    print(f"Total MC Moves: {len(mcmoves)}")

    return np.sum(mcmoves)/3600/24, len(tstarts)-1
