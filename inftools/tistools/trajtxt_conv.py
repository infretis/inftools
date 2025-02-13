import argparse
from typing import Annotated
import typer

def trajtxt_conv(
    i: Annotated[str, typer.Option("-i")],
    o: Annotated[str, typer.Option("-o")],
    r: Annotated[int, typer.Option("-r", help="Subcycle ratio, subcycle xtc /subsycle trr")],
    ):
    """inft trjcat cannot be directly used when xtc is kept
    while delete_old=true during a gromacs infretis sim.

    Here we provide a quick traj.txt converter, that also consders
    if xtc frames are saved more frequently."""

    import numpy as np

    # Assert and data
    assert r >= 1
    data = np.loadtxt(i, dtype="str")
    data_out = []
    files = []

    # First get unique path files
    for file in data[:, 1]:
        if file not in files:
            files.append(file)

    # iterate over files
    for file in files:

        # get first, last and direction of frames in file
        file_idx_s = int(data[data[:, 1] == file][0][2])
        file_idx_e = int(data[data[:, 1] == file][-1][2])
        it = int(data[data[:, 1] == file][0][-1])

        # fill new data list
        for idx in np.arange(file_idx_s*r, file_idx_e*r + it, it):
            data_out.append([file.replace(".trr", ".xtc"), idx, it])

    with open(o, "w") as write:
        write.write("# comment\n")
        write.write("# comment\n")
        for idx, line in enumerate(data_out):
            write.write(f"{idx:10.0f}    {line[0]}{line[1]:12.0f}{line[2]:7.0f}\n")

    print(f"trajtxt_conv Done! with file {o}")
