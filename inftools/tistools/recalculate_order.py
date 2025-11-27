from typing import Annotated, Tuple
from enum import Enum
import typer

# Disable automatic underscore -> hyphen in CLI names
# typer.main.get_command_name = lambda name: name

# for Choices CLI parameter
class Format(str, Enum):
    ONE = "trr"
    TWO = "xyz"
    THR = "traj"
    FOU = "rkf"



def recalculate_order(
    toml: Annotated[str, typer.Option("-toml")] = "infretis.toml",
    traj: Annotated[str, typer.Option("-traj")] = "traj.xyz",
    out: Annotated[str, typer.Option("-out", help="the output of the analysis")] = "order_rec.txt",
    format: Annotated[str, typer.Option("-format", case_sensitive=False, help="the file format of the trajectory (.traj is ase format)")]= "trr",
    box: Annotated[Tuple[float, float, float], typer.Option("-box", help="xyz only; box dimensions in angstrom (e.g. 30 30 30)")] = None,
    ):
    """
    Recalculate the orderparamter from a trajectory file.

    TODO:

    * Read velocity information, and set vel_rev in config
    """

    import os
    import numpy as np
    import tomli
    import struct

    import MDAnalysis as mda

    from infretis.classes.engines.cp2k import read_xyz_file
    from infretis.classes.engines.gromacs import read_gromos96_file
    from infretis.classes.orderparameter import create_orderparameter
    from infretis.classes.system import System
    from inftools.analysis.gromacs import read_trr_file

    if os.path.exists(out):
        raise FileExistsError(f"File {out} exists, aborting")

    with open(toml, "rb") as toml_file:
        toml_dict = tomli.load(toml_file)
    orderparameter = create_orderparameter(toml_dict)
    # interfaces = toml_dict["simulation"]["interfaces"]

    if not traj.endswith(format):
        msg = f"[ WARNING ] {traj} may not be a {format} formatted"\
        " file. Use the -format flag to be sure you get some output"
        print(msg)
    if format == "xyz":
        if box is None:
            raise ValueError("Provide a box size for xyz files.")
        traj = read_xyz_file(traj)
        box = np.array(box, dtype=float).flatten()
    elif format == "traj":
        from ase.io.trajectory import Trajectory
        traj = Trajectory(traj)
    elif format == "trr":
        traj = read_trr_file(traj)
    elif format == "g96":
        _, xyz, vel, box = read_gromos96_file(traj)
        traj = [[xyz, vel, box]]
    elif format == "rkf":
        from scm.plams.trajectories.rkffile import RKFTrajectoryFile
        traj = RKFTrajectoryFile(traj)
        traj.store_mddata()
        mol = traj.get_plamsmol()
    else:
        u = mda.Universe(traj, format = format)
        traj = u.trajectory


    with open(out, "w") as writefile:
        writefile.write("# step\torder\n")
        for i, frame in enumerate(traj):
            if format == "xyz":
                pos=np.vstack((frame["x"], frame["y"], frame["z"])).T
            elif format == "traj":
                pos = frame.positions
                box = np.diag(frame.cell.array)
            elif format == "trr":
                pos=frame[1]["x"]
                box=np.diag(frame[1]["box"])
            elif format == "g96":
                pos = frame[0]
                box = frame[2]
            elif format == "rkf":
                pos, box = traj.read_frame(i, molecule=mol)
                box = np.diag(box)
            else:
                pos = u.atoms.positions
                box = u.dimensions[:3]

            system = System()
            system.config = (traj, i)
            system.pos = pos
            system.box = box

            op = orderparameter.calculate(system)
            line = f"{i} " + " ".join([f"{opi}" for opi in op]) + "\n"
            writefile.write(line)

    print(f"[ INFO ] Orderparameter values written to {out}\n")

    return
