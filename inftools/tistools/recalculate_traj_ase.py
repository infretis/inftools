from typing import Annotated
import typer

def recalculate_traj_ase(
    path: Annotated[str, typer.Option("-path", help="Path to the path's path number, usually in load/{pn}")],
    toml: Annotated[str, typer.Option("-toml")] = "infretis.toml",
    out: Annotated[str, typer.Option("-out", help="The output of the analysis")] = "order_rec.txt",
):
    """
    Recalculate the order parameter from a sampled infretis path
    using ASE .traj files instead of MDAnalysis-readable trajectories.
    """

    import os
    import numpy as np
    import tomli
    from ase.io.trajectory import Trajectory
    from infretis.classes.orderparameter import create_orderparameter
    from infretis.classes.system import System

    # --- Load configuration ---
    with open(toml, "rb") as toml_file:
        toml_dict = tomli.load(toml_file)
    orderparameter = create_orderparameter(toml_dict)

    # --- Load trajectory index ---
    tdata = np.loadtxt(os.path.join(path, "traj.txt"), dtype=str)
    files = set(tdata[:, 1])

    # --- Load ASE trajectories ---
    unis = {}
    for f in files:
        traj_path = os.path.join(path, "accepted", f)
        if not os.path.isfile(traj_path):
            raise FileNotFoundError(f"Trajectory file not found: {traj_path}")
        unis[f] = Trajectory(traj_path)

    # --- Handle CP2K boxes (for xyz trajectories) ---
    if "xyz" in [i.split('.')[-1] for i in files]:
        from infretis.classes.engines.cp2k import read_cp2k_box
        cp2k_inp = os.path.join(toml_dict["engine"]["input_path"], "cp2k.inp")
        box, pbc = read_cp2k_box(cp2k_inp)
        for traj in unis.values():
            for atoms in traj:
                atoms.set_cell(box)
                atoms.set_pbc(True)

    # --- Prepare output file ---
    out_file = os.path.join(path, out)
    for i in range(1, 1000):
        if not os.path.isfile(out_file):
            break
        out_file = os.path.join(path, f"order_rec_{i}.txt")

    # --- Recalculate order parameters ---
    with open(out_file, "w") as writefile:
        writefile.write("# step\torder\n")
        for step, fname, index, vel in tdata:
            traj = unis[fname]
            idx = int(index)
            atoms = traj[idx]

            system = System()
            system.config = (fname, index)
            system.pos = atoms.get_positions()
            system.box = atoms.get_cell().lengths()

            op = orderparameter.calculate(system)
            line = f"{step} " + "\t".join(map(str, op)) + "\n"
            writefile.write(line)

    print(f"[ INFO ] Orderparameter values written to {out_file}\n")

