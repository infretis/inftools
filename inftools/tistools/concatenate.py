from typing import Annotated

import typer


def trjcat(
    out: Annotated[
        str, typer.Option("-out", help="output trajectory name, format")
    ],
    traj: Annotated[
        str, typer.Option("-traj", help="path to the traj.txt file")
    ],
    selection: Annotated[
        str,
        typer.Option("-selection", help="atom selection to write to output"),
    ] = "all",
    engine: Annotated[
        str,
        typer.Option("-engine", help="what engine to use for concatenation"),
    ] = "mda",
    traj_format: Annotated[
        str, typer.Option("-format", help="trajectory format, can be inferred")
    ] = None,
    topology: Annotated[
        str, typer.Option("-topology", help="a structure file with atom info")
    ] = None,
    centersel: Annotated[
        str,
        typer.Option(
            "-centersel",
            help="center this atom selction in the middle of the box",
        ),
    ] = "",
):
    """
    Concatenate the trajectories from an infretis simulation to a single file.
    """
    import pathlib

    import numpy as np

    if engine == "mda":
        import MDAnalysis as mda
        from MDAnalysis import transformations as trans
    elif engine == "ase":
        from ase.io.trajectory import Trajectory
        if topology:
            import MDAnalysis as mda

    traj = pathlib.Path(traj)
    out = pathlib.Path(out)
    if not traj.exists():
        raise FileNotFoundError(f"No such file {traj}")
    if out.exists():
        raise FileExistsError(f"File {out} allready exists.")
    if topology:
        topology = pathlib.Path(topology)
        if not topology.exists():
            raise FileNotFoundError(f"No such file {topology}")

    if engine == "ase" and out.suffix != ".traj" and not topology:
        raise ValueError("Supply a topology to write MDA supported formats")

    if engine == "ase" and selection != "all" and not topology:
        raise ValueError("Supply a topology to write specific atom selections")


    traj_dir = traj.parents[0]
    traj_file_arr, index_arr = np.loadtxt(
        str(traj), usecols=[1, 2], comments="#", dtype=str, unpack=True
    )
    index_arr = index_arr.astype(int)

    # trajectory transformations
    if centersel and engine != "mda":
        raise ValueError("centersel is only suppoerted with mda")

    # read the trajectories
    U = {}
    for traj_file in np.unique(traj_file_arr):
        traj_fpath = traj_dir / "accepted" / traj_file
        if not traj_fpath.exists():
            raise FileNotFoundError(f"\n No such file {traj_fpath.resolve()}")
        if engine == "mda":
            if topology is not None:
                u = mda.Universe(
                    topology,
                    str(traj_fpath),
                    format=traj_format,
                    refresh_offsets=True,
                )
            else:
                u = mda.Universe(
                    str(traj_fpath), format=traj_format, refresh_offsets=True
                )
            U[traj_file] = u
            if centersel:
                tmp_sel = u.select_atoms(centersel)
                workflow = [trans.center_in_box(tmp_sel), trans.wrap(u.atoms)]
                u.trajectory.add_transformations(*workflow)
            n_atoms = u.select_atoms(selection).n_atoms
        elif engine == "ase":
            # create mda Universe to write ASE trajectory in MDA suppported formats
            U[traj_file] = Trajectory(str(traj_fpath))
            if topology is not None:
                mda_universe = mda.Universe(topology)
                n_atoms = mda_universe.select_atoms(selection)
            else:
                n_atoms = U[traj_file][0].positions.shape[0]

    # write the concatenated trajectory
    if engine == "mda":
        with mda.Writer(str(out), n_atoms) as wfile:
            import tqdm
            for traj_file, index in tqdm.tqdm(zip(traj_file_arr, index_arr)):
                u = U[traj_file]
                ag = u.select_atoms(selection)
                u.trajectory[index]
                wfile.write(ag.atoms)

    elif engine == "ase" and not topology:
        out = Trajectory(str(out), "w")
        for traj_file, index in zip(traj_file_arr, index_arr):
            # traj object
            u = U[traj_file]
            atoms = u[index]
            out.write(atoms)
    else:
        with mda.Writer(str(out), n_atoms) as wfile:
            for traj_file, index in zip(traj_file_arr, index_arr):
                u = U[traj_file]
                atoms = u[index]
                mda_universe.atoms.positions = atoms.positions
                ag = mda_universe.select_atoms(selection)
                wfile.write(ag.atoms)

    print("\nAll done!")
