from typing import Annotated, Tuple
import typer

# Disable automatic underscore -> hypen in CLI names
typer.main.get_command_name = lambda name: name

import argparse
import glob
import os
import shutil
import subprocess
from types import SimpleNamespace

import numpy as np
import tomli
from infretis.classes.orderparameter import Puckering


def check_indices(
        sdf: Annotated[str, typer.Option("-sdf", help="The .sdf file of your molecule (e.g. mol.sdf)")],
        idx: Annotated[Tuple[int, int, int, int, int, int], typer.Option("-idx", help="The ordered indices of your molecule (e.g. 2 5 11 8 1 0)")],
    ):
    """Calculate the theta and phi angle for an .sdf file given a set of indices."""

    orderparameter = Puckering(index=idx)

    subprocess.run(f"obabel -isdf {sdf} -oxyz -O .temp.xyz", shell=True)
    x = np.loadtxt(".temp.xyz", skiprows=2, usecols=[1, 2, 3])
    subprocess.run("rm .temp.xyz", shell=True)

    system = SimpleNamespace(
        pos=x,
        box=[np.inf, np.inf, np.inf],
    )
    op = orderparameter.calculate(system)

    print(f"\nTheta = {op[0]:.3f} degrees")
    print(f"Phi = {op[1]:.3f} degrees")

def concatenate(
        out: Annotated[str, typer.Option("-out", help="the outfile trajectory name (e.g. md-traj.xyz)")],
        tpr: Annotated[str, typer.Option("-tpr", help="the .tpr file (e.g. ../gromacs_input/topol.tpr)")],
        path: Annotated[str, typer.Option("-path", help="the filepath to the path_nr folder (e.g. run3/46/)")],
    ):
    "Reverse and concatenate trajectories from RETIS simulations."
    import MDAnalysis as mda
    from MDAnalysis.analysis.align import alignto
    from MDAnalysis.lib.mdamath import make_whole

    # parser.add_argument(
    #     "--selection",
    #     help="The selection, e.g. 'index 1 2 3 4' or 'resname MOL0' \
    #        (default 'resname MOL0')",
    #     default="resname MOL0",
    # )

    selection = "resname MOL0"


    traj_file_arr, index_arr = np.loadtxt(
        f"{path}/traj.txt",
        usecols=[1, 2],
        comments="#",
        dtype=str,
        unpack=True,
    )
    traj_file_arr = [f"{path}/accepted/{traj_i}" \
            for traj_i in traj_file_arr]
    # traj_file_arr = np.char.replace(traj_file_arr,"trr","xtc")
    index_arr = index_arr.astype(int)

    U = {}
    # reference for rmsd alignment
    for traj_filename in traj_file_arr:
        if traj_filename.endswith(".g96"):
            traj_file_arr.remove(traj_filename)
    ref = mda.Universe(tpr, traj_file_arr[0]).select_atoms(selection)
    make_whole(ref.atoms)
    for traj_file in np.unique(traj_file_arr):
        print(f"Reading {traj_file} ...")
        subprocess.run(
            f"printf '1\n0\n' | gmx trjconv \
            -f {traj_file} -o tmp.{traj_file.split('/')[-1]} \
            -pbc whole -center -s {tpr}",
            shell=True,
        )
        if not os.path.exists(traj_file):
            exit(f"Could not find file {traj_file}.?")

        U[traj_file] = mda.Universe(
                tpr, f"tmp.{traj_file.split('/')[-1]}"
                )

    with mda.Writer(
        out, U[traj_file].select_atoms(selection).n_atoms
    ) as wfile:
        for traj_file, index in zip(traj_file_arr, index_arr):
            u = U[traj_file]
            ag = u.select_atoms(selection)
            make_whole(ag)
            u.trajectory[index]
            alignto(ag, ref)
            wfile.write(ag.atoms)

    for traj_file in np.unique(traj_file_arr):
        subprocess.run(f"rm tmp.{traj_file.split('/')[-1]}", shell=True)

    subprocess.run(f"perl -pi -e 'chomp if eof' {out}", shell=True)

    print("\nAll done!")
    print(f"Trajectory written to {out}.")

def generate_openff_topology(
        sdf: Annotated[str, typer.Option("-sdf", help="The .sdf file of your molecule (e.g. mol.sdf)")],
    ):
    """Generate an OpenFF topology for the puckering example."""
    from openff.interchange import Interchange
    from openff.toolkit import ForceField, Molecule
    from openff.units import unit
    # load the molecule
    mol = Molecule.from_file(sdf)
    topology = mol.to_topology()
    topology.box_vectors = unit.Quantity([2.1, 2.1, 2.1], unit.nanometer)
    # Load the OpenFF 2.1.0 forcefield called "Sage"
    sage = ForceField("openff-2.1.0.offxml")
    out = Interchange.from_smirnoff(force_field=sage, topology=topology)
    out.to_gro("gromacs_input/mol.gro")
    out.to_top("gromacs_input/mol.itp")

    with open("gromacs_input/topol.top", "w") as writefile:
        with open("gromacs_input/mol.itp") as readfile:
            for line in readfile:
                if "[ moleculetype ]" in line:
                    writefile.write("; Include tip3p water topology\n")
                    writefile.write('#include "amber99.ff/ffnonbonded.itp"\n')
                    writefile.write('#include "amber99.ff/tip3p.itp"\n\n')
                writefile.write(line)

def initial_path_from_iretis(
    traj: Annotated[str, typer.Option("-traj", help="The path to the folder containing the trajectories (e.g. ../iretis0/trajs/)")],
    toml: Annotated[str, typer.Option("-toml", help="The .toml input file for reading the interfaces (e.g. 'infretis.toml')")],
        ):
    "Generate initial paths for an infretis simulation using paths from an earlier infretis simulation"

    # read interfaces from .toml file
    with open(toml, "rb") as toml_file:
        toml_dict = tomli.load(toml_file)
    interfaces = toml_dict["simulation"]["interfaces"][:-1]

    out = {}  # ensemble - traj_idx

    trajs = glob.glob(f"{traj}/*")  # folder to trajectories
    np.argsort([int(f.split("/")[-1]) for f in trajs])
    trajs = sorted(trajs, key=os.path.getctime)

    # iterate backwards to get decorrelated paths
    for traj in trajs[::-1]:
        x = np.loadtxt(f"{traj}/order.txt", usecols=[0, 1])
        # zero minus
        if x[1, 1] < interfaces[0]:
            if 0 not in out.keys():
                out[0] = traj

        # 0+ intf
        else:
            omax = np.max(x[:, 1])
            valid_in = False
            for i, interface in enumerate(interfaces):
                if omax > interface:
                    valid_in = i + 1
            if valid_in and valid_in not in out.keys():
                out[valid_in] = traj

    # if we miss some lower ensembles, add to
    # them the paths from the higher ensembles
    for i in range(len(interfaces) + 1):
        if i not in out.keys():
            for j in range(i + 1, len(interfaces) + 1):
                if j in out.keys():
                    out[i] = out[j]
                    print(f"[INFO] Added path from ens{j} to ens{i}")


    # Check if we have paths in all ensembles
    for i in range(len(interfaces) + 1):
        assert (
            i in out.keys()
        ), f"Did not find any paths in ensemble {i}\
    that cross the corresponding interface"

    loaddir = "load"
    if os.path.exists(loaddir):
        exit(
            f"\nDirectory {loaddir}/ exists. Will not overwrite.\
    \nRename or delete it manually. Aborting."
        )
    else:
        os.mkdir(loaddir)

    for i, traj in zip(out.keys(), out.values()):
        shutil.copytree(traj, f"{loaddir}/{i}")

    print("\nAll done! Created folder load/ with new initial paths.")

def initial_path_from_md(
        trr: Annotated[str, typer.Option("-trr", help="The .trr trajectory file")],
        toml: Annotated[str, typer.Option("-toml", help="The .toml input for reading the interfaces")],
        order: Annotated[str, typer.Option("-order", help="The order file corresponding to the trajectory")],
        ):
    "Generate initial paths for an infretis simulation from an equilibrium run."

    predir = "load"
    if os.path.exists(predir):
        exit(
            f"\nDirectory {predir}/ exists."
            + " Will not overwrite."
            + " Rename or remove it and try again."
        )
    else:
        os.mkdir(predir)

    traj = trr  # trajectory  file
    order = np.loadtxt(order)  # order file

    # read interfaces from .toml file
    with open(toml, "rb") as toml_file:
        toml_dict = tomli.load(toml_file)
    interfaces = toml_dict["simulation"]["interfaces"]
    if len(interfaces) == 0:
        print(f"\nNo interfaces defined in '{toml}'!")
        return

    for i in range(len(interfaces)):
        dirname = os.path.join(predir, str(i))
        accepted = os.path.join(dirname, "accepted")
        trajfile = os.path.join(accepted, "traj.trr")
        orderfile = os.path.join(dirname, "order.txt")
        trajtxtfile = os.path.join(dirname, "traj.txt")
        print(f"Making folder: {dirname}")
        os.makedirs(dirname)
        print(f"Making folder: {accepted}")
        os.makedirs(accepted)
        print(
            "Writing trajectory {} and order {} and trajfile {}".format(
                trajfile, orderfile, trajtxtfile
            )
        )

        # minus ensemble
        if i == 0:
            idx = (order[:, 1] > interfaces[0]).astype(int)
            grad = idx[1:] - idx[:-1]
            # hopping above interface0 grad = 1
            above = np.where(grad == 1)[0]
            # select an ending point where the path hops
            # above interface0
            end = above[-2] + 1  # we want the ending point
            # select a starting point where the path hops
            # below interface0, and has to precede the end point
            # only look at the array up til end point
            below = np.where(grad[:end] == -1)[0]
            start = below[-1]
            iterator = [i for i in range(start, end + 1)]
            print("=" * 10)
            print(iterator)
            print(order[iterator, 1])

        # plus ensembles
        else:
            idx = (order[:, 1] > interfaces[0]).astype(int)
            grad = idx[1:] - idx[:-1]
            # hopping above interface0 grad = 1
            above = np.where(grad == 1)[0]
            # select a starting point where the path hops
            # above interface0. Dont select last point
            # as we may not jump below again after that
            start = above[-2]
            # select an ending point where the path hops
            # below interface0
            # truncate where wee look for this point
            below = np.where(grad[: above[-1]] == -1)[0]
            end = below[-1] + 1  # only look at the array up til end point
            iterator = [i for i in range(start, end + 1)]
            print("=" * 10)
            print(iterator)
            print(order[iterator, 1])

            # check if valid path for wire-fencing
            idx = np.where(
                (order[iterator, 1] >= interfaces[i - 1])
                & (order[iterator, 1] <= interfaces[i])
            )[0]
            if len(idx) == 0 and i > 1:
                # no points between interface i-1 and i
                idx = np.where(
                    (order[iterator, 1] >= interfaces[i - 1])
                    & (order[iterator, 1] <= interfaces[i + 1])
                )[0]
                exit("Invalid path for wf!!")

        with open(".frames.ndx", "w") as index_file:
            index_file.write("[ frames ]\n")
            for idxi in iterator:
                index_file.write(f"{idxi+1}\n")

        cmd = f"gmx trjconv -f {traj} -o {trajfile} -fr .frames.ndx"
        print(cmd)
        subprocess.run(cmd, shell=True)

        # write order file
        N = len(iterator)
        np.savetxt(
            orderfile,
            np.c_[order[:N, 0], order[iterator, 1:]],
            header=f"{'time':>10} {'order':>15}",
            fmt=["%10.d"] + ["%15.4f" for fmti in order[0, 1:]],
        )
        np.savetxt(
            trajtxtfile,
            np.c_[
                [str(i) for i in range(N)],
                ["traj.trr" for i in range(N)],
                [str(i) for i in range(N)],
                [str(1) for i in range(N)],
            ],
            header=f"{'time':>10} {'trajfile':>15} {'index':>10} {'vel':>5}",
            fmt=["%10s", "%15s", "%10s", "%5s"],
        )

    print("\nAll done! Created folder load/ containing the initial paths.")

def plot_order(
        traj: Annotated[str, typer.Option("-traj", help="The path to the folder containing the trajectories (e.g. 'run1/load/')")],
        toml: Annotated[str, typer.Option("-toml", help="The .toml input file for reading the interfaces (e.g. 'infretis.toml')")],
        xy: Annotated[Tuple[int, int], typer.Option("-xy", help="The indices of the columns to plot (default 0 1)")] = (0, 1),
        skip: Annotated[bool, typer.Option("-skip" , help="skip initial load paths")] = False,
    ):
    "Plot the order parameter of all paths from an infretis simulation."
    import matplotlib.pyplot as plt
    # read interfaces from the .toml file
    with open(toml, "rb") as toml_file:
        toml_dict = tomli.load(toml_file)
    interfaces = toml_dict["simulation"]["interfaces"]

    # get the filenames of all created paths
    paths = glob.glob(f"{traj}/*/order.txt")

    # sort filenames by time of path creation
    sorted_paths = sorted(paths, key=os.path.getctime)

    # plotting stuff, modify by your needs
    f, a = plt.subplots()

    # add horisontal lines for interfaces
    for interface in interfaces:
        a.axhline(interface, c="k", lw=0.5)

    if 2 in xy:
        lw = 0

    else:
        lw = 1

    # plot all paths, modify by your needs
    if skip:
        sorted_paths = sorted_paths[len(interfaces):]
    for path in sorted_paths:
        x = np.loadtxt(path)
        if x.shape[1] > 2:
           if x[-1, 1] > interfaces[-1]:
                print(
                    f"The path in {path} is reactive with \
        phi={x[-1,2]:.2f}! \U0001F389 \U0001F938 \U0001F483"
                )
           #    continue # continues to next iteration in loop
        a.plot(
            x[:, xy[0]],
            x[:, xy[1]],
            c="C0",
            marker="o",
            markersize=2.5,
            lw=lw,
        )

    plt.show()
