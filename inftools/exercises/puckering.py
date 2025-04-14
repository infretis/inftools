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
import pathlib

import numpy as np
import tomli
import tomli_w
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
    traj: Annotated[str, typer.Option("-traj", help="The path to the trajectory folder (e.g. run0/)")],
    toml: Annotated[str, typer.Option("-toml", help="The .toml input file for reading interfaces")],
    out_dir: Annotated[str, typer.Option("-out_dir", help="The out directory containing the new paths")] = "load",
    restart: Annotated[str, typer.Option("-restart", help=".toml to restart file for reading last active paths")] = "",
    out_toml: Annotated[str, typer.Option("-out_toml", help="Output file if interfaces and shooting_moves are chagned")] = "",
    keep_all_active: Annotated[bool, typer.Option(help = "If active paths are no longer valid, add new interfaces")] = False,
    active_path_dir: Annotated[str, typer.Option(help = "Directory to the active paths ('-traj' if not given)")] = "",
    ):
    """Pick out initial paths from an earlier infretis simulation.

    It is assumed that the paths present in the '-traj' folder where generated
    with identical interfaces[0] as specified in the '-toml' file.
    We assume shooting_moves = ['sh', 'sh', 'wf', ...,  'wf']

    If a restart.toml file is supplied with '-restart', the last active paths
    from the previous simulation can be used. If some of these are not valid
    with the new interfaces, additional interfaces can be added to accommodate
    these paths by giving invoking the '--keep_all_active' flag. The active
    paths are valuable, since they are the most decorrelated paths from the
    initial paths and discarding them is wasteful. If '--keep_all_active' paths
    is set, new interfaces may be added from the '-restart' .toml file in which
    the active paths where valid. It will also add as many more 'wf' moves to
    the shooting_moves as required.

    TODO:
        * Ideally, to save space, we should simply move them instead of copying.
        However, that requires us to use symlinks because if we move one path that
        is going to be used again (cloned), then we have to point to that path.
        * Check that selected paths are actually valid wf paths using infretis
        functions
    """
    out_dir = pathlib.Path(out_dir)
    toml = pathlib.Path(toml)
    if not active_path_dir:
        # just use the most recent changed run dir
        active_path_dir = sorted(glob.glob(f"{traj}"), key = os.path.getctime)[-1]
    active_path_dir = pathlib.Path(active_path_dir)
    if not active_path_dir.exists():
        raise ValueError(f"No such directory: {active_path_dir}")
    if out_toml:
        out_toml = pathlib.Path(out_toml)

    if out_dir.exists():
        raise ValueError(
            f"Directory {out_dir.resolve()} exists. Will not overwrite. "
            "Rename or delete it manually. Aborting."
        )
    else:
        os.mkdir(out_dir)

    # read interfaces from .toml file
    with open(toml, "rb") as toml_file:
        toml_dict = tomli.load(toml_file)

    if restart:
        with open(restart, "rb") as toml_file:
            restart_dict = tomli.load(toml_file)

    interfaces = toml_dict["simulation"]["interfaces"]
    sh_m = toml_dict["simulation"]["shooting_moves"]
    trajs = [pathlib.Path(i) for i in glob.glob(f"{traj}/*")]
     # sort so most recent trajectory files are considered first
    trajs = sorted(trajs, key=os.path.getctime)[::-1]

    # try to read active paths
    if restart and restart_dict.get("current", False):
        active_paths = restart_dict["current"].get("active", False)
        if not active_paths:
            print(f"* {restart} has 'current' section but no 'active' section.")
        else:
            active_paths = [active_path_dir / str(ap) for ap in active_paths]
            print("* Active paths found. Trying to use previous active paths.")

    out = {}  # ensemble:traj_idx in load folder

    if restart and active_paths:
        # pick active [0-] and [0+] paths, since interfaces[0] doesn't change
        out[0] = active_paths[0]
        out[1] = active_paths[1]

        # first load order.txt to calculate maxop
        omax_active_paths = []
        for traj in active_paths:
            order_file = traj / "order.txt"
            # ensure we are reading an actual path directory
            if not order_file.exists():
                print(f"* Did not find {order_file.resolve()}")
                continue
            x = np.loadtxt(order_file, usecols=[0, 1])
            omax = np.max(x[:, 1])
            omax_active_paths.append(omax)

        # sort wrt to increasing maxop, but exclude [0-] and [0+] path
        sorted_a_idx = np.argsort(omax_active_paths[2:])
        # sorted active (sa) paths, but [0-] and [0+] stay in 0 and 1
        sa_omax = omax_active_paths[:2] + [omax_active_paths[2:][i] for i in sorted_a_idx]
        sa_paths = active_paths[:2] + [active_paths[2:][i] for i in sorted_a_idx]
        # put the paths (sorted wrt max op) in increasing ensembles, add
        # interfaces from restart if needed
        for traj, omax in zip(sa_paths[2:], sa_omax[2:]):
            is_valid_path = False
            valid_in = 0
            for i, interface in enumerate(interfaces[:-1]):
                if omax > interface:
                    # valid_in = 2 correspsonds to ensemble [1+]
                    valid_in = i + 1
                    if valid_in not in out.keys():
                        # stop on first valid path for each ensemble
                        is_valid_path = True
                        out[valid_in] = traj
                        break
                else:
                    break
            if not is_valid_path:
                print(f"* Previous active path {traj} is not valid with the"
                f" new interfaces (omax {omax} valid in {valid_in} intf {interface})")
                if keep_all_active:
                    # if a path is not valid, e.g. in [i+], (i > 0) add a new
                    #interface halfway between interfaces[i-1] and omax
                    new_intf = np.round(interfaces[valid_in - 1] + (omax - interfaces[valid_in-1])/2, 4)
                    print(f"* keep_all_active flag set. Adding a new interface"
                    f" at {new_intf}")
                    out[valid_in+1] = traj
                    interfaces = interfaces[:valid_in] + [new_intf] + interfaces[valid_in:]
                    if sh_m[valid_in] == "sh":
                        # sh sh wf should not be sh wf sh wf
                        sh_m = sh_m + ["wf"]
                    else:
                        sh_m = sh_m[:valid_in] + ["wf"] + sh_m[valid_in:]

    if out_toml and keep_all_active:
        if out_toml.resolve() == toml.resolve():
            print(f"* {toml} and {out_toml} are the same file.")
            new_toml_name = toml.resolve().parent / "old_infretis.toml"
            toml.rename(new_toml_name)
            print(f"* Renamed {toml} to {new_toml_name}")

        toml_dict["simulation"]["interfaces"] = interfaces
        toml_dict["simulation"]["shooting_moves"] = sh_m
        with open(out_toml, "wb") as wfile:
            tomli_w.dump(toml_dict, wfile)
            print(f"* Modified interfaces of {out_toml}")

    # Now clone paths to fill up the remaining ensembles
    np.random.shuffle(trajs) # We pick paths by random
    max_valid_in = -np.inf
    global_maxop = -np.inf
    for traj in trajs:
        order_file = traj / "order.txt"
        if not order_file.exists():
            print(f"* Did not find {order_file.resolve()}")
            continue
        x = np.loadtxt(order_file, usecols=[0, 1])
        # zero minus path (can include lambda_minus_one)
        if x[1,1] < interfaces[0]:
            if 0 not in out.keys():
                out[0] = traj

        # 0+ intf
        else:
            omax = np.max(x[:, 1])
            valid_in = False
            for i, interface in enumerate(interfaces[:-1]):
                if omax > interface:
                    valid_in = i + 1
            if valid_in > max_valid_in:
                max_valid_in = valid_in
            if omax > global_maxop:
                global_maxop = omax
            while valid_in:
                # if paths allready exist in high ensembles
                if valid_in in out.keys():
                    valid_in -= 1
                else:
                    out[valid_in] = traj
                    valid_in = 0

    # if we miss some lower ensembles, add to
    # them the paths from the higher ensembles
    for i in range(len(interfaces) + 1):
        if i not in out.keys():
            for j in range(i + 1, len(interfaces) + 1):
                if j in out.keys():
                    out[i] = out[j]
                    print(f"* copied path from ens{j} to ens{i}")


    # Check if we have paths in all ensembles
    for i in range(len(interfaces)):
        assert (
            i in out.keys()
        ), f"* Did not find any paths in ensemble {i}\
    that cross the corresponding interface"


    for i, traj in zip(out.keys(), out.values()):
        shutil.copytree(traj, out_dir / str(i))

    print(f"\nAll done! Created folder {out_dir} with new initial paths.")

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
