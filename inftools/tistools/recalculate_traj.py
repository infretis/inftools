from typing import Annotated, Tuple
import typer

def recalculate_traj(
    path: Annotated[str, typer.Option("-path", help="path to the path's path number, usually in load/{pn}")] ,
    toml: Annotated[str, typer.Option("-toml")] = "infretis.toml",
    out: Annotated[str, typer.Option("-out", help="the output of the analysis")] = "order_rec.txt",
    top: Annotated[str, typer.Option("-top", help="topology file passed to mdanalysis")] = "",
    mda_format: Annotated[str, typer.Option("-format", help = "The format of the trajectory file (e.g. 'xyz', 'xtc', 'trr')")] = "",
    ):
    """
    Recalculate the orderparamter from a sampled infretis path.
    """

    import os
    import numpy as np
    import tomli

    import MDAnalysis as mda

    from infretis.classes.engines.cp2k import read_xyz_file
    from infretis.classes.engines.gromacs import read_gromos96_file
    from infretis.classes.orderparameter import create_orderparameter
    from infretis.classes.system import System
    from inftools.analysis.gromacs import read_trr_file

    with open(toml, "rb") as toml_file:
        toml_dict = tomli.load(toml_file)
    orderparameter = create_orderparameter(toml_dict)
    tdata = np.loadtxt(os.path.join(path, "traj.txt"), dtype=str)
    files = set(tdata[:, 1])

    if not mda_format:
        mda_format = list(files)[0].split(".")[-1]
        print(f"Assuming '{mda_format}' trajectory format based on extension.")

    if not top:
        top = None

    # assume all trajs are mda readable
    unis = {i:mda.Universe(top, os.path.join(path, "accepted", i), format = mda_format) for i in files}

    # for cp2k: assume no NPT is ran with cp2k, st. default box is set for all
    if "xyz" in [i.split('.')[-1] for i in files]:
        from infretis.classes.engines.cp2k import read_cp2k_box
        cp2k_inp = os.path.join(toml_dict["engine"]["input_path"], "cp2k.inp")
        box, pbc = read_cp2k_box(cp2k_inp)
        for uni in unis.values():
            uni.dimensions = list(box) + [90]*3

    # write to new file
    out_file = os.path.join(path, out)
    for i in range(1, 1000):
        if not os.path.isfile(out_file):
            break
        out_file = os.path.join(path, f"order_rec_{i}.txt")

    with open(out_file, "w") as writefile:
        writefile.write("# step\torder\n")
        for step, fname, index, vel in tdata:
            system = System()
            system.config = (fname, index)
            system.pos = unis[fname].trajectory[int(index)].positions
            system.box = unis[fname].trajectory[int(index)].dimensions[:3]
            op = orderparameter.calculate(system)
            line = f"{step} " + "\t".join([f"{opi}" for opi in op]) + "\n"
            writefile.write(line)
    print(f"[ INFO ] Orderparameter values written to {out_file}\n")

