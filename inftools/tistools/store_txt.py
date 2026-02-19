from typer import Option as Opt
from typing import Annotated as And


def store_txt(
    name: And[str, Opt("-name", help="archive file name, name.h5")] = "archive",
    load: And[str, Opt("-load", help="The path to the folder containing the trajectories")] = "load",
    order: And[bool, Opt("-order", help="Set if order.txt is to be saved")] = False,
    energy: And[bool, Opt("-energy", help="Set if energy.txt is to be saved")] = False,
    traj: And[bool, Opt("-traj", help="Set if traj.txt is to be saved")] = False,
    ):
    """Save/Overwrite text files into a h5 file."""
    import importlib.util
    import os
    import pathlib
    import numpy as np
    if not any([order, energy, traj]):
        print("order, energy and traj are all False, so we will not make h5 file.")
        return

    if importlib.util.find_spec("h5py") is None:
        print("h5py is not installed")
        return
    else:
        import h5py

    filename = pathlib.Path(name).with_suffix('.h5')
    # if os.path.exists(filename):
    #     print(f"{name} already exists!")
    #     return 

    path = pathlib.Path(load)
    types = [
        (order, "order.txt"), 
        (energy, "energy.txt"),
        (traj, "traj.txt")
    ]
    dtype_traj = [
        ('time', 'u4'),
        ('trajfile', 'S30'),
        ('index', 'u4'),
        ('vel', 'i1')
    ]

    with h5py.File(filename, "a") as h5:
        for save, ttype in types:
            if not save:
                continue
            dtype = dtype_traj if "traj" in ttype else np.float64
            tfiles = path.rglob(ttype)
            for idx, fname in enumerate(sorted(tfiles, key=lambda p: int(p.parent.name))):
                arr = np.genfromtxt(
                    fname,
                    invalid_raise=False,
                    encoding="utf-8",
                    dtype = dtype,
                )

                pn = fname.parent.name
                grp = h5.require_group(f"{pn}")
                if ttype in grp:
                    del grp[ttype]
                grp.create_dataset(
                    ttype,
                    data=arr,
                    compression="gzip",
                    chunks=True,
                )
            print(f"saved {idx} {ttype} files")
