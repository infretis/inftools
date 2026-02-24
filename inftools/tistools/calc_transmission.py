from typer import Option as Opt
from typing import Annotated as And


def calc_transmission(
    ts: And[float,  Opt("-ts", help="The proposed transition state value")],
    dim: And[int,  Opt("-dim", help="The CV in order.txt to calculate the transmission coeff for. Note: Reactive trajs are still defined by the original order parameter.")] = 1,
    toml: And[str,  Opt("-toml", help="The infretis.toml file")] = "infretis.toml",
    pweights: And[str,  Opt("-pweights", help="To get the equilibrium average. This file can be obtained from 'inft get_path_weights'")] = "path_weights.txt",
    traj: And[str,  Opt("-traj", help="location to the order.txt files")] = "load",
    h5: And[str,  Opt("-h5", help="location to the archive.h5 files")] = None,
    nskip: And[int,  Opt("-nskip", help="Skip initial order.txt files")] = 100,
    ):
    """Estimates the transimssion coefficient for a specific CV.
    
    Defined as in Figure 8 in 10.1021/acs.jctc.5c01814
    """
    import numpy as np
    import os
    import matplotlib.pyplot as plt
    from inftools.misc.infinit_helper import read_toml

    def cross_cnt(x, y, ts0):
        """Count positive crossings with ts0"""
        crossings = 0
        for y0, y1 in zip(y[:-1], y[1:]):
            if y1 > ts0 and ts0 > y0:
                crossings += 1
        return crossings

    # get last interface
    intfb = read_toml(toml)["simulation"]["interfaces"][-1]

    # load pdata
    data = np.loadtxt(pweights)[nskip:]

    # consider only trajectories that has crossed the transition state
    data = data[data[:, 1] >= ts]

    if h5 is not None:
        if not os.path.exists(h5):
            print("h5 file does not exist!")
            return
        import h5py
        h5file = h5py.File(h5, "r")

    pweights, preacts, pcrosses = [], [], []
    for pdata in data:
        pn, max_op, pweight = pdata
        
        # get data from .h5 or txt..
        ppath_txt = os.path.join(traj, str(int(pn)), "order.txt")
        ppath_h5  = f"{int(pn)}/order.txt"
        if h5 and ppath_h5 in h5file:
            order = h5file[ppath_h5]
        elif os.path.exists(ppath_txt):
            order = np.loadtxt(ppath_txt)
        else:
            continue

        crosses = cross_cnt(order[:, 0], order[:, dim], ts)
        react = int(max_op > intfb)

        preacts.append(react)
        pcrosses.append(crosses)
        pweights.append(pweight)

    pweights = np.array(pweights)
    preacts = np.array(preacts)
    pcrosses = np.array(pcrosses)

    # transmission coefficient
    tcoeff = np.sum(pweights*preacts)/np.sum(pweights*pcrosses)
    print("transmission coefficient:", tcoeff)
