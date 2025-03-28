import typer
from typing import Annotated

def get_path_weights(
        toml: Annotated[str, typer.Option("-toml", help = "The .toml file")] = "infretis.toml",
        data: Annotated[str, typer.Option("-data", help = "The infretis_data.txt file")] = "infretis_data.txt",
        out: Annotated[str, typer.Option("-out", help = "Output .txt of the path weights")] = "path_weights.txt",
        nskip: Annotated[int, typer.Option("-nskip", help = "Skip the first nskip entries of the infretis_data.txt file")] = 0,
        plot: Annotated[bool, typer.Option("-plotP", help="If true plot the binless crossing probability")] = False,
        ):
    """Calculate the unbiased weight for paths in the plus ensembles.

    Can be used to calculate observables as <O> = sum(wi*Oi) or for predictive power.
    """
    import numpy as np
    import tomli
    import os

    if os.path.exists(out):
        raise ValueError(f"Out file {out} exists!")

    # load toml and infretis_data
    with open(toml, "rb") as rfile:
        toml = tomli.load(rfile)
    interfaces = toml["simulation"]["interfaces"]
    data = np.loadtxt(data, dtype=str)

    # drop nskip
    data = data[nskip:]
    # we only need non-zero paths
    non_zero_paths = (data[:,3]=="----")
    data[data=="----"] = "0.0"
    D = {}

    D["pnr"] = data[non_zero_paths,0:1].astype(int)
    D["len"] = data[non_zero_paths,1:2].astype(int)
    D["maxop"] = data[non_zero_paths,2:3].astype(float)
    D["path_f"] = data[non_zero_paths,4:3+len(interfaces)].astype(float)
    D["path_w"] = data[non_zero_paths,4+len(interfaces):3+2*len(interfaces)].astype(float)

    w = D["path_f"]/D["path_w"]
    w[np.isnan(w)]=0
    # Need to scale w such that sum equals the number of (fractional) samples n
    # where n=sum(D['path_f'],axis=0)
    w = w/np.sum(w, axis=0)*np.sum(D["path_f"],axis=0)

    ploc_unscaled = np.zeros(len(interfaces))
    ploc_wham = np.zeros_like(ploc_unscaled)
    ploc_unscaled[0] = 1.0
    ploc_wham[0] = 1.0
    for i,intf_p1 in enumerate(interfaces[1:]):
        h1 = D["maxop"]>=intf_p1
        # nmr of paths crossing lambda_i for each ensemble
        nj = np.sum(w[:,:i+1], axis=0)
        # nmr of paths crossing lambda_i+1 for each ensemble
        njl = np.sum(h1*w[:,:i+1], axis=0)
        ploc_unscaled[i+1] = njl[i]/nj[i]
        ploc_wham[i+1] = np.sum(njl)/np.sum(nj/ploc_wham[:i+1])

    # the unbiased path weights
    A = np.zeros_like(D["maxop"])
    Q = 1/np.cumsum(nj/ploc_wham[:-1])

    for j,pathnr in enumerate(D["pnr"][:, 0]):
        # unbiased weight of each path
        K = min(np.where(D["maxop"][j]>interfaces)[0][-1], len(interfaces)-2)
        A[j] = Q[K]*np.sum(w[j])

    print(f"\nAll done! Weights saved to {out}.")
    np.savetxt(out, np.c_[D["pnr"],D["maxop"],A],
            header="path_nr\tmax_op\tweight",
            fmt=["%8d","%8.4f", "%16.8e"])

    # plot the binless crossing probability
    if plot:
        import matplotlib.pyplot as plt
        idx = np.argsort(D["maxop"].flatten())
        maxop_sorted = D["maxop"][idx].flatten()
        weight_sorted = A[idx].flatten()
        sumw = np.sum(weight_sorted)
        res_y = [1.0]
        res_x = [interfaces[0]]
        for i,moi in enumerate(maxop_sorted):
            res_y.append(np.sum(weight_sorted[i:])/sumw)
            res_x.append(moi)
        plt.plot(res_x, res_y)
        plt.yscale("log")
        for intf in interfaces:
            plt.axvline(intf,c="k")
        plt.show()
