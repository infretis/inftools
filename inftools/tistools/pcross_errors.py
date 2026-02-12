from typing import Annotated as Atd
from typer import Option as Opt

def pcross_errors(
    pw_file: Atd[str, Opt("-pw", help="path_weights.txt file")] = "path_weights.txt",
    toml: Atd[str, Opt("-toml", help=".toml file to read interfaces")] = "infretis.toml",
    nbins: Atd[int, Opt("-nbins", help="number of bins along lambda, dont pick to large values (slow)")] = 20,
    nblock: Atd[int, Opt("-nblock", help="minimum number of blocks for error analysis, see wham")] = 5,
    plot: Atd[bool, Opt(help="plot results")] = True,
    ref: Atd[str, Opt("-ref", help="Reference wham/ folder for comparison") ] = "",
    ):
    """Estimate error-bars of a Pcross curve from a path weights file.

    The error in a property can be calculated from the running average/estimate
    of the specific property, as illustrated in Appendix A of:

    W. Vervust, D. T. Zhang, A. Ghysels, S. Roet, T. S. van Erp, E. Riccardi,
    J. Comput. Chem. 2024, 45(15), 1224. https://doi.org/10.1002/jcc.27319

    Here, we calculate the running average of the crossing probability at -nbins
    locations along the order parameter, starting at 'interfaces[0]' and ending at
    'interfaces[-1]', which are read from -toml. This is done by using the
    WHAM weights wi for each path i and the maximum order of the same path.
    """
    import tomli
    import numpy as np
    from inftools.analysis.rec_error import rec_block_errors
    if plot:
        import matplotlib.pyplot as plt
        f,ax = plt.subplots(2,2,figsize=(6.5,6))
        a0,a1,a2,a3 = ax[0,0],ax[0,1],ax[1,0],ax[1,1]

    # load path weights
    x = np.loadtxt(pw_file)
    maxop = x[:,1]
    pw = x[:,-1]

    # plot reference values from Wham
    if ref:
        P0 = np.loadtxt(f"{ref}/Pcross.txt")
        runavP0 = np.loadtxt(f"{ref}/runav_Pcross.txt")
        errP0 = np.loadtxt(f"{ref}/errPtot.txt")

    with open(toml, "rb") as rfile:
        config = tomli.load(rfile)

    intf0 = config["simulation"]["interfaces"][0]
    intf1 = config["simulation"]["interfaces"][-1]

    #exclude
    dlambda = np.linspace(intf0, intf1, nbins)
    bin_counts = np.digitize(maxop, dlambda)

    # identiy matrix is used to select rows that match bin locations
    M = np.tri(nbins+1, dtype=int)[bin_counts]

    # average Pcross at each binpoint
    Pcross = np.sum(M*x[:,-1:],axis=0)

    # all paths should have a maxop >= interfaces[0]
    # if not the below assertion fails
    assert Pcross[0] == Pcross[1], "Do some paths have maxop<interfaces[0]?"
    Pcross = Pcross[1:]

    # now calculate running average Pcross at each binpoint
    cumsum = np.cumsum(x[:,-1]).reshape(-1,1)
    runavPcross = np.cumsum(M*x[:,-1:],axis=0)/cumsum

    # allocate array for the errors
    errs = np.zeros(nbins)
    for i,pi in enumerate(runavPcross.T[1:]):
        half_average_error, statistical_ineficiency, relative_error = rec_block_errors(pi, nblock)
        errs[i] = half_average_error
        if plot:
            a1.plot(relative_error,lw=0.75)
            a3.plot(pi)

    if ref:
        print("="*40,)
        print(f"[INFO] Values from the reference {ref} might deviate if the wham command had different options.")
        print("[INFO] E.g., -nblock (or -nskip from 'inft get_path_weights')")

    if plot:
        a0.plot(dlambda, Pcross,marker="o",markersize=2.5)
        a0.fill_between(dlambda, Pcross*(1-errs), Pcross*(1+errs),alpha=0.5,color="C1")
        a2.plot(dlambda[1:], errs[1:], marker="o")
        # set axis labels, etcs
        a0.set(yscale = "log", xlabel = "x, lambda", ylabel = "P(x|x0) +/- error in P(x|x0)", title="Pcross interval")
        a1.set(xlabel="block length - 1", ylabel = "relative error P(xi|x0)",title="Block errors at bin loc")
        a2.set(xlabel="x, lambda", ylabel="relative error",title="Error vs lambda")
        a3.set(xlabel="nr of accepted plus paths", ylabel="P(xi|x0)",title="runav Pcross at bin loc",yscale="log")
        # plot reference values from wham directory 'ref'
        if ref:
            a0.plot(P0[:,0],P0[:,-1],c="k",label="reference",ls="--")
            a3.plot(runavP0[:,-1],c="k",label="reference",ls="--")
            a1.plot(errP0[:,0]-1,errP0[:,-1],label="reference",c="k",ls="--")
            a0.legend()
        f.tight_layout()
        plt.show()
    else:
        return dlambda, Pcross, errs
