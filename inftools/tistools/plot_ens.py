from typing import Annotated as Atd
import typer



def plot_ens(
    toml: Atd[str, typer.Option("-toml")] = "restart.toml",
    data: Atd[str, typer.Option("-data")] = "",
    skip: Atd[bool, typer.Option("-skip", help="skip initial load paths")] = False,
    save: Atd[ str, typer.Option("-save", help="save with scienceplots.") ] = "no",
    pp: Atd[ bool, typer.Option("-pp", help="partial paths version") ] = False,
    cap: Atd[ int, typer.Option("-cap", help="max paths plotted per ens") ] = 100,
    time: Atd[float, typer.Option("-time", help="divide dt by an amount") ] = 1,
    load: Atd[str, typer.Option("-load",
        help = "the path directory, reads load_dir from toml if not given") ] = "",
):
    """Plot sampled ensemble paths with interfaces"""
    import os
    import matplotlib.pyplot as plt
    import numpy as np
    import tomli

    from inftools.tistools.max_op import COLS
    if save != "no":
        import scienceplots
        plt.style.use("science")
        plt.figure(figsize=(14, 10))

    # Read toml info
    with open(toml, "rb") as toml_file:
        toml = tomli.load(toml_file)
    intf = toml["simulation"]["interfaces"]
    if not toml["output"].get("data_file", False) and not data:
        exit("Supply a infretis_data.txt file with -data")
    elif data:
        print(f"Using {data}")
        datafile = data
    else:
        datafile = toml["output"]["data_file"]
    if not load:
        load_dir = toml["simulation"]["load_dir"]
    else:
        load_dir = load

    plt.title("intfs: " + " ".join([str(i) for i in intf]))
    plt.axhline(intf[0], ls="--", color="k", alpha=0.5)
    plt.axhline(intf[-1], ls="--", color="k", alpha=0.5)

    cut = 5 if pp else 3

    acclen = 1
    for ens in list(range(len(intf))):
        cnt = 0
        if ens not in (0, len(intf) - 1):
            plt.axhline(intf[ens], color=f"{COLS[(ens+1)%len(COLS)]}", alpha=1.0)

        with open(datafile) as read:
            for line in read:
                if "#" in line:
                    continue
                rip = line.rstrip().split()
                frac = rip[cut : cut + len(intf)]
                if "-" in frac[ens]:
                    continue
                pnum = rip[0]
                pnum_f = f"{load_dir}/{pnum}/order.txt"
                if not os.path.isfile(pnum_f) or int(pnum) < len(intf):
                    continue
                data = np.loadtxt(pnum_f)
                plt.plot((data[:, 0] + acclen)/time, data[:, 1], color=f"{COLS[ens%len(COLS)]}")
                acclen += len(data[:, 0])
                cnt += 1
                if cnt == cap:
                    break

    plt.xlabel(r"Time")
    plt.ylabel(r"Order Parameter")
    if save != "no":
        plt.savefig(save)
    else:
        plt.show()
