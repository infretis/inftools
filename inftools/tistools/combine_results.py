from typing import Annotated
import typer

def combine_data(tomls: Annotated[list[str], typer.Option("-tomls", help="tomls for all simulations")],
                 datas: Annotated[list[str], typer.Option("-datas", help="data files for all simulations")],
                 skip: Annotated[list[int], typer.Option("-skip", help="skip initial lines for simulations.")] = [100],
                 out: Annotated[str, typer.Option("-out", help="name for output .txt/toml file.")] = "combo"
):
    """Combine different infretis simulations.

    A relatively convoluted task given that we have to merge
    simulation ensembles together from small to big (if I
    assume correctly).

    Require input of their individual toml and data files,
    in their respective order, e.g.

    -tomls sim1 -tomls sim2 -datas data1 -datas data2.

    -skip can be either one value, or specific skip value
    must be specified for each simulation data.
    """

    import numpy as np
    import tomli
    import tomli_w

    from inftools.misc.data_helper import data_reader


    # do some initial checks
    assert len(set(tomls)) == len(tomls) == len(set(datas)) == len(datas)

    # check that len(skip) is either 1 (so same skip for all),
    # or 1 for each sim
    assert len(skip) in (1, len(set(tomls)))
    if len(set(tomls)) > 1 and len(skip) == 1:
        skip = skip*len(set(tomls))

    # initialize some variables
    sims = {}
    intfa, intfb, intfs = [], [], []

    # loop through individual simulations
    for idx, (toml, data) in enumerate(zip(tomls, datas)):
        sims[idx] = {"toml": toml, "data": data}
        with open(toml, "rb") as rfile:
            sims[idx]["conf"] = tomli.load(rfile)
            sims[idx]["intf"] = sims[idx]["conf"]["simulation"]["interfaces"]
            sims[idx]["paths"] = data_reader(data)

            # Initial column placements.
            sims[idx]["cols"] = list(range(len(sims[idx]["intf"])))

            # store their intfa and intfb for assertion
            intfa.append(sims[idx]["intf"][0])
            intfb.append(sims[idx]["intf"][-1])
            intfs += sims[idx]["intf"]

    # Check that sim interfaces start and end the same place.
    assert len(set(intfa)) == 1
    assert len(set(intfb)) == 1

    # create col_order, i in (0, 1) are the 0^- and 0^+ ensembles
    tot_intfs = list(sorted(set(intfs)))
    col_order = {intf: i for i, intf in enumerate(tot_intfs[1:-1], 2)}

    # write combo toml
    with open(f"{out}.toml", "wb") as f:
        tomli_w.dump({"simulation":{"interfaces": tot_intfs}}, f)

    # create combo col for individual sims
    datalines = [[] for _ in range(len(sims))]
    for idx in sims.keys():
        sim_intfs = sims[idx]["intf"]
        cols = list(range(len(sims[idx]["intf"])))
        recol = [0, 1] + [col_order[i] for i in sim_intfs[1:-1]]
        col_dic = {i: j for i, j in zip(cols, recol)}

        for line_num, path in enumerate(sims[idx]["paths"]):
            if line_num < skip[idx]:
                continue
            frac, weig = [], []

            # create new pdic with combo cols:
            pcols = {col_dic[i]:value for i, value in path["cols"].items()}

            # fill up with weights or "----"
            for col in range(len(tot_intfs)):
                if col in pcols:
                    frac.append(pcols[col][0])
                    weig.append(pcols[col][1])
                else:
                    frac.append("----")
                    weig.append("----")
            string = ""
            string += f"\t{path['pn']}\t"
            string += f"{path['len']}" + "\t"
            string += f"{path['max_op']}" + "\t"
            datalines[idx].append(
                    string + "\t".join(frac) + "\t" + "\t".join(weig) + "\t\n"
            )

    # Zip data from individual simulations together for
    # more accurate block error estimates
    with open(f"{out}.txt", "w") as write:
        data_lens = [len(i) for i in datalines]
        data_tot = sum(data_lens)
        modulos = [np.round(data_tot/i) for i in data_lens]
        idx = 0
        while sum([len(i) for i in datalines]) > 0:
            towrite = []
            for mod, data in zip(modulos, datalines):
                if idx%mod == 0 and data:
                    towrite.append(data.pop(0))

            for line in towrite:
                write.write(line)
            idx += 1

    print("The following command can now be ran:")
    print("inft wham -toml combo.toml -data combo.txt -folder combo")
