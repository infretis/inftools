from typing import Annotated as Atd
from typer import Option as Opt

def pcross_errors(
    toml: Atd[str, Opt("-toml", help=".toml file to read interfaces")] = "infretis.toml",
    data: Atd[str, Opt("-data", help="the infretis_data.txt file")] = "infretis_data.txt",
    overw: Atd[bool, Opt("-O", help="Overwrite existing temporary files")] = False,
    ):
    """This script strips certain ensembles from the infretis_data.txt file.

    Given an infretis_data.txt file, this script create N_interfaces-1
    additional data files and corresponding toml files in a way that allows
    estimating the pcross error at each interface location.

    For example, given a simulation with interfaces = [0,1,2,3], this script
    creates 3 new data and toml files:
    * file 1: data with ensemble [0-] [0+] and toml interfaces at [0,1]
    * file 2: data with ensemble [0-] [0+] [1+] and toml interfaces at [0,1,2]
    * file 3: data with ensemble [0-] [0+] [1+] [2+] and toml interfaces at [0,1,2,3]

    This allows estimating the error in the Ptot at interface i using 'inft wham'.

    Note: This file requires you to run 'inft wham' N_interfaces-1 times.
    """
    import os

    import numpy as np
    import tomli
    import tomli_w
    import tqdm

    # load toml and infretis_data
    with open(toml, "rb") as rfile:
        toml = tomli.load(rfile)
    interfaces = np.array(toml["simulation"]["interfaces"])
    data = np.loadtxt(data, dtype=str)


    # check that we do not overwrite any data that we should not overwrite
    for i, intfi in enumerate(interfaces[1:]):
        outf = f"infretis_{i}_interface_data.txt"
        if os.path.exists(outf) and not overw:
            raise ValueError(f"{outf} exists! Use -O option to overwrite")
        outt = f"infretis_{i}_interface.toml"
        if os.path.exists(outt) and not overw:
            raise ValueError(f"{outt} exists! Use -O option to overwrite")


    for i, intfi in enumerate(interfaces[1:]):
        outf = f"infretis_{i}_interface_data.txt"
        outt = f"infretis_{i}_interface.toml"
        with open(outf, "w") as wfile:
            idx = np.where(intfi>=interfaces)[0]
            relevant_columns = [0,1,2] + [3 + i for i in idx] + [3 + len(interfaces) + i for i in idx]
            #print(intfi, relevant_columns)
            with open(outt, "wb") as wtoml:
                relevant_interfaces = [interfaces[0]] + list(interfaces[idx][1:])
                #print(relevant_interfaces)
                toml["simulation"]["interfaces"] = relevant_interfaces
                tomli_w.dump(toml, wtoml)
            data_clean = data[:, relevant_columns]
            for data_line in data_clean:
                string = "\t" + "\t".join(data_line) + "\t\n"
                wfile.write(string)
            # print(string)
            print(f"inft wham -data {outf} -toml {outt} -folder wham_interface{i}")
