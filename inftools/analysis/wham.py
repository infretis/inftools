from typing import Annotated as Atd
import typer
from typer import Option as Opt
import importlib

import importlib.util

import sys


def wham(
    toml: Atd[str, Opt("-toml", help="The infretis .toml file")] = "infretis.toml",
    data: Atd[str, Opt("-data", help="The infretis_data.txt file")] = "infretis_data.txt",
    nskip: Atd[int, Opt("-nskip", help="Number of lines to skip in infretis_data.txt")] = 100,
    lamres: Atd[float, Opt("-lamres", help="Resolution along the orderparameter, (intf1-intf0)/10)")] = None,
    nblock: Atd[int, Opt("-nblock", case_sensitive=False, help="Minimal number of blocks in the block-error analysis")] = 5,
    folder: Atd[str, Opt("-folder", help="Output folder")] = "wham",
    fener: Atd[bool, Opt("-fener", help="If set, calculate the conditional free energy. See Wham_")] = False,
    nbx: Atd[int, Opt("-nbx", help="Number of bins in x-direction when calculating the free-energy")] = 100,
    nby: Atd[int, Opt("-nby", help="Same as -nbx but in y-direction")] = None,
    minx: Atd[float, Opt("-minx", help="Minimum orderparameter value in the x-direction when calculating FE")] = 0.0,
    maxx: Atd[float, Opt("-maxx", help="Maximum orderparameter value in the x-direction when calculating FE")] = 100.0,
    miny: Atd[float, Opt("-miny", help="Same as -minx but in y-direction")] = None,
    maxy: Atd[float, Opt("-maxy", help="Same as -maxx but in y-direction")] = None,
    xcol: Atd[int, Opt("-xcol", help="What column in order.txt to use as x-value when calculating FE")] = 1,
    ycol: Atd[int, Opt("-ycol", help="Same as -xcol but for y-value")] = None,
    ):
    """Run Titus0 wham script."""
    import os

    # import tomli
    from inftools.analysis.Wham_Pcross import run_analysis
    # P = importlib.import_module("inftools.analysis.Wham_Pcross")
    # run_analysis = getattr(P, "run_analysis")
    # run_analysis = lazy_import("inftools.analysis.Wham_Pcross.run_analysis")


    inps = {
        "toml": toml,
        "data": data,
        "nskip": nskip,
        "lamres": lamres,
        "nblock": nblock,
        "fener": fener,
        "folder": folder,
        "histo_stuff":{
            "nbx":nbx, "minx":minx, "maxx":maxx, "xcol":xcol,
            "nby":nby, "miny":miny, "maxy":maxy, "ycol":ycol,
            }
    }
    # load input:
    if os.path.isfile(inps["toml"]):
        with open(inps["toml"], mode="rb") as read:
            config = tomli.load(read)
    else:
        print("No toml file, exit.")
        return
    inps["intfs"] = config["simulation"]["interfaces"]
    if "tis_set" in config["simulation"]:
        inps["lm1"] = config["simulation"]["tis_set"].get("lambda_minus_one", None)
    else:
        inps["lm1"] = None

    if inps["lamres"] is None:
        inps["lamres"] = (inps["intfs"][1] - inps["intfs"][0]) / 10

    if inps["fener"]:
        inps["trajdir"] = config["simulation"].get("load_dir", "load")

    run_analysis(inps)
