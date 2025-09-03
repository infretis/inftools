import argparse
from typing import Annotated
import typer

COLS = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
        '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
        '#bcbd22', '#17becf']


def plot_max_op(
    ensp: Annotated[int, typer.Option("-ensp")],
    data: Annotated[str, typer.Option("-data")] = "infretis_data.txt",
    toml: Annotated[str, typer.Option("-toml")] = "infretis.toml",
    weights: Annotated[bool, typer.Option("-weights")] = False,
    ):
    """Plots the max order parameter for paths in an certain ensemble."""
    from inftools.misc.tomlreader import infretis_data_reader
    import numpy as np
    import matplotlib.pyplot as plt
    import tomli

    # read interfaces from .toml file
    with open(toml, "rb") as toml_file:
        toml_dict = tomli.load(toml_file)
    intfs = toml_dict["simulation"]["interfaces"]

    plt.axhline(intfs[0], color=f'k')
    plt.axhline(intfs[-1], color=f'k')
    plt.axhline(intfs[ensp], color=f'C0')
    plt.axhline(intfs[ensp+1], color=f'C1')
    print('# interface', intfs[ensp])

    ens_dic = infretis_data_reader(data)
    data = ens_dic[ensp+1]
    for i in range(2):
        we = np.array(data['w'])/np.array(data['haw'])
        we = we/max(we)
        maxidx = np.argmax(we)
        data['op'].pop(maxidx)
        data['haw'].pop(maxidx)
        data['w'].pop(maxidx)
    we = np.array(data['w'])/np.array(data['haw'])
    we = (we/max(we))**2

    if weights:
        plt.scatter(list(range(len(data['op']))), data['op'], color=f'C0', alpha=we)
    else:
        plt.scatter(list(range(len(data['op']))), data['op'], color=f'C0')

    plt.show()
