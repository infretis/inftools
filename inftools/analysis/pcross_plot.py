import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.ticker import LogLocator

def plot_combined_pcross(wham_folder):
    """
    Generate a figure with three subplots:
    (a) Pcross vs lambda with interface markers
    (b) runav_Pcross vs path count
    (c) Error vs block length for WHAM method

    Parameters:
    - wham_folder (str): Path to the directory containing input WHAM files.
    - output_file (str): Filename to save the generated plot.
    """

    # --- Matplotlib configuration ---
    rcParams.update({
       # "text.usetex": True,
       # "font.family": "serif",
       # "font.serif": ["Times"],
        "axes.labelsize": 10,
        "font.size": 10,
        "legend.fontsize": 10,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "lines.linewidth": 1,
        "lines.markersize": 2
    })

    # --- Load interface positions from file ---
    with open(os.path.join(wham_folder, "interfaces.txt"), 'r') as f:
        first_line = f.readline().strip()
    interfaces = [float(x) for x in first_line.split('\t')[1:]]

    # --- Set up figure with 3 vertically stacked subplots ---
    fig, axs = plt.subplots(3, 1, figsize=(3.35, 6))

    # ======================
    # Plot (a): Pcross
    # ======================
    data = pd.read_csv(os.path.join(wham_folder, "Pcross.txt"), sep=r'\s+')
    ax = axs[0]
    ax.plot(data['#lam'], data.iloc[:, 1], label=data.columns[1])

    for x in interfaces:
        ax.axvline(x=x, color='black', linestyle='--', alpha=0.4)

    ax.set_xlabel(r'$\lambda$ (Å)')
    ax.set_ylabel(r'$P_A(\lambda|\lambda_A)$')
    ax.set_yscale('log')
    ax.tick_params(axis='both')

    # ======================
    # Plot (b): runav_Pcross
    # ======================
    data = pd.read_csv(os.path.join(wham_folder, "runav_Pcross.txt"), sep=r'\s+')
    crossing_prob = data.iloc[-1, 1]
    ax = axs[1]
    ax.plot(data['#counter'], data.iloc[:, 1])
    ax.set_xlabel('Accepted path count')
    ax.set_ylabel(r'$P_A(\lambda_B|\lambda_A)$')
    ax.set_yscale('log')
    ax.tick_params(axis='both')

    # ======================
    # Plot (c): errPtot
    # ======================
    data = pd.read_csv(os.path.join(wham_folder, "errPtot.txt"), sep=r'\s+', comment='#', header=None)
    data.columns = ['Block', 'PM', 'WHAM']
    ax = axs[2]
    ax.plot(data['Block'], data['WHAM'])

    second_half = data['WHAM'].iloc[len(data)//2:]
    mean_val = second_half.mean()
    x_vals = data['Block'].iloc[len(data)//2:]
    ax.plot(x_vals, [mean_val] * len(x_vals), color='red', linestyle='--', zorder=3)

    x_annotate = x_vals.iloc[len(x_vals)//2]
    ax.annotate(f'{mean_val:.2f}', xy=(x_annotate, mean_val), xytext=(0, 30),
                textcoords='offset points', ha='center',
                arrowprops=dict(arrowstyle='->', color='red'))

    ax.set_xlabel('Block length')
    ax.set_ylabel(r'Relative error $P_A(\lambda_B|\lambda_A)$')
    ax.tick_params(axis='both')

    axs[0].set_title(f"Pcross = {crossing_prob:.2e} ± {mean_val * 100:.0f} \%")

    # ======================
    # Save figure
    # ======================
    plt.subplots_adjust(hspace=0.5)
    filename = os.path.join(wham_folder, "Combined_Pcross.png")
    plt.savefig(filename, dpi=500, bbox_inches='tight')
    plt.close()
