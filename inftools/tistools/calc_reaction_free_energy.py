from typing import Annotated, Optional
import typer

def calc_reaction_free_energy(
    wham1: Annotated[str, typer.Option("-wham1", help="The wham folder of the A->B reaction")],
    wham2: Annotated[str, typer.Option("-wham2", help="The wham folder of the B->A reaction")],
    out_unit: Annotated[float, typer.Option("-out_unit", help="Multiplicative factor for the free energy values. Default is 1 (k_B*T)")] = 1.0,
    out_unit_name: Annotated[str, typer.Option("-out_unit_name", help="Name of the unit in which the free energy is returned. Default is 'k_B*T'.")] = "k_B*T",
    fn_out: Annotated[Optional[str], typer.Option("-out", help="Name of the free energy output file, default is None (no file is saved)")] = None,
    plot: Annotated[Optional[str], typer.Option("-plot", help="File where a plot of the free energy can be saved, default is None (no plot is saved)")] = None,
):
    """
    Given two infRETIS simulations for an A->B process and the reverse B->A process, 
    for which a conditional free energy has been computed with WHAM analysis,
    merge the two conditional free energies into the actual free energy.
    The binning of the two free energies should be the same.
    """
    from pathlib import Path
    import numpy as np
    from inftools.misc.free_help import recursive_simpson, plot_FE

    assert not (fn_out is None and plot is None), "Please provide out and/or plot arguments"
    if out_unit != 1.0 and out_unit_name == "k_B*T":
        print("WARNING: it seems you have set out_unit but not the respective out_unit_name. The values and labels in the output files WILL NOT MATCH.")
    elif out_unit == 1.0 and out_unit_name != "k_B*T":
        print("WARNING: it seems you have set out_unit_name but not the respective out_unit. The values and labels in the output files WILL NOT MATCH.")

    path_wham_f = Path.cwd() / wham1
    path_wham_b = Path.cwd() / wham2
    if not path_wham_f.exists() or not path_wham_b.exists():
        raise FileNotFoundError("One or both the provided wham folders do not exist!")
    
    fn_prob_f = path_wham_f / "histo_probability.txt"
    fn_prob_b = path_wham_b / "histo_probability.txt"
    if not fn_prob_f.exists() or not fn_prob_b.exists():
        raise FileNotFoundError("I cannot locate the probability histograms in at least one of the wham folders. Have you requested a free energy calculation?")

    prob_f = np.loadtxt(fn_prob_f)
    prob_b = np.loadtxt(fn_prob_b)
    assert prob_f.shape == prob_b.shape, "The probability histograms have non-matching shapes!"

    cvs = []
    x_f = np.loadtxt(path_wham_f / "histo_xval.txt")
    x_b = np.loadtxt(path_wham_b / "histo_xval.txt")
    assert np.array_equal(x_f, x_b), "The conditional free energies must be computed on the same collective variable bins!"
    cvs.append(x_f)

    if prob_f.ndim == 2:
        y_f = np.loadtxt(path_wham_f / "histo_yval.txt")
        y_b = np.loadtxt(path_wham_b / "histo_yval.txt")
        assert np.array_equal(y_f, y_b), "The conditional free energies must be computed on the same collective variable bins!"
        cvs.append(y_f)
    elif prob_f.ndim > 2:
        raise NotImplementedError("Free energies along more than 2 collective variables are not supported at the moment.")

    # extract the rate constants
    k_f = np.loadtxt(path_wham_f / "runav_rate.txt")[-1,3]
    k_b = np.loadtxt(path_wham_b / "runav_rate.txt")[-1,3]
    K_eq = k_f / k_b # equilibrium constant

    # rescale the B->A conditional probability to satisfy the equilibrium constant
    prob_b *= (K_eq / (recursive_simpson(prob_b, cvs) / recursive_simpson(prob_f, cvs)))

    # get the total free energy 
    prob_tot = prob_f + prob_b 
    fe = - out_unit * np.log(prob_tot)
    fe -= np.min(fe) # set the minimum to zero

    if fn_out is not None:
        np.savetxt(Path.cwd() / fn_out, fe, header=f"Free energy values [{out_unit_name}]")
        print(f"Saved free energy values in {fn_out}")
    
    if plot is not None:
        plot_FE(Path.cwd() / plot, cvs, fe, out_unit_name)
        print(f"Saved free energy plot in {plot}")
