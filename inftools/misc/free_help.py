from scipy.integrate import simpson

def recursive_simpson(data, axes):
    """
    Recursively apply Simpson's rule to all dimensions of data.
    """
    assert data.ndim == len(axes)
    if data.ndim == 1:
        return simpson(data, axes[0])
    else:
        integ = simpson(data, axes[-1], axis=-1)
        return recursive_simpson(integ, axes[:-1])


def plot_FE(fn_out, cvs, fe, fe_units):
    import matplotlib.pyplot as plt
    
    fig, ax = plt.subplots()
    if len(cvs) == 1:
        ax.plot(cvs[0], fe)
        ax.set_xlabel(r"$q_1$")
        ax.set_ylabel(f"Free energy ({fe_units})")

    elif len(cvs) == 2:
        levels = 20
        contour = ax.contourf(cvs[0], cvs[1], fe.T, cmap=plt.get_cmap("viridis"), levels=levels)
        ax.contour(cvs[0], cvs[1], fe.T, colors='black', alpha=0.5, linewidths=0.75, zorder=100, levels=levels)
        cbar = fig.colorbar(contour, ax=ax)
        cbar.set_label(f"Free energy ({fe_units})")
        ax.set_xlabel(r"$q_1$")
        ax.set_ylabel(r"$q_2$")

    else:
        raise NotImplementedError("Free energies along more than 2 collective variables are not supported at the moment.")
    fig.savefig(fn_out)