from typing import Annotated
import typer

def estimate_interfaces2(
    i: Annotated[str, typer.Option("-i", help="The Pcross.txt file")],
    num: Annotated[int, typer.Option("-num", help="Number of interfaces")] = None,
    ploc: Annotated[float, typer.Option("-ploc", help="Place interfaces to match ploc")] = None,
    i0: Annotated[float, typer.Option("-i0", help="Position of the first interface")] = None,
    iN: Annotated[float, typer.Option("-iN", help="Position of the last interface")] = None,
    plot: Annotated[bool, typer.Option("-plot", help="Time from nskip")] = False,
    ):
    """
    Estimate the interface positions from Pcross.txt by linearizing the curve.

    TODO:
        * how to handle large drops in the curve? e.g. more than pL?
        * when we have little data, the curve will be very blocky

    """
    import numpy as np
    from inftools.misc.infinit_helper import (
        smoothen_pcross,
        estimate_interface_positions)
    if plot:
        import matplotlib.pyplot as plt

    x = np.loadtxt(i)
    x0 = x[:,0]
    p0 = x[:,1]

    # trim endpoints if i0 or iN given
    iN = x0[-1] if iN is None else iN
    i0 = x0[0] if i0 is None else i0
    last_idx = np.where(x0 <= iN)[0][-1]
    first_idx = np.where(x0 >= i0)[0][0]
    x = x0[first_idx:last_idx]
    p = p0[first_idx:last_idx]
    p = p/p[0]

    # linearalize the curve
    x, p = smoothen_pcross(x, p)
    if ploc is not None and num is None:
        interfaces = estimate_interface_positions(x, p, ploc)
        interfaces = list(interfaces) + [iN]

    elif num is not None and ploc is None:
        ploc = np.exp(np.log(p[-1]) / (num - 1)) - 1e-6
        interfaces = estimate_interface_positions(x, p, ploc)
        interfaces = list(interfaces) + [iN]

    else:
        raise ValueError("Cant have both (or none) of -num and -ploc set!")

    ploc = np.exp(np.log(p[-1]) / (len(interfaces) - 1))
    print(f"[INFO] There are now {len(interfaces)} interfaces.")
    print(f"[INFO] Calculated local crossing probability = {ploc}")

    print("interfaces = [", ", ".join([f"{itf:.6f}" for itf in interfaces]) + "]")

    if plot:
        f, a = plt.subplots(figsize=(8, 4))
        a.plot(x0, p0, markersize=3, label = f"raw Pcross")
        a.plot(x, p, label = "processed Pcross")
        for inter in interfaces:
            a.axvline(inter, c="k", lw=1)
        a.set(yscale="log")
        plt.legend()
        plt.show()



def estimate_interfaces(
    i: Annotated[str, typer.Option("-i", help="The Pcross.txt file")],
    num: Annotated[int, typer.Option("-num", help="Number of interfaces")],
    order: Annotated[int, typer.Option("-order", help="Polynomial order")],
    i0: Annotated[float, typer.Option("-i0", help="Position of the first interface")] = None,
    iN: Annotated[float, typer.Option("-iN", help="Position of the last interface")] = None,
    plot: Annotated[bool, typer.Option("-plot", help="Time from nskip")] = False,

    ):
    """
    Estimate interfaces from Pcross.txt

    A simple script to recaluclate the interface positions
    from an existing crossing proability file.

    A polynomial of
    order <order> is fit in a semi-log plot.
    The first <i0> and last <iN> interface postions are set by the user,
    as well as the number of interfaces <N>.
    As the crossing probaiblity should decrease or remain constant for
    increasing x, we can play with the hyperparameter <alpha> that
    penalizes derivatives greater than zero.
    """
    import matplotlib.pyplot as plt
    import numpy as np
    from scipy.optimize import curve_fit, minimize
    x = np.loadtxt(i)

    if i0 is None:
        i0 = x[0, 0]  # set first interface
    else:
        i0 = float(i0)

    if iN is None:
        iN = float(x[-1, 0])  # set last interface
    else:
        iN = float(iN)
    N = num  # number of interfaces
    alpha = 0  # hyperparameter to penalize positive derivative of pcross, e.g. 50

    last_idx = np.where(x[:, 0] <= iN)[0][-1]
    first_idx = np.where(x[:, 0] >= i0)[0][0]

    y_fit = np.log(x[first_idx:last_idx, 1]) - np.log(x[last_idx, 1])  # y[-1]=0
    x_fit = x[first_idx:last_idx, 0] / (x[last_idx, 0] - x[first_idx, 0])  # x_fit in range (0,1)
    shift = x_fit[0]
    x_fit -= shift

    # first points in pcross should match
    f0 = y_fit[0]

    def fnc(x, *p, f0=f0):
        """

        An Nth order polynomial that crosses fnc(0)=f0 and fnc(1)=0.

        Parameters
        ----------
        x : array
            x-values.
        *p : array
            polynomial coefficients.
        f0 : the crossing point with the x-axis, optional
            DESCRIPTION. The default is f0.
        alpha : hyperparameter, weight of postitive derivative in objective fnc
            DESCRIPTION. The default is alpha.

        Returns
        -------
        y : array
            polynomial evaluated at x.

        """
        y = f0 + -f0 * x
        for order, pi in enumerate(p):
            y += -pi * x + pi * x ** (order + 2)
        return y


    def dy_dp(x, *p, f0=f0):
        dy = np.zeros((x.shape[0], len(p)))  # dyi_dp0, dyi_dp1, ..., dyi_dpN
        for order, pi in enumerate(p):
            dy[:, order] = -x + x ** (order + 2)
        return dy


    def dy_dx(x, *p, f0=f0):
        dy = np.zeros(x.shape[0])  # dyi_dxi
        dy = -f0
        for order, pi in enumerate(p):
            dy += -pi + (order + 2) * pi * x ** (order + 1)
        return dy


    def of(p, x_fit, y_fit, f0, alpha):
        y = fnc(x_fit, *p, f0=f0)
        dy = dy_dx(x_fit, *p, f0=f0)
        idx = np.where(dy > 0)[0]
        return np.sqrt(np.sum((y - y_fit) ** 2)) + alpha * np.sum(dy[idx])

    # some initial guess of parameters
    popt, pcov = curve_fit(fnc, x_fit, y_fit, p0=np.ones(order))

    # actual optimization with penalizing derivatives
    res = minimize(of, popt, args=(x_fit, y_fit, f0, alpha), method="Nelder-mead")

    x_eval = np.linspace(0, 1, 10000)
    y_eval = fnc(x_eval, *res.x)

    # transform back to original coordiantes
    y_plot = x[last_idx, 1] * np.exp(y_eval)
    x_plot = x_eval * (x[last_idx, 0] - x[first_idx, 0]) + x[first_idx, 0]

    # divide by y_plot[0] here in case intf0 is not equal to default
    # then shift back for plotting later
    y0 = y_plot[0]
    y_plot = y_plot/y0
    pcross = y_plot[-1]  # total crossing probability

    pt = np.exp(np.log(pcross) / (N - 1))  # local crossing probability
    print(f"Plocal = {pt}")

    interfaces = [i0]  # first interface
    id_interfaces = [0]
    for i in range(N - 2):
        tmp_id = np.where(y_plot < pt ** (i + 1))[0][0]
        id_interfaces.append(tmp_id)
        posx = x_plot[tmp_id]
        interfaces.append(posx)

    interfaces.append(iN)  # add last interface
    id_interfaces.append(len(x_plot) - 1)
    interfaces = np.round(np.array(interfaces), 9)
    print("interfaces = [", ", ".join([f"{itf:.04f}" for itf in interfaces]) + "]")

    if plot:
        f, a = plt.subplots(figsize=(8, 4))
        a.plot(x[:, 0], x[:, 1], marker="o", markersize=3)
        a.plot(x_plot, y_plot*y0)
        for inter in interfaces:
            a.axvline(inter, c="k", lw=1)
        a.set(yscale="log")
        plt.show()
