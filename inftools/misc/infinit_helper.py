import tomli
import tomli_w
import shutil
import subprocess

import pathlib as pl
import numpy as np

PHRASES = [
    ["Infinit mode:", "Engaging endless loops with ∞RETIS",],
    ["Infinitely fast from A to B","with ∞RETIS"],
    ["Infinit mode:", "Engaging endless loops with ∞RETIS",],
    ["Installing the infinite improbability drive ...", "∞RETIS is allready installed!",],
    ["Solving the time warp", "in digital chemical discoveries",],
    ["Performing infinite swaps", "no molecules are harmed in this exchange!",],
    ["Asynchronous mode enabled", "because who has time for synchronization?",],
    ["Pushing for transitions", "molecules, please keep your seatbelts fastened!",],
    ["Propagating forwards and backwards in time", " please hold on to your atoms!",],
    ["Fluxing through rare events",  "with ∞RETIS",],
    ["Taking molecular strolls in parallel", "with ∞RETIS",],
    ["Shooting through the void", "with ∞RETIS",],
    ["Infinit /ˈɪnfɪnət/ (adjective)", "limitless or endless in space, extent, or size"]
]

def set_default_infinit(config):
    """
    TODO:
        check before startup that no folders runi exist.
        Dont move infretis.toml to runi, because if we remove we delete
            the toml. Dont move any infretis.toml, just copy.
        infretis_data_i.txt file is not update in config at runtime

    """
    interfaces = config["simulation"]["interfaces"]
    if not 0 <= config["infinit"]["skip"] < 1:
        raise ValueError("skip should be a fraction between 0 and 1!")
    assert config["infinit"]["pL"] > 0
    cstep = config["infinit"]["cstep"]
    if cstep == -1:
        assert len(interfaces) == 2, "Define 2 interfaces!"
    #config["simulation"]["interfaces"] = [interfaces[0], interfaces[-1]]
    #config["current"]["active"] = [0,1]

    steps_per_iter = config["infinit"]["steps_per_iter"]
    config["infinit"]["steps_per_iter"] = steps_per_iter

    config["infinit"]["cstep"] = cstep
    assert cstep < len(steps_per_iter), "Nothing to do"
    if not np.all(np.array(steps_per_iter[max(cstep,0):])
            >= config["runner"]["workers"]):
        raise ValueError("The number of infretis steps in steps_per_iter"
                " has to be larger or equal to the number of workers!")
    assert config["output"]["delete_old"] == False
    assert config["output"].get("delete_old_all", False) == False
    # update check here based on maximum op value
    # else we need less workers, or lower lamres
    if cstep > 0:
        print(f"Restarting infinit from iteration {cstep}.")

    return config["infinit"]

def read_toml(toml):
    toml = pl.Path(toml)
    if toml.exists():
        with open(toml, "rb") as rfile:
            config = tomli.load(rfile)
        return config
    else:
        return False

def write_toml(config, toml):
    with open(toml, "wb") as wfile:
        tomli_w.dump(config, wfile)

def run_infretis_ext(steps):
    c0 = read_toml("infretis.toml")
    c1 = read_toml("restart.toml")
    if c1 and c0["infinit"]["cstep"] == c1["infinit"]["cstep"] and len(c0["simulation"]["interfaces"])==len(c1["simulation"]["interfaces"]) and np.allclose(c0["simulation"]["interfaces"],c1["simulation"]["interfaces"]):
        print("Running with restart.toml")
        c1["simulation"]["steps"] = steps
        # might have updated steps_per_iter
        c1["infinit"]["steps_per_iter"] = c0["infinit"]["steps_per_iter"]
        write_toml(c1, "restart.toml")
        subprocess.run("infretisrun -i restart.toml", shell = True)
    else:
        print("Running with infretis.toml")
        c0["simulation"]["steps"] = steps
        write_toml(c0, "infretis.toml")
        subprocess.run("infretisrun -i infretis.toml", shell = True)



def update_toml_interfaces(config):
    """Update the interface positions from crossing probability.

    It is based on the linearization of the crossing probability and
    the fact that we want equal local crossing probabilities betewen
    interfaces.
    """
    config1 = read_toml("restart.toml")
    x, p = calc_pcross(
            config1["output"]["data_file"],
            config1["simulation"]["interfaces"],
            config1["infinit"]["lamres"],
            int(config1["infinit"]["skip"]*config1["infinit"]["cstep"],
            )

    if 'x' in config1["infinit"]:
        p0 = config1["infinit"]["p"]
        p0 = np.pad(p0, (0, len(x)-len(p0)))
        w_acc = config1["infinit"]["w_acc"]
        w_n = config1["simulation"]["steps"]

        for idx, ip in enumerate(p):
            p[idx] = ip*w_n/(w_n+w_acc) + p0[idx]*w_acc/(w_n+w_acc)

    x, p = smoothen_pcross(x, p) # also remove NaN and 0 Pcross

    # don't place interfaces above cap
    intf_cap = config["simulation"]["tis_set"].get(
            "interface_cap", config["simulation"]["interfaces"][-1]
            )
    last_point = np.where(x>intf_cap)[0]
    if len(last_point)>0:
        x = x[:last_point[0]]
        p = p[:last_point[0]]

    n = config1["runner"]["workers"]

    # save x, p for next round
    config["infinit"]["x"] = x.tolist()
    config["infinit"]["p"] = p.tolist()
    config["infinit"]["w_acc"] = config1["infinit"].get("w_acc", 0) + config1["simulation"]["steps"]

    Ptot = p[-1]
    num_ens = config["infinit"].get("num_ens", False)
    if num_ens:
        pL = Ptot**(1/(num_ens-2))
    else:
        pL = max(0.3, Ptot**(1/(2*n)))
    config["infinit"]["prev_Pcross"] = pL
    interfaces = estimate_interface_positions(x, p, pL)
    intf = list(interfaces) + config["simulation"]["interfaces"][-1:]
    config["simulation"]["interfaces"] = intf
    config["simulation"]["shooting_moves"] = sh_moves = ["sh", "sh"] + ["wf" for i in range(len(intf)-2)]

def update_folders():
    config = read_toml("infretis.toml")
    old_dir = pl.Path(config["simulation"]["load_dir"])
    new_dir = pl.Path(f"run{config['infinit']['cstep']}")
    if new_dir.exists():
        msg = (f"{str(new_dir)} allready exists! Infinit does not "
                + "overwrite.")
        print(msg)
        if not old_dir.exists():
            print("Did not find {old_dir}.")
            return False

        return True

    shutil.move(old_dir, new_dir)
    return False

def update_toml(config):
    config0 = read_toml("infretis.toml")
    config0["simulation"]["interfaces"] = config["simulation"]["interfaces"]
    config0["simulation"]["shooting_moves"] = config["simulation"]["shooting_moves"]
    config0["infinit"] = config["infinit"]
    shutil.copyfile("infretis.toml", f"infretis_{config['infinit']['cstep']}.toml")
    write_toml(config0, "infretis.toml")

def calc_pcross(ifile, lambda_interfaces, lamres, nskip):
    # the Cxy values of [0+] are stored in the i0plus-th
    # column (first coulumn is counted as column nr 0)
    i0plus = 4

    lambdaA = lambda_interfaces[0]
    lambdaB = lambda_interfaces[-1]
    # number of interfaces
    nintf = len(lambda_interfaces)
    # number of plus-ensembles [i+]
    nplus_ens = nintf - 1
    # The sum of fractional sampling occurrences for each [i+] ensemble
    eta = np.zeros(nplus_ens)
    # Generate a list of values
    lambda_values = np.array([
        i * lamres
        for i in range(round(lambdaA / lamres), round(lambdaB / lamres) + 1)
    ])
    # v_alpha: to become the total crossing probability based on WHAM
    v_alpha = np.zeros(len(lambda_values))
    # This entry is not changed anymore
    v_alpha[0] = 1


    # Open the file
    with open(ifile) as file:
        # Initialize the matrix
        matrix = []
        # Loop through each line in the file
        for line in file:
            # Ignore comment lines
            if line.startswith("#"):
                continue
            # Remove any trailing whitespace or newlines
            line = line.strip()
            # Split the line into a list of values
            values = line.split()
            # Replace any "----" values with 0.0
            values = [float(x) if x != "----" else 0.0 for x in values]
            # Add the row to the matrix
            matrix.append(values)

    # delete the first nskip entries
    del matrix[:nskip]

    # Format of data-file is not yet unweighted with HA-weights
    # column for [0-] ensemble
    i0min = i0plus - 1
    # sum of Pxy for all y being [0-], [0+], [1+] etc
    sumPxy = [0.0] * nintf
    # sum of Pxy after weighting with inverse HA-weights
    sumPxy_afterw = [0.0] * nintf
    for x in matrix:
        for y in range(nintf):
            # index of value that requires unweighting
            y1 = i0min + y
            # index of HA-weight
            y2 = y1 + nintf
            if x[y2] > 0:  # non-zero weight
                sumPxy[y] += x[y1]  # sum before unweighting
                x[y1] /= x[y2]
                sumPxy_afterw[y] += x[y1]  # sum after unweighting
            elif x[y1] > 0:
                print("Division by zero for path x=", x)
                exit()
            # sumPxy and suminvw remain the same.
            # if x[y1]==x[y2]==0, giving 0/0, x[y1] remains zero.
            # Next:
            # we used to store the HA-weights
            # store the running sum of sumPxy in the matrix where
            # running SumPxy values are needed for
            # computing the running WHAM averages
            # HA-weights are no longer needed after this step, while
            x[y2] = sumPxy[y]
    # loop over x completed
    # dvide each column by the average of the inverse HA-weight
    # This gives more comparible eta[y] values between ensembles.
    # For instance if [0+] is based on shooting and [i+] is based on WF
    # This allows to use the standard WHAM procedure
    for y in range(nintf):
        AvinvwHA = sumPxy_afterw[y] / sumPxy[y]
        y1 = i0min + y  # index of [0-], [0+], [1+] etc
        for x in matrix:
            x[y1] /= AvinvwHA

    # This function gives back the WHAM based crossing
    # probabilities at the interfaces and the Q-factor
    # to normailze the v_alpha vector
    def WHAM_PQ(npe, interf, res, eta, v_alpha):
        P = [0.0] * npe  # Crossin probability at lambda_0, lambda_1,
        #  .... lambda_{n-1}: P_A(\lambda_i|\lambda_0) based on WHAM
        Q = [0.0] * npe  # The Q-factor for Whamming (JCTC 2015 Lervik et al)
        invQ = [0.0] * npe  # The inverse Q-factor
        # Initial values
        P[0] = 1.0
        invQ[0] = eta[0]
        if invQ[0] == 0:
            return P, Q
        Q[0] = 1 / invQ[0]
        # Solve other values using recursive relations
        lambdaA = interf[0]
        for i in range(1, npe):  # Loop over all nterfaces after [0+]
            lambda_i = interf[i]
            alpha = round(
                (lambda_i - lambdaA) / res
            )  # the alpha index corresponding to lambda_i
            P[i] = v_alpha[alpha] * Q[i - 1]
            if P[i] == 0:
                return P, Q
            invQ[i] = invQ[i - 1] + (eta[i] / P[i])
            Q[i] = 1 / invQ[i]
        # add the final value of lambda_B
        P.append(v_alpha[-1] * Q[nplus_ens - 1])
        return P, Q

    # Loop over all paths x in the matrix to:
    # Compute eta[i] for all interfaces except the last lambda_B
    # Create the v_alpha() vector. Yet, without the proper normalization
    # Create the local crossing probabilities p_loc[.... ]
    # which is a matrix showing the local crossing probabilities
    # for [0+], [1+] etc
    # Create running average of local crossing probability
    # P_A(\lambda_{i+1}|\lambda_i): run_av_ploc
    # The structure of the data file is the following
    # x[0]=path-index, x[1]=path length, x[2]=lambda-max
    for x in matrix:
        lambdamax = x[2]
        for i in range(nplus_ens):
            Cxy_index = i0plus + i
            Cxy = x[Cxy_index]
            eta[i] += Cxy  # increase eta[i]
            # Determine lower and upper bound for
            # increasing v_alpha- and p_loc values
            lambda_i = lambda_interfaces[i]
            alpha_max = int(np.floor((lambdamax - lambdaA) / lamres))
            alpha_min = round((lambda_i - lambdaA) / lamres)
            # Note: lambda_i-lambdaA)/lamres is an integer as
            # lambda_i and lambdaA should be commensurate with lamres
            if alpha_max > len(v_alpha) - 1:
                alpha_max = len(v_alpha) - 1  # -1 as we start counting from 0
            alpha_min += 1  # v(alpha) at the interface lambda_1, lambda_2
            # etc are determined by the previous [0+], [1+] etc
            for alpha in range(alpha_min, alpha_max + 1):
                v_alpha[alpha] += Cxy

    # Normalize the p_loc arrays such that each local
    # crossing probability of ensemble [i+] starts with 1 at lambda_i
    _, Q = WHAM_PQ(nplus_ens, lambda_interfaces, lamres, eta, v_alpha)

    # we need now a loop ovber all alpha values and determine K(alpha)
    # It is however faster to loop over K(alpha) indexes and split the
    # alpha-loop into parts that have the same K(alpha) value
    for Kalpha in range(nplus_ens):  # Loop over all pluse-ensembles after [0+]
        lambda_i = lambda_interfaces[Kalpha]
        lambda_next = lambda_interfaces[Kalpha + 1]
        alpha_max = round((lambda_next - lambdaA) / lamres)
        # Note again (lambda_next-lambdaA)/lamres is an integer
        # as interfaces should be commensurate with lamres resolution
        alpha_min = (
            round((lambda_i - lambdaA) / lamres) + 1
        )  # +1 because K(lambda) refers to index K with lambda_K < lambda
        # (not less or equal)
        for alpha in range(alpha_min, alpha_max + 1):
            v_alpha[alpha] *= Q[Kalpha]
    # v_alpha now represents the total crosing probability based on WHAM
    return lambda_values, v_alpha



def print_logo(step: int = 0):
    from rich.console import Console
    from rich.text import Text

    console = Console()
    art = Text()
    # add some initial states, a [0+] path and [0-] path,
    # some A-B paths, etcetera etcetera
    if step >= len(PHRASES) or step == -1:
        step = np.random.choice(len(PHRASES))
    str1,str2 = PHRASES[step]
    #art.append("\n")
    art.append("         o             ∞       ____          \n", style="bright_blue")
    art.append("__________\\________o__________/____\\_________\n", style="bright_blue")
    art.append("           \\      /          o      o        \n", style="bright_blue")
    art.append("            \\____/                           \n", style="bright_blue")
    art.append(str1,style="bright_cyan")
    art.append(f"{'_'*(45-len(str1))}\n", style="bright_blue")
    art.append("_____________________________________________\n", style="bright_blue")
    art.append(" _          __  _         _  _               \n", style="bold light_cyan")
    art.append("(_) _ __   / _|(_) _ __  (_)| |_             \n", style="bold bright_yellow")
    art.append("| || '_ \\ | |_ | || '_ \\ | || __|            \n", style="bold bright_magenta")
    art.append("| || | | ||  _|| || | | || || |_             \n", style="bold bright_green")
    art.append("|_||_| |_||_|  |_||_| |_||_| \\__|            \n", style="bold white")
    art.append("______________\\______________________________\n", style="bright_blue")
    art.append("   ∞           o                o            \n", style="bright_blue")
    #art.append("             o                    ",style="bright_blue")
    art.append(f"{str2:>45}", style="bright_cyan")
    console.print(art)

def smoothen_pcross(x, p):
    x_out = x[p!=0]
    p = p[p!=0]
    p = np.log10(p)

    grad = p[1:] - p[:-1]
    cliff_points = np.where(grad!=0)[0]
    # the foot of a cliff is cliff_point + 1
    foot_points = cliff_points + 1
    foot_points = np.hstack((0,foot_points))
    # add last point, its a cliff
    cliff_points = np.hstack((cliff_points,p.shape[0]-1))
    p_out = np.zeros((2, p.shape[0]))
    prev_idx = 0
    for idx in cliff_points:
        delta = idx - prev_idx
        x_ = np.linspace(0, 1, delta + 1)
        p_out[0, prev_idx:idx] = p[prev_idx] - x_[:-1]*(p[prev_idx] - p[idx])
        prev_idx = idx
    p_out[0,idx:] = p[idx]

    prev_idx = 0
    for idx in foot_points:
        delta = idx - prev_idx
        x_ = np.linspace(0, 1, delta + 1)
        p_out[1, prev_idx:idx] = p[prev_idx] - x_[:-1]*(p[prev_idx] - p[idx])
        prev_idx = idx
    p_out[1,idx:] = p[idx]
    return x_out, 10**np.mean(p_out, axis=0)

def estimate_interface_positions(x, p, pL):
    """Place new interfaces based on Pcross.

    The interfaces are placed *evenly* such that the local crossing probability
    is at *least* pL, meaning it is a bit higher than pL most of the time.
    """
    i = 0
    interfaces = [0]
    n_intf = int(np.log(p[-1])/np.log(pL)) + 1
    pL_new = p[-1]**(1/n_intf)
    for intf in range(n_intf):
        idx = np.where(p/p[i] >= pL_new)[0]
        if len(idx) == len(p):
            break
        i = idx[-1]
        if p[-1]/p[i] > pL_new:
            break
        # when we drop more than pL**2 (but not more than pL**3)
        if i in interfaces:
            i = i+1
        interfaces.append(i)
    return x[interfaces]

class LightLogger:
    def __init__(self, fname):
        self.fname = str(pl.Path(fname).resolve())

    def log(self, msg):
        with open(self.fname, "a") as wfile:
            wfile.write(msg + "\n")

