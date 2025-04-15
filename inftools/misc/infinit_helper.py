import tomli
import tomli_w
import shutil
import subprocess

import pathlib as pl
import numpy as np
from inftools.tistools.path_weights import get_path_weights

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
    xp = get_path_weights(
        toml = "restart.toml",
        data = config1["output"]["data_file"],
        nskip = int(config1["infinit"]["cstep"]*config1["infinit"]["skip"]),
        outP = "last_infretis_pcross.txt",
        out = "last_infretis_path_weigths.txt",
        overw = True,
    )
    x = xp[:,0]
    p = xp[:,1]


    if 'x' in config1["infinit"]:
        p0 = config1["infinit"]["p"]
        p0 = np.pad(p0, (0, len(x)-len(p0)))
        w_acc = config1["infinit"]["w_acc"]
        w_n = config1["simulation"]["steps"]

        for idx, ip in enumerate(p):
            p[idx] = ip*w_n/(w_n+w_acc) + p0[idx]*w_acc/(w_n+w_acc)

    x, p = smoothen_pcross(x, p) # also remove NaN and 0 Pcross

    # don't place interfaces above cap or last interface
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
        pL = max(config["infinit"]["pL"], Ptot**(1/(2*n)))
    config["infinit"]["prev_Pcross"] = pL
    interfaces = estimate_binless_interface_positions(x, p, pL)
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
    is at *least* pL, meaning it is a bit higher than pL most of the time
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

def estimate_binless_interface_positions(x, p, pL):
    """Estimate binless interfaces that are equally spaced wrt pL, meaning
    we only **approximately** have a local crossing probability of at least
    pL.

    We interpolate pcross (x) vs orderp (y), such that interp(0.5) gives the
    orderp value that corresponds to P=0.5 (we actually interp -log[pcross]).

    Notes:
        * If the crossing probability drops more than N*pL in one step, we add N-1
        interfaces between the previous and the current interface.

        * The last interface is not added (the state B interface)

    """
    # estimate how many interfaces we need
    n_intf = int(np.log(p[-1])/np.log(pL))
    # p needs to be increasing so use -log(p[::-1]) here, meaning we have to
    ip = np.interp([-np.log(pL**(i+1)) for i in range(n_intf)], -np.log(p), x)
    print("=1"*5,ip, -np.log(p))

    # we have N interfaces on top of each other if there are large drops (>pL)
    intf, N = np.unique(ip, return_counts = True)

    # add the first interface as well
    intf = [x[0]] + list(intf)

    # add N-1 interfaces between values where drops are large
    print("=2"*5, intf)
    for i in range(len(intf)-1):
        Ni = N[i]
        if Ni > 1:
            i0 = intf[i]
            i1 = intf[i+1]
            delta = (i1-i0)/Ni
            intf = intf[:i+1] + [i0 + delta*(j+1) for j in range(Ni-1)] + intf[i+1:]
            print(Ni, i0, i1, delta, intf)
        print(Ni)

    return intf

class LightLogger:
    def __init__(self, fname):
        self.fname = str(pl.Path(fname).resolve())

    def log(self, msg):
        with open(self.fname, "a") as wfile:
            wfile.write(msg + "\n")

