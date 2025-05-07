import tomli
import tomli_w
import shutil
import subprocess

import pathlib as pl
import numpy as np
from inftools.tistools.path_weights import get_path_weights
from inftools.tistools.combine_results import combine_data

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
    # current infinit step
    cstep = config1["infinit"]["cstep"]
    tomls = []
    datas = []
    skip = []
    # use inft combine_data with previous combo.txt (if it exists)
    # and current infretis_data.txt file
    if pl.Path(f"combo_{cstep-1}.toml").exists() and pl.Path(f"combo_{cstep-1}.txt").exists():
        tomls += [f"combo_{cstep-1}.toml"]
        datas += [f"combo_{cstep-1}.txt"]
        skip += [0]
    combine_data(
            tomls = tomls + ["restart.toml"],
            datas = datas + [config1["output"]["data_file"]],
            out=f"combo_{cstep}",
            skip= skip + [int(config1["current"]["cstep"]*config1["infinit"]["skip"])],
            )
    # calculate crossing probability for interface estimation
    xp = get_path_weights(
        toml = f"combo_{cstep}.toml",
        data = f"combo_{cstep}.txt",
        nskip = 0,
        outP = "last_infretis_pcross.txt",
        out = "last_infretis_path_weigths.txt",
        overw = True,
        )
    x = xp[:,0]
    p = xp[:,1]

    # don't place interfaces above cap or last interface
    intf_cap = config["simulation"]["tis_set"].get(
            "interface_cap", config["simulation"]["interfaces"][-1]
            )
    last_point = np.where(x>intf_cap)[0]
    if len(last_point)>0:
        x = x[:last_point[0]]
        p = p[:last_point[0]]

    n = config1["runner"]["workers"]

    Ptot = p[-1]
    num_ens = config["infinit"].get("num_ens", False)
    if num_ens:
        pL = Ptot**(1/(num_ens-2))
    else:
        pL = max(config["infinit"]["pL"], Ptot**(1/(2*n)))
    interfaces, pL_used = estimate_interface_positions(x, p, pL)
    config["infinit"]["prev_Pcross"] = pL_used
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

def estimate_interface_positions(x, p, pL):
    """Estimate (binless) interfaces that are equally spaced wrt pL, meaning
    we only **approximately** have a local crossing probability of at least
    pL.

    We interpolate pcross (x) vs orderp (y), such that interp(0.5) gives the
    orderp value that corresponds to P=0.5 (we actually interp -log[pcross]).

    Notes:
        * The last interface is not added (the state B interface), such that the
        probability to reach state B (or the interfa cap) remains pL.

    Returns:
        the interfaces estimated spaced approximately pL apart
        the actual pL separation between interfaces

    """
    # estimate how many interfaces we need
    n_intf = int(np.log(p[-1])/np.log(pL))

    # the actual pL_n_intf to use
    pL_n_intf = p[-1]**(1/n_intf)

    # p needs to be increasing so use -log(p[::-1]) here, meaning we have to
    intf = np.interp([-np.log(pL_n_intf**(i+1)) for i in range(n_intf-1)], -np.log(p), x)

    # add the first interface as well
    intf = [x[0]] + list(intf)

    return intf, pL_n_intf

class LightLogger:
    def __init__(self, fname):
        self.fname = str(pl.Path(fname).resolve())

    def log(self, msg):
        with open(self.fname, "a") as wfile:
            wfile.write(msg + "\n")

