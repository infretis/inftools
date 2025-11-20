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
    interfaces = np.array(config["simulation"]["interfaces"])
    # set some default values first
    pL = config["infinit"].get("pL", 0.3)
    assert 0 <  pL < 1, "pL must be greater than 0 and less than 1"
    config["infinit"]["pL"] = pL

    err_msg = "'nskip' in [infinit] is replaced with 'skip', which is a fraction"
    assert "nskip" not in config["infinit"].keys(), err_msg

    skip = config["infinit"].get("skip", 0.1)
    assert 0 <= skip < 1, "'skip' should be a fraction between 0 and 1!"
    config["infinit"]["skip"] = skip

    # smallest lamres possible is 1e-5 (based on infretis_data.txt resolution)
    lamres = config["infinit"].get("lamres", 0.001)
    assert lamres >= 1e-5, "'lamres' must be >= 1e-5"
    config["infinit"]["lamres"] = lamres

    cstep = config["infinit"].get("cstep",-1)
    config["infinit"]["cstep"] = cstep
    if cstep == -1:
        assert len(interfaces) == 2, "Define 2 interfaces!"
        err_msg = "No 'initial_conf' specified in [infinit]"
        assert "initial_conf" in config["infinit"].keys(), err_msg
        init_conf = pl.Path(config["infinit"]["initial_conf"])
        err_msg = f"'initial_conf' {init_conf.resolve()} does not exist!"
        assert init_conf.exists(), err_msg

    # check that interfaces are multiples of lamres
    # but skip for first iteration
    if cstep not in [-1,0]:
        rounded_intf = np.round(np.round(interfaces/lamres)*lamres, decimals=10)
        err_intf = np.where(rounded_intf != interfaces)[0]
        err_msg = (
            f"Interfaces {interfaces[err_intf]} are not multiples of "
            f"lamres={lamres}"
            )
        assert len(err_intf)==0, err_msg


    steps_per_iter = config["infinit"]["steps_per_iter"]
    config["infinit"]["steps_per_iter"] = steps_per_iter

    assert cstep < len(steps_per_iter), "Nothing to do"
    if not np.all(np.array(steps_per_iter[max(cstep,0):])
            >= config["runner"]["workers"]):
        raise ValueError("The number of infretis steps in steps_per_iter"
                " has to be larger or equal to the number of workers!")

    if config["output"]["delete_old"]:
        keep_max = config["output"].get("keep_maxop_trajs", False)
        assert keep_max, "Cant have delete_old=true and keep_maxop_trajs=false"
    assert config["output"].get("delete_old_all", False) == False
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
    """Run infretis as a subprocess.

    Returns True if successful run, else False.
    """
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
    # check if successful run
    c1 = read_toml("restart.toml")
    if not c1:
        print(" *** Did not find restart.toml after infretisrun.")
        print(" *** Something prevented infretis from starting")
        return False
    completed_steps = c1["current"]["cstep"]
    c1_infinit_cstep = c1["infinit"]["cstep"]
    c0_infinit_cstep = c0["infinit"]["cstep"]
    if not c1_infinit_cstep == c0_infinit_cstep:
        print(f" *** restart.toml {c1_infinit_cstep} != {c0_infinit_cstep} infretis.toml")
        print(" *** 'cstep' in [infinit] differ between restart.toml and infretis.toml")
        return False
    if not c1["infinit"]["steps_per_iter"][c1_infinit_cstep] == completed_steps:
        print(" *** infretisrun did not complete all steps.")
        return False

    return True

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
    if pl.Path(f"combo_{cstep}.toml").exists() and pl.Path(f"combo_{cstep}.txt").exists():
        tomls += [f"combo_{cstep}.toml"]
        datas += [f"combo_{cstep}.txt"]
        skip += [0]
    combine_data(
            tomls = tomls + ["restart.toml"],
            datas = datas + [config1["output"]["data_file"]],
            out=f"combo_{cstep+1}",
            skip= skip + [int(config1["current"]["cstep"]*config1["infinit"]["skip"])],
            )
    # calculate crossing probability for interface estimation
    xp = get_path_weights(
        toml = f"combo_{cstep+1}.toml",
        data = f"combo_{cstep+1}.txt",
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

    # in some cases we might not have enough data to construct a Pcross
    # curve, resulting in p having elements that are 0
    zero_idx = np.where(p>0)[0]

    # if all elements in pcross are 0 continue with same interfaces
    # else remove the elements that are 0 and continue as usual
    if len(zero_idx) == 0 or zero_idx[-1] == 0:
        print("*** Could not construct Pcross curve, all elements are 0.")
        print("*** Continuing with the same interfaces")
        config["infinit"]["prev_Pcross"] = 0.0
        return

    p = p[:zero_idx[-1] + 1]
    x = x[:zero_idx[-1] + 1]

    Ptot = p[-1]
    num_ens = config["infinit"].get("num_ens", False)
    if num_ens:
        interfaces, pL_used = estimate_interface_positions(x, p, num_ens=num_ens)
    else:
        pL = max(config["infinit"]["pL"], Ptot**(1/(2*n)))
        interfaces, pL_used = estimate_interface_positions(x, p, pL=pL)
    config["infinit"]["prev_Pcross"] = pL_used
    intf = list(interfaces) + config["simulation"]["interfaces"][-1:]
    # round interfaces to lambda resolution, and avoid precision errors
    lamres = config["infinit"]["lamres"]
    intf_tmp = np.round(np.round(np.array(intf[1:-1])/lamres)*lamres, decimals=10)
    # remove duplicates if any appear due to rounding of interfaces
    intf_tmp = list(np.unique(intf_tmp))
    # if we suddenly have less workers than interfaces, just return the
    # 'non-rounded' interfaces
    if len(intf_tmp) + 1 < n:
        intf_tmp = intf[1:-1]
        print("* Not rounding interfaces to lamres. Would give n_workers > n_ens.")
    elif len(intf_tmp) < len(intf)-2:
        print(f"* Rounding interfaces to lamres, there are now {len(intf_tmp)+1} plus ensembles, and not {num_ens}")

    config["simulation"]["interfaces"] =intf[:1] +  intf_tmp + intf[-1:]
    config["simulation"]["shooting_moves"] = sh_moves = ["sh", "sh"] + ["wf" for i in range(len(intf)-2)]

def update_toml(config):
    config0 = read_toml("infretis.toml")
    config0["simulation"]["interfaces"] = config["simulation"]["interfaces"]
    config0["simulation"]["shooting_moves"] = config["simulation"]["shooting_moves"]
    config0["infinit"] = config["infinit"]
    shutil.copyfile("infretis.toml", f"infretis_{config['infinit']['cstep']}.toml")
    write_toml(config0, "infretis.toml")

def rename_file(old_file, new_file):
    old_file = pl.Path(old_file)
    new_file = pl.Path(new_file)
    old_file.rename(new_file)

def update_actives_toml(out):
    config0 = read_toml("infretis.toml")
    config1 = read_toml("restart.toml")
    config0["current"] = config1["current"]
    # remove/change some entries in the 'current' section
    config0["current"]["frac"] = {}
    config0["current"]["cstep"] = 0
    if "restarted_from" in config0["current"].keys():
        config0["current"].pop("restarted_from")
    traj_num = config0["current"]["traj_num"]
    # rename active paths for next round. They are not sorted wrt.
    # ensembles so well do this here
    active = [-1 for i in out.values()]
    for i, path_old in out.items():
        new_path_nr = traj_num + i
        active[i] = new_path_nr
        path_new = path_old.parent/f"{new_path_nr}"
        path_new.symlink_to(path_old.resolve(), target_is_directory=True)
    # traj_num should be 1 larger than largest path nr
    config0["current"]["traj_num"] = max(active) + 1
    # update active path list from the path numbers picked with new interfaces
    config0["current"]["active"] = active
    config0["current"]["size"] = len(active)
    # new data file to be written is added to the infretis.toml
    data_file = pl.Path(f"infretis_data_{config0['infinit']['cstep']+1}.txt")
    config0["output"]["data_file"] = str(data_file)

    # avoid overwriting any data files when running infretis with a 'data_file'
    # given, so we rename any existing infretis_data.txt files here
    if data_file.exists():
        for i in range(1000):
            new_data_file = data_file.parent / (data_file.name + f".{i}")
            if not new_data_file.exists():
                data_file.rename(new_data_file)
                print(f"* File renamed from {data_file} to {new_data_file}")
                break

    # add interface column to infretis_data file
    with open(config0["output"]["data_file"], "a") as wfile:
        line = "#intf: " + ",".join(str(i) for i in config0["simulation"]["interfaces"]) + "\n"
        wfile.write(line)
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

def estimate_interface_positions(x, p, pL=0.3, num_ens=False):
    """Estimate (binless) interfaces that are equally spaced wrt pL, meaning
    we only **approximately** have a local crossing probability of at least
    pL. Alternatively, the number of ensembles can be supplied.

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
    if not num_ens:
        num_ens = int(np.log(p[-1])/np.log(pL))

    # the actual pL_n_intf to use
    pL_num_ens = p[-1]**(1/num_ens)

    # p needs to be increasing so use -log(p[::-1]) here, meaning we have to
    intf = np.interp([-np.log(pL_num_ens**(i+1)) for i in range(num_ens-1)], -np.log(p), x)

    # add the first interface as well
    intf = [x[0]] + list(intf)

    return intf, pL_num_ens

class LightLogger:
    def __init__(self, fname):
        self.fname = str(pl.Path(fname).resolve())

    def log(self, msg):
        with open(self.fname, "a") as wfile:
            wfile.write(msg + "\n")

