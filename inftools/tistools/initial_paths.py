import os
from typing import Annotated

import typer
from inftools.misc.infinit_helper import *

# export _TYPER_STANDARD_TRACEBACK=1

def generate_zero_paths(
    conf: Annotated[str, typer.Option("-conf", help="The name (not the path) of the initial configuration to propagate from. inftools will look in the input folder specified in the .toml file.")],
    toml: Annotated[str, typer.Option("-toml",)] = "infretis.toml",
    config : Annotated[str, typer.Option(hidden = True)] = None,
    ):
    """Generate initial paths for the [0-] and [0+]

    ensembles by propagating a single configuration forward
    and backward in time until it crosses the lambda0 interface.
    These can be used to e.g. push the system up the barrier using
    multiple infRETIS simulations."""
    import pathlib as pl

    import numpy as np
    from infretis.classes.engines.factory import create_engines
    from infretis.classes.orderparameter import create_orderparameters
    from infretis.classes.path import Path, paste_paths
    from infretis.classes.repex import REPEX_state
    from infretis.classes.system import System
    from infretis.setup import setup_config


    # make a directory we work from
    tmp_dir = pl.Path("temporary_load/")
    tmp_dir.mkdir(exist_ok = False)
    load_dir = pl.Path("load/")
    load_dir.mkdir(exist_ok = False)

    initial_configuration = conf

    # maximal length of initial paths

    config0 = read_toml(toml)
    config0["runner"]["workers"] = 1
    write_toml(config0, "zero_paths.toml")
    # infretis parameters
    config = setup_config("zero_paths.toml")
    maxlen = config["simulation"]["tis_set"]["maxlength"]
    state = REPEX_state(config, minus=True)

    # setup ensembles
    state.initiate_ensembles()
    state.engines, state.engine_occ = create_engines(config)
    create_orderparameters(state.engines, config)

    # initial configuration to start from
    system0 = System()
    engine_key = list(state.engines.keys())[0]
    engine = state.engines[engine_key][0]
    engine.exe_dir = str(tmp_dir.resolve())
    if "dask" in config.keys():
        wmdrun = config["dask"]["wmdrun"][0]
    else:
        wmdrun = config["runner"]["wmdrun"][0]
    engine.set_mdrun(
        {"wmdrun": wmdrun, "exe_dir": engine.exe_dir}
    )
    system0.set_pos((os.path.abspath(initial_configuration), 0))
    system0.order = engine.calculate_order(system0)
    engine.rgen = np.random.default_rng()
    engine.modify_velocities(system0, config["simulation"]["tis_set"])

    # empty paths we will fill forwards in time in [0-] and [0+]
    path0 = Path(maxlen=maxlen)
    path1 = Path(maxlen=maxlen)

    # propagate forwards from the intiial configuration
    # note that one of these does not get integrated because
    # the initial phasepoint is either below or above interface 0
    print("Propagating in ensemble [0-]")
    status0, message0 = engine.propagate(path0, state.ensembles[0], system0)
    system0.set_pos((os.path.abspath(initial_configuration), 0))
    system0.order = path0.phasepoints[0].order
    print("Propagating in ensemble [0+]")
    status1, message1 = engine.propagate(path1, state.ensembles[1], system0)

    # we did only one integration step in ensemble 0 because
    # we started above interface 0
    if path0.length == 1:
        print("Re-propagating [0-] since we started above lambda0")
        system0.set_pos((engine.dump_config(path1.phasepoints[-1].config), 0))
        system0.order = path1.phasepoints[-1].order
        path0 = Path(maxlen=maxlen)
        status0, message0 = engine.propagate(path0, state.ensembles[0], system0)

    # or we did only one integration step in ensemble 1 because
    # we started below interface 0
    elif path1.length == 1:
        print("Re-propagating [0+] since we started below lambda0")
        system0.set_pos((engine.dump_config(path0.phasepoints[-1].config), 0))
        system0.order = path0.phasepoints[-1].order
        path1 = Path(maxlen=maxlen)
        status1, message1 = engine.propagate(path1, state.ensembles[1], system0)

    else:
        raise ValueError("Something fishy!\
                Path lengths in one of the ensembles != 1")

    # backward paths
    path0r = Path(maxlen=maxlen)
    path1r = Path(maxlen=maxlen)

    print("Propagating [0-] in reverse")
    status0, message0 = engine.propagate(
        path0r, state.ensembles[0], path0.phasepoints[0], reverse=True
    )

    print("Propagating [0+] in reverse")
    status1, message1 = engine.propagate(
        path1r, state.ensembles[1], path1.phasepoints[0], reverse=True
    )

    print("Done! Making load/ dir")
    # make load directories
    pathsf = [path0, path1]
    pathsr = [path0r, path1r]
    for i in range(2):
        dirname = load_dir / str(i)
        accepted = dirname / "accepted"
        orderfile = dirname / "order.txt"
        trajtxtfile = dirname / "traj.txt"
        dirname.mkdir()
        accepted.mkdir()
        # combine forward and backward path
        path = paste_paths(pathsr[i], pathsf[i])
        # save order paramter
        order = np.array([pp.order for pp in path.phasepoints])
        # return max op of [0+] path
        if i == 1:
            max_op = np.max(order[:,0])
        order = np.hstack((np.arange(len(order)).reshape(-1,1), np.array(order)))
        fmt = ["%d"] + ["%12.6f" for i in range(order.shape[1]-1)]
        np.savetxt(str(orderfile), order, fmt=fmt)
        N = len(order)
        # save traj.txt
        np.savetxt(
            str(trajtxtfile),
            np.c_[
                [str(i) for i in range(N)],
                [pp.config[0].split("/")[-1] for pp in path.phasepoints],
                [pp.config[1] for pp in path.phasepoints],
                [-1 if pp.vel_rev else 1 for pp in path.phasepoints],
            ],
            header=f"{'time':>10} {'trajfile':>15} {'index':>10} {'vel':>5}",
            fmt=["%10s", "%15s", "%10s", "%5s"],
        )
        # copy paths
        for trajfile in np.unique(
            [pp.config[0].split("/")[-1] for pp in path.phasepoints]
        ):
            traj_path = tmp_dir / trajfile
            traj_path.rename(accepted / trajfile)
    return max_op


def infinit(
    toml: Annotated[str, typer.Option("-toml", help="Path to .toml")] = "infretis.toml",
    log: Annotated[str, typer.Option("-log", help="File for logging output")] = "infretis_init.log",
    ):
    """The infretis initial path generator."""

    from inftools.exercises.puckering import initial_path_from_iretis

    # Based on the YouTube series:
    # https://www.youtube.com/watch?v=mW9tC2A7COs&list=PL5dSi5eZMe1iN_Uz8pTph6i8AGXhVUZIj&index=24

    # Lecture 04:
    # define the grid spacing for lambda values. All interface positions are
    # are rounded off to this vlue

    # skip this fraction of initial paths for analysis (when estimating new intf?)

    # compute efficiency measure for present set of interfaces and optimal set
    # (estimated from WHAM). If efficiency is worse than some factor, we update

    # estimated lower bound for local crossing probability

    # njumps Lp/Ls where Lp is average len full path, and Ls average len sub traj
    # lambda_cap to avoid A -> B trajs. E.g. 10% B->B paths in wf. Alternatively,
    # palce lambda_cap in half between lambdaN and lambda(N-1)

    # Lecture 05: set lambda0 lambdaN


    # TODO: restart, log file instead of print
    log = LightLogger(log)

    # we need among others parameters set in [infinit]
    config = read_toml(toml)
    # get the infinit settings from 'config' and set default parameters
    iset = set_default_infinit(config)

    if iset["cstep"]  == -1:
        log.log("Generating zero paths ...")
        init_conf = pl.Path(iset["initial_conf"]).resolve()
        max_op = generate_zero_paths(str(init_conf), toml = toml)
        log.log(f"Done with zero paths! Max op: {max_op}\n")
        iset["cstep"] = 0
        # for placing interfaces if we start with more than 1 worker
        intf = config["simulation"]["interfaces"]
        d_lambda = max_op - intf[0]
        nworkers = config["runner"]["workers"]
        lamres0 = d_lambda/nworkers
        # new interfaces to use for first infretis sim
        intf = [intf[0]] + [intf[0] + lamres0*(i+1) for i in range(nworkers-1)] + [intf[1]]
        sh_moves = ["sh", "sh"] + ["wf" for i in range(len(intf)-2)]

        # create symlink to load/1 path Nworker-1 times
        load_dir = pl.Path("load")
        load0 = load_dir / "1"
        for i in range(1,nworkers):
            loadn = load_dir / str(i + 1)
            #loadn.symlink_to(load0.relative_to(loadn.parent), target_is_directory=True)
            shutil.copytree(load0, loadn)

        c0 = read_toml(toml)
        c0["infinit"] = iset
        c0["simulation"]["interfaces"] = intf
        c0["simulation"]["shooting_moves"] = sh_moves
        write_toml(c0, "infretis.toml")

    log.log('Running infretis initialization "infinit" ...')
    if not pl.Path("infretis.toml").exists():
        print("Writing infretis.toml")
        c0 = read_toml(toml)
        c0["infinit"] = iset
        write_toml(c0, "infretis.toml")
    print_logo(step = -1)
    for iretis_steps in iset["steps_per_iter"][iset["cstep"]:]:
        log.log(f"Step {iset['cstep']}: Running infretis")
        success = run_infretis_ext(iretis_steps)
        if not success:
            print(f" *** infinit exiting loop at cstep={iset['cstep']}")
            return 1
        log.log("Updating interfaces.")
        # print(config)
        update_toml_interfaces(config)
        msg = "interfaces = ["
        msg += ", ".join([str(intf) for intf in config["simulation"]["interfaces"]])
        msg += "]"
        log.log(msg)
        #log.log("Moving and writing files.")
        #has_load = update_folders()
        iset["cstep"] += 1
        update_toml(config)
        out = initial_path_from_iretis("load", "infretis.toml", restart = "restart.toml", active_path_dir=f"load", return_pathnr = True)
        # update infretis.toml to be a restart.toml
        update_actives_toml(out)
        # rename restart file
        rename_file("restart.toml", f"restart_{iset['cstep']}.toml")
