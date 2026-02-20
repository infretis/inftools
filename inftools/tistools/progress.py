from typer import Option as Opt
from typing import Annotated as And


def progress(
    toml: And[str,  Opt("-toml", help="The infretis.toml file")] = "infretis.toml",
    restart: And[str,  Opt("-restart", help="The restart toml file")] = "restart.toml",
    sim: And[str,  Opt("-sim", help="The sim.log file")] = "sim.log",
    poll: And[float,  Opt("-poll", help="The poll interval in seconds")] = 1.0,

    ):
    """Plots a progress bar in the terminal for an ongoing (or completed) infretis simulation
    this function depends on a python package (rich).
    """

    import time
    import importlib
    import os
    import re
    from inftools.misc.infinit_helper import read_toml

    if not os.path.exists(toml):
        print(f"{toml} does not exist!")
        return

    if not os.path.exists(restart):
        print(f"{restart} does not exist!")
        return

    if not os.path.exists(sim):
        print(f"{sim} does not exist!")
        return

    if importlib.util.find_spec("rich") is None:
        print("rich is not installed")
        return
    else:
        from rich.progress import Progress, BarColumn, TextColumn, TimeRemainingColumn, TimeElapsedColumn
        from rich.live import Live
        from rich.spinner import Spinner
        from rich.text import Text
        from rich.table import Table
        from rich.console import Group

    SHOOTING_RE = re.compile(
        r"\[INFO\]: shooting (sh|wf) in ensembles: ([\d ]+?) with paths: ([\d ]+?) and worker: (\d+)"
    )

    def parse_sim_log_tail(filepath, n_lines=500):
        worker_state = {}
        try:
            with open(filepath, "r") as f:
                lines = f.readlines()[-n_lines:]
            for line in lines:
                m = SHOOTING_RE.search(line)
                if m:
                    shoot_type, ensembles, paths, worker = m.groups()
                    worker_state[int(worker)] = {
                        "type": shoot_type,
                        "ensembles": ensembles.strip(),
                        "paths": paths.strip(),
                    }
        except FileNotFoundError:
            pass
        return worker_state

    def build_layout(progress, worker_state, num_workers, spinners):
        renderables = [progress]
        for w in range(num_workers):
            state = worker_state.get(w)
            row = Table.grid(padding=(0, 1))
            if state:
                text = Text()
                text.append(f"worker {w}:  ", style="bold cyan")
                text.append("move: ", style="dim")
                text.append(f"({state['type']}) ", style="yellow")
                text.append("ensembles: ", style="dim")
                text.append(state["ensembles"], style="green")
                text.append("  paths: ", style="dim")
                text.append(state["paths"], style="green")
                row.add_row(spinners[w], text)
            else:
                text = Text()
                text.append(f"worker {w}  ", style="bold cyan")
                text.append("idle", style="dim")
                row.add_row(Text("  "), text)
            renderables.append(row)
        return Group(*renderables)

    config = read_toml(toml)
    tsteps = config["simulation"]["steps"]
    num_workers = config["runner"]["workers"]

    last_mtime = None
    cstep = 0
    spinners = [Spinner("dots", style="green") for _ in range(num_workers)]

    progress = Progress(
        TextColumn("[bold blue]{task.description}"),
        BarColumn(),
        TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
        TextColumn("Step {task.completed}/{task.total}"),
        TimeElapsedColumn(),
        TimeRemainingColumn(),
    )
    task = progress.add_task("infretis running...", total=tsteps)

    with Live(refresh_per_second=4) as live:
        while cstep < tsteps:
            try:
                mtime = os.path.getmtime(restart)
                if mtime != last_mtime:
                    last_mtime = mtime
                    # sometimes we are reading incomplete toml files
                    try:
                        cstep = read_toml(restart)["current"]["cstep"]
                    except:
                        print("boing")
                        pass
                    progress.update(task, completed=cstep)
            except FileNotFoundError:
                pass

            worker_state = parse_sim_log_tail(sim)
            live.update(build_layout(progress, worker_state, num_workers, spinners))
            time.sleep(poll)
