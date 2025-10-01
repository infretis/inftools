from typing import Annotated

import numpy as np
import tomli
import tomli_w
import typer

from inftools.misc.data_helper import data_reader
from collections import defaultdict
from ase.io import read, write


def trainingset(
                data: Annotated[str, typer.Option("-data", help="data file for the respective simulation")],
                input_dir: Annotated[str, typer.Option("-input_dir", help="directory where the data files are located.")] = "load/",
                toml: Annotated[str, typer.Option("-toml", help="toml for simulation")] = "infretis.toml",
                n_frames: Annotated[int, typer.Option("-n_frames", help="Number of total frames collected from the run.")] = 100,
                out: Annotated[str, typer.Option("-out", help="name for output .xyz and _report.txt file.")] = "trainingset",
                mode: Annotated[str, typer.Option("-mode", help="mode: 'interfaces' or or 'r'")] = 'interfaces',
                skip: Annotated[list[int], typer.Option("-skip", help="skip initial lines for simulations.")] = [0]
):
    """
    Create a training set from trajectory data. Currently supports selection based on interfaces.
    Note:
    Only ASE-supported trajectory formats are supported
    (e.g., .traj, multi-frame .xyz, NetCDF)
    """
     
    print('test')
    sim = {}
    sim = {"toml": toml, "data": data}
    with open(toml, "rb") as rfile:
        sim["conf"] = tomli.load(rfile)
        sim["intf"] = sim["conf"]["simulation"]["interfaces"]
        sim["paths"] = data_reader(data)


    print(sim["intf"])

    # Group paths by the pair of interfaces between which their max_op lies
    if mode == 'interfaces':

        grouped_paths = defaultdict(list)
        for path in sim["paths"]:
            cols = path["cols"]
            for key in cols:
                grouped_paths[key].append(path)
        # print(grouped_paths[8])

        # Select round(n_frames/groups) paths from each group
        n_groups = len(grouped_paths.keys())
        n_select = int(np.ceil(n_frames / n_groups))
        # Randomly select n_select unique paths per group, ensuring no path appears in more than one group
        all_selected = set()
        selected_paths = {}
        for key in grouped_paths.keys():
            group = [p for p in grouped_paths[key] if id(p) not in all_selected]
            n_to_select = min(n_select, len(group))
            chosen = list(np.random.choice(group, size=n_to_select, replace=False))
            selected_paths[key] = chosen
            all_selected.update(id(p) for p in chosen)
        # print(selected_paths)

        # for each group, we need to read the order.txt file and select one frame
        # between the two interfaces, for group 0, it must be below interface 0.
        all_frames = {}
        for key in sorted(grouped_paths.keys()):
            upper_interface = sim["intf"][key]
            lower_interface = None
            if key > 1:
                lower_interface = sim["intf"][key - 1]
            for path in selected_paths[key]:
                order_file = f'{input_dir}/{path["pn"]}/order.txt'
                with open(order_file, "r") as f:
                    lines = f.readlines()
                    #Skip all lines with # at the beginning
                    lines = [line for line in lines if not line.startswith("#")]
                    # Select the frame for the current group
                    selected_frame = None
                    matching_frames = []
                    op = {}
                    for line in lines:
                        # second column is the order parameter, always a float 
                        cols = line.split()
                        if len(cols) > 1:
                            order_param = float(cols[1])
                            if lower_interface is None and order_param < upper_interface:
                                matching_frames.append(cols[0])
                                op[cols[0]] = order_param
                            elif lower_interface is not None:
                                if order_param > lower_interface and order_param <= upper_interface:
                                    matching_frames.append(cols[0])
                                    op[cols[0]] = order_param
                    if matching_frames:
                        selected_frame = np.random.choice(matching_frames)
                        all_frames[(path["pn"])] = (selected_frame, op[selected_frame], key)
                    # print(f"Selection group {key}: path {path['pn']} frame {selected_frame}")

        for key in sorted(all_frames.keys()):
            traj_file = f'{input_dir}/{key}/traj.txt'
            with open(traj_file, "r") as f:
                lines = f.readlines()
                # Skip all lines with # at the beginning
                lines = [line for line in lines if not line.startswith("#")]
                # Select the frame for the current group
                selected_frame = all_frames[key][0]
                for line in lines:
                    cols = line.split()
                    if cols[0] == selected_frame:
                        all_frames[key] = (cols[1], cols[2], all_frames[key][1], all_frames[key][2])  # x,y,z
                        break

        # collect all selected frames 
        frame_list = []
        report_lines = []
        for key in sorted(all_frames.keys()):
            traj_name = all_frames[key][0]
            frame_idx = all_frames[key][1]
            order_param = all_frames[key][2]
            ensemble = all_frames[key][3]
            trajectory = f'{input_dir}/{key}/accepted/{traj_name}'

            traj_path = trajectory
            # Read the trajectory file using ASE
            atoms_frame = read(traj_path, index=frame_idx)
            frame_list.append(atoms_frame)
            report_lines.append(f"{key}, {frame_idx}, {traj_name}, {order_param}, {ensemble}")

        # Write all selected frames to the output .xyz file
        write(out + ".xyz", frame_list)

        # Write the report file
        with open(out + "_report.txt", "w") as report_file:
            report_file.write("path, frame_index, trajectory_file, order_parameter, ensemble\n")
            for line in report_lines:
                report_file.write(line + "\n")
        print(f"Written {len(frame_list)} frames to {out}.xyz")
        print(f"Written report to {out}_report.txt")