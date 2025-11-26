from typing import Annotated
import typer


def center_periodic(
    i: Annotated[str, typer.Option("-i", help="input")],
    o: Annotated[str, typer.Option("-o", help="output")],
    c: Annotated[float, typer.Option("-c", help="The length of a cubic cell")],
    idx: Annotated[
        int, typer.Option("-idx", help="The particle idxes at center")
    ],
):
    """This command re-centers an xyz trajectory to idx."""
    from inftools.misc.xyz_help import calc_center


    atoms = None
    cnt = 0
    center = [0, 0, 0]
    with open(i) as read:
        with open(o, "w") as write:
            while True:
                # read frame header
                header = [read.readline(), read.readline()]
                if all(i == "" for i in header):
                    break
                atoms = int(header[0].rstrip())

                # read xyz
                xyzs = []
                for _ in range(atoms):
                    line = read.readline().split()
                    xyzs.append([line[0]] + [float(i) for i in line[1:]])

                # center xyzs
                if cnt == 0:
                    center = xyzs[idx][1:]
                xyzs_c = calc_center(center, xyzs, c)

                # write frame
                for line_w in header + xyzs_c:
                    write.write(line_w)
                cnt += 1
