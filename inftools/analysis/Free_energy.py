import os
import numpy as np
import matplotlib.pyplot as plt
import tarfile
import io


def build_member_map(tarfile_path, suffix="order.txt"):
    """
    Returns dict: 'label/order.txt' -> TarInfo
    """
    with tarfile.open(tarfile_path, "r") as tf:
        members = [
            m for m in tf.getmembers()
            if m.isfile() and m.name.endswith(f"/{suffix}")
        ]

    member_map = {}
    for m in members:
        rel = "load/" + "/".join(m.name.split("/")[-2:])  # label/order.txt
        member_map[rel] = m

    return member_map


def extract(trajfile, xcol, ycol=None, member_map={}):
    # Read and process the file
    if not member_map:
        traj = np.loadtxt(trajfile)
    else:
        m = member_map[trajfile]
        with tarfile.open("infretis_data.tar", "r") as tf:
            with tf.extractfile(m) as fh:
                traj = np.loadtxt(io.TextIOWrapper(fh))

    data = traj[1:-1, xcol] # remove first and last frames
    if ycol is not None:
        data = np.vstack((data, traj[1:-1, ycol]))
    return data


def update_histogram(data, factor, histogram, Minx, Miny, dx, dy):
    if Miny is not None and dy is not None:
        x = data[0]
        y = data[1]

        ix = ((x - Minx) / dx).astype(int)
        iy = ((y - Miny) / dy).astype(int)

        np.add.at(histogram, (ix, iy), factor)

    else:
        x = data if data.ndim == 1 else data[:,0] # make sure x is one dimensional
        ix = ((x - Minx) / dx).astype(int)
        np.add.at(histogram, ix, factor)

    return histogram



def calculate_free_energy(trajlabels, WFtot, Trajdir, outfolder, histo_stuff, interfaces):
    print("We are now going to perform the Landau Free Energy calculations")
    Nbinsx, Nbinsy = histo_stuff["nbx"], histo_stuff["nby"]
    Maxx, Minx = histo_stuff["maxx"], histo_stuff["minx"]
    Maxy, Miny = histo_stuff["maxy"], histo_stuff["miny"]
    xcol, ycol = histo_stuff["xcol"], histo_stuff["ycol"]

    if any(var is None for var in [Nbinsy, Maxy, Miny, ycol]):
        none_vars = [name for name, var in zip(["nby", "maxy", "miny", "ycol"], [Nbinsy, Maxy, Miny, ycol]) if var is None]
        assert all(var is None for var in [Nbinsy, Maxy, Miny, ycol]), \
            f"The following variables are None and should be set: {', '.join(none_vars)}"
    if Nbinsy is not None:
        histogram = np.zeros((Nbinsx, Nbinsy))
        histogram_comm = np.zeros((Nbinsx, Nbinsy))
        dy = (Maxy - Miny) / Nbinsy
        yval = [Miny + 0.5 * dy + i * dy for i in range(Nbinsy)]
    else:
        histogram = np.zeros(Nbinsx)
        histogram_comm = np.zeros(Nbinsx)
        dy = None
        yval = None
    dx = (Maxx - Minx) / Nbinsx
    xval = [Minx + 0.5 * dx + i * dx for i in range(Nbinsx)]

    member_map = {}
    tarfile_path = "infretis_data.tar"
    if os.path.exists(tarfile_path):
        member_map = build_member_map(tarfile_path)

    for label, factor in zip(trajlabels, WFtot):
        trajfile = Trajdir + "/" + str(label) + "/order.txt"
        data = extract(trajfile, xcol, ycol, member_map)
        data0 = np.loadtxt(trajfile)
        histogram = update_histogram(data, factor, histogram, Minx, Miny, dx, dy)

        # calculate committor, here each path is determined to be reactive
        # or not based on the initial interfaces and order.txt files
        # 2 for reactive, 1 for unreactive for now.
        reactive = 2 if data0[-1, 1] > interfaces[-1] else 1
        histogram_comm = update_histogram(data, factor*reactive, histogram_comm, Minx, Miny, dx, dy)

    # calculate the committor projection
    histogram_comm[histogram_comm != 0] /= histogram[histogram_comm != 0]

    # set non-assigned values (0) to np.inf and shift actual values between [0, 1]
    histogram_comm[histogram_comm == 0] = np.inf
    histogram_comm -= 1

    # normalize such that the highest value equals 1
    max_value = np.max(histogram)
    histogram /= max_value

    np.savetxt(os.path.join(outfolder, "histo_xval.txt"), xval)
    if not yval is None:
        np.savetxt(os.path.join(outfolder, "histo_yval.txt"), yval)
    np.savetxt(os.path.join(outfolder, "histo_probability.txt"), histogram)
    np.savetxt(os.path.join(outfolder, "histo_committor.txt"), histogram_comm)

    histogram = -np.log(histogram)  # get Landau free energy in kBT units
    np.savetxt(os.path.join(outfolder, "histo_free_energy.txt"), histogram)
