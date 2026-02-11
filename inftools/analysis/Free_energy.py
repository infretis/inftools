import os
import numpy as np
import tqdm
import h5py

def extract(trajfile, xcol, ycol=None):
    # Read and process the file
    traj = np.loadtxt(trajfile)
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


def calculate_free_energy(trajlabels, WFtot, Trajdir, outfolder, histo_stuff):
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
        dy = (Maxy - Miny) / Nbinsy
        yval = [Miny + 0.5 * dy + i * dy for i in range(Nbinsy)]
    else:
        histogram = np.zeros(Nbinsx)
        dy = None
        yval = None
    dx = (Maxx - Minx) / Nbinsx
    xval = [Minx + 0.5 * dx + i * dx for i in range(Nbinsx)]

    if os.path.exists("order.h5"):
        print("Using order.h5 for whamming")
        with h5py.File("order.h5", "r") as h5f:
            for label, factor in tqdm.tqdm(zip(trajlabels, WFtot)):
                data = h5f[f"arrays/{label:06d}"]
                if ycol is not None:
                    data = data[1:-1,[xcol,ycol]].T
                else:
                    data = data[1:-1, xcol]

                histogram = update_histogram(data, factor, histogram, Minx, Miny, dx, dy)
    else:
        for label, factor in zip(trajlabels, WFtot):
            trajfile = Trajdir + "/" + str(label) + "/order.txt"
            data = extract(trajfile, xcol, ycol)
            histogram = update_histogram(data, factor, histogram, Minx, Miny, dx, dy)


    # normalize such that the highest value equals 1
    max_value = np.max(histogram)
    histogram /= max_value

    np.savetxt(os.path.join(outfolder, "histo_xval.txt"), xval)
    if not yval is None:
        np.savetxt(os.path.join(outfolder, "histo_yval.txt"), yval)
    np.savetxt(os.path.join(outfolder, "histo_probability.txt"), histogram)
    
    histogram = -np.log(histogram)  # get Landau free energy in kBT units
    np.savetxt(os.path.join(outfolder, "histo_free_energy.txt"), histogram)
