import os
import numpy as np

def extract(trajfile, xcol, ycol):
    # Read and process the file
    with open(trajfile) as file:
        # Use a list comprehension to extract non-comment lines
        lines = [line for line in file if not line.startswith("#")]
        # Use another list comprehension to extract values from columns
        data = [
            (float(line.split()[xcol]), float(line.split()[ycol]))
            for line in lines
        ]
        # delete first and last timeslices
        del data[0]
        del data[-1]
        # Unzip the data into separate lists
    return data


def update_histogram(xy, factor, histogram, Minx, Miny, dx, dy):
    for timeslice in xy:
        x = timeslice[0]
        y = timeslice[1]
        ix = int((x - Minx) / dx)
        iy = int((y - Miny) / dy)
        histogram[ix, iy] += factor
    return histogram


def printhisto(xval, yval, histogram, ofile):
    with open(ofile, "w") as file:
        # Find the maximum length of any value in the histogram
        strlength = 25
        ystr = [str(y).rjust(strlength) for y in yval]
        row = " " * (strlength - 3) + r"x\y" + "|" + " ".join(ystr)
        file.write(row + "\n")
        row = "-" * len(row)
        file.write(row + "\n")

        # Write the x-values and histogram values in subsequent rows
        for i in range(histogram.shape[0]):
            x = xval[i]
            values_for_x = histogram[i, :]
            formatted_values = [
                str(val).rjust(strlength) for val in values_for_x
            ]
            row = str(x).rjust(strlength) + "|" + " ".join(formatted_values)
            file.write(row + "\n")


def calculate_free_energy(trajlabels, WFtot, Trajdir, outfolder, histo_stuff):
    print("We are now going to perform the Landau Free Energy calculations")
    print(
        "Check Free_energy.py and modify to your needs"
        + " as it contains some hard coded pytrhon script."
    )
    Nbinsx, Nbinsy = histo_stuff["nbx"], histo_stuff["nby"]
    Maxx, Minx = histo_stuff["maxx"], histo_stuff["minx"]
    Maxy, Miny = histo_stuff["maxy"], histo_stuff["miny"]
    xcol, ycol = histo_stuff["xcol"], histo_stuff["ycol"]

    histogram = np.zeros((Nbinsx, Nbinsy))
    dx = (Maxx - Minx) / Nbinsx
    dy = (Maxy - Miny) / Nbinsy
    xval = [Minx + 0.5 * dx + i * dx for i in range(Nbinsx)]
    yval = [Miny + 0.5 * dy + i * dy for i in range(Nbinsy)]
    # mid-points of the bins
    for label, factor in zip(trajlabels, WFtot):
        trajfile = Trajdir + "/" + str(label) + "/order.txt"
        xy = extract(trajfile, xcol, ycol)
        histogram = update_histogram(xy, factor, histogram, Minx, Miny, dx, dy)
    # normalize such that the highest value equals 1
    max_value = np.max(histogram)
    histogram /= max_value
    #printhisto(xval, yval, histogram, f"{outfolder}/histogram_output.txt")
    np.savetxt(os.path.join(outfolder, "histo_xval.txt"), xval)
    np.savetxt(os.path.join(outfolder, "histo_yval.txt"), yval)
    np.savetxt(os.path.join(outfolder, "histo_probability.txt"), histogram)
    histogram = -np.log(histogram)  # get Landau free energy in kBT units
    np.savetxt(os.path.join(outfolder, "histo_free_energy.txt"), histogram)
    #printhisto(xval, yval, histogram, "free_energy_output.txt")
