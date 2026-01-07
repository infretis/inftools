import os
import numpy as np
import matplotlib.pyplot as plt

def extract(trajfile, xcol, ycol=None):
    # Read and process the file
    traj = np.loadtxt(trajfile)
    data = traj[1:-1, xcol] # remove first and last frames
    if ycol is not None:
        data = np.vstack((data, traj[1:-1, ycol]))
    return data


def update_histogram(data, factor, histogram, Minx, Maxx, Miny, dx, dy):
    if Miny is not None and dy is not None:
        x = data[0]
        y = data[1]

        ix = ((x - Minx) / dx).astype(int)
        iy = ((y - Miny) / dy).astype(int)

        np.add.at(histogram, (ix, iy), factor)

    else:
        x = data if data.ndim == 1 else data[:,0] # make sure x is one dimensional
        mask = (x >= Minx) & (x <= Maxx)
        ix = ((x[mask] - Minx) / dx).astype(int)
        ix[ix == len(histogram)] = len(histogram) - 1
        np.add.at(histogram, ix, factor)

    return histogram


def calculate_free_energy(trajlabels, WFtot, Trajdir, outfolder, histo_stuff, lm1, lA, lB, xi, sym):
    print("We are now going to perform the Landau Free Energy calculations")
    Nbinsx, Nbinsy = histo_stuff["nbx"], histo_stuff["nby"]
    Maxx, Minx = histo_stuff["maxx"], histo_stuff["minx"]
    Maxy, Miny = histo_stuff["maxy"], histo_stuff["miny"]
    xcol, ycol = histo_stuff["xcol"], histo_stuff["ycol"]
    
    # Set these optiongs when histograms will be stitched together later
    setnbins = histo_stuff["setnbins"] 
    Minbinx, Maxbinx = histo_stuff["minbx"], histo_stuff["maxbx"]
    Binw = histo_stuff["binw"]

    # if lm1 is being used the Minx and Maxx will be set to lm1 and lB.
    if lm1 is not None:
        Minx = lm1
        Maxx = lB

    if any(var is None for var in [Nbinsy, Maxy, Miny, ycol]):
        none_vars = [name for name, var in zip(["nby", "maxy", "miny", "ycol"], [Nbinsy, Maxy, Miny, ycol]) if var is None]
        assert all(var is None for var in [Nbinsy, Maxy, Miny, ycol]), \
            f"The following variables are None and should be set: {', '.join(none_vars)}"
    if Nbinsy is not None:
        histogram = np.zeros((Nbinsx, Nbinsy))
        dy = (Maxy - Miny) / Nbinsy
        yval = Miny + np.arange(Nbinsy)*dy + 0.5 * dy

    else:
        histogram = np.zeros(Nbinsx)
        dy = None
        yval = None
    
    if setnbins:
        # Bins are determined with the given binwidth (Binw) and the given range (Minbinx, Maxbinx)
        # Very useful when stitching multiple histograms together, e.g. MSM of membrane.
        dx = Binw
        Nbinsx = (Maxbinx - Minbinx) / dx
        histogram =np.zeros(int(Nbinsx))
        unw_histogram = np.zeros(int(Nbinsx))
        Minx = Minbinx
        binedges = Minx + np.arange(Nbinsx+1)*dx
        xval = Minx + np.arange(Nbinsx)*dx + 0.5 * dx
    
    else:
        dx = (Maxx - Minx) / Nbinsx
        binedges = Minx + np.arange(Nbinsx+1)*dx
        xval = Minx + np.arange(Nbinsx)*dx + 0.5 * dx

    for label, factor in zip(trajlabels, WFtot):
        # Collect trajectory date and create histogram
        trajfile = Trajdir + "/" + str(label) + "/order.txt"
        data = extract(trajfile, xcol, ycol)
        histogram = update_histogram(data, factor, histogram, Minx, Maxx, Miny, dx, dy)

    index_lA = None
    if xi is not None:
        # Corrects histogram untill lambda A with xi 
        index_lA = np.where(np.array(xval) > lA)[0][0]
        histogram[:index_lA] = histogram[:index_lA] / xi

    if not yval is None:
        np.savetxt(os.path.join(outfolder, "histo_yval.txt"), yval)
        
    np.savetxt(os.path.join(outfolder, "histo_probability.txt"), np.column_stack((xval, histogram)), 
                header=f"lm1={lm1}, lA={lA}, lB={lB}, first_bin_after_lA_index={index_lA}", comments="# ")
    
    np.savetxt(os.path.join(outfolder, "histo_binedges.txt"), (binedges), 
                header=f"lm1={lm1}, lA={lA}, lB={lB}, first_bin_after_lA_index={index_lA}", comments="# ")  

    np.savetxt(os.path.join(outfolder, "unweighted_histo_probability.txt"), np.column_stack((xval, unw_histogram)), 
                header=f"lm1={lm1}, lA={lA}, lB={lB}, first_bin_after_lA_index={index_lA}", comments="# ") 
    
    # normalize such that the highest value equals 1
    max_value = np.max(histogram)    
    conditional_free_energy = -np.log(histogram/max_value)  # get Landau free energy in kBT units
    np.savetxt(os.path.join(outfolder, "histo_free_energy.txt"), np.column_stack((xval, conditional_free_energy)),
                header=f"lm1={lm1}, lA={lA}, lB={lB}, first_bin_after_lA_index={index_lA}", comments="# ") 

    if Nbinsy is None:
        plt.plot(xval, histogram, marker='o', markersize=3, linewidth=1)
        plt.xlabel("Order parameter (Å)")
        plt.ylabel("Probability")
        plt.grid()
        plt.tight_layout()
        plt.savefig(os.path.join(outfolder, "histogram.png"), dpi=300)
        plt.close()

        plt.plot(xval, conditional_free_energy, marker='o', markersize=3, linewidth=1)
        plt.xlabel("Order parameter (Å)")
        plt.ylabel("Free energy (kBT)")
        plt.grid()
        plt.tight_layout()
        plt.savefig(os.path.join(outfolder, "free_energy_conditional.png"), dpi=300)
        plt.close()

    if Nbinsy is None and sym:
        # If sym true, calculates and plots the backwards conditional and symmetrized free energy
        histogram = np.append(histogram, np.zeros(index_lA, dtype=histogram.dtype))        
        histo_sym = histogram + histogram[::-1]
        max_value = np.max(histo_sym)
        free_energy_sym = -np.log(histo_sym/max_value)
        xval_sym = [Minx + 0.5 * dx + i * dx for i in range(Nbinsx+index_lA)]
        plt.plot(xval, conditional_free_energy, label="Conditional A→B", marker='o', markersize=3, linewidth=1)
        plt.plot(xval_sym[index_lA:], conditional_free_energy[::-1], label="Conditional B→A", marker='o', markersize=3, linewidth=1)
        plt.plot(xval_sym, free_energy_sym, label="Symmetrized", marker='o', markersize=3, linewidth=1)
        np.savetxt(os.path.join(outfolder, "histo_free_energy_sym.txt"), np.column_stack((xval_sym, free_energy_sym)), 
                   header=f"lm1={lm1}, lA={lA}, lB={lB}, first_bin_after_lA_index={index_lA}", comments="# ")        
        plt.xlabel("Order parameter (Å)")
        plt.ylabel("Free energy (kBT)")
        plt.grid()
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(outfolder, "free_energy_symmetrized.png"), dpi=300)
        plt.close()

