import numpy as np
from phonopy import Phonopy

from pyscat.green import Green
from pyscat.calc.dos import cal_dos_f2

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# --------------------------------
# input parameters
# --------------------------------
nc_scell = [[2,0,0], [0,2,0], [0,0,2]]
primat = [[0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]]
#primat = [[0, 0.25, 0.25], [0.25, 0, 0.25], [0.25, 0.25, 0]]

POS_PURE, FORCE_PURE = "../POSCAR", "../FORCE_SETS"

def output(FNAME, f1, d1, f2, d2):
    ofs = open(FNAME, "w")
    ofs.write("# normal DOS\n")
    ofs.write("#Freq[THz]  DOS\n")
    for iw in range(len(f1)):
        ofs.write("{:10.7f} {:10.7f}\n".format(f1[iw], d1[iw]))
    ofs.write("\n")
    #
    ofs.write("# green DOS\n")
    ofs.write("#Freq[THz]  DOS\n")
    for iw in range(len(f2)):
        ofs.write("{:10.7f} {:10.7f}\n".format(f2[iw], d2[iw]))
    ofs.write("\n")

def draw(f1, d1, f2, d2):
    fig, ax1 = plt.subplots()
    plt.xlabel("Frequency (THz)")
    plt.ylabel("DOS (a.u.)")
    ax1.plot(f1, d1, '-', label="normal")
    ax1.plot(f2, d2, '.', c='orange', label="G_0")
    ax1.legend()
    #plt.show()
    FNAME = "dos_pure.png"
    plt.savefig(FNAME)
    print("Output", FNAME)

def main():
    
    mesh = [5, 5, 5]
    
    ndiv_integral = 200
    
    green = Green(mesh, ndiv=ndiv_integral, 
            fposcar=POS_PURE, primat=primat, ncells=nc_scell,
            fforce=FORCE_PURE)
    green.set_qmesh()
    
    Nfreq = 10
    frequencies = np.linspace(0., 11., Nfreq)
    dos_green = np.zeros_like(frequencies)
    for i in range(Nfreq):
         
        green.cal_green_function(frequency=frequencies[i])
        
        dos_green[i] = green.dos_green
        print("{:10.4f} {:15.10f}".format(frequencies[i], dos_green[i]))
    
    # -- calculate DOS using IWs
    freqs = np.linspace(0., 11., 200)
    weights = np.ones(len(green.grid.qs))
    dos_normal = cal_dos_f2(green.grid.thm, weights, f2s=freqs**2)
    
    # -- draw
    output("dos.txt", freqs, dos_normal, frequencies, dos_green)
    draw(freqs, dos_normal, frequencies, dos_green)

if __name__ == "__main__":
    main()

