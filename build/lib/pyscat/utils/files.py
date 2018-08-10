import numpy as np

OUTFILE = "SCATTERING_RATE"

#def initialize_output(
#        natoms, nelements, cell,
#        mesh, qpoints, f2s):
#    ofs = open(OFILE, "a")
#    ofs.write("# ----------------------------------------------\n")
#    ofs.write("# Phonon scattering rate due to an impurity\n")
#    ofs.write("# calculated by phtmat\n")
#    ofs.write("# ----------------------------------------------\n")
#    ofs.write("## GENERAL information\n")
#    ofs.write("#SYSTEM\n")
#    ofs.write("{:d} {:d}\n",format(natoms, nelements))
#    ofs.write("#END SYSTEM\n")
#    
#    ofs.write("#KPOINT\n")
#    ofs.write("{:d} {:d} {:d}\n".format(mesh[0], mesh[1], mesh[2]))
#    ofs.write("{:d}\n".format(len(irr_qpoints)))
#    
#    for iq, q in enumerate(qpoints):
#        ofs.write("{:4d}: ".format(iq))
#        for j in range(3):
#            ofs.write("{:13.7e} ".format(q[j]))
#        ofs.write("\n".format())
#    ofs.write("#END KPOINT\n")
#    
#    #ofs.write("#SMEARING\n".format())
#    #ofs.write("#END SMEARING\n".format())
#    
#    ofs.write("##END GENERAL information\n")
#    
#    # --- phonon
#    frequencies = np.sqrt(abs(f2s)) * np.sign(f2s)
#    ofs.write("##Phonon Frequency\n")
#    ofs.write("#K-point (irreducible), Branch, Omega (THz)\n")
#    for iq, q in enumerate(qpoints):
#        ofs.write("{:4d} ".format(iq))
#        for im, freq in enumerate(frequencies[iq]):
#            ofs.write("{:4d} {:4d} {:15.7f}\n".format(
#                iq+1, im+1, freq))
#    ofs.write("##END Phonon Frequency\n")
#    
#    ofs.close()


def initialize_scattering_file(mesh, nq, OFILE=OUTFILE):
    """Make file of scattering rates
    Parameters
    ------------
    mesh : array, integer, shape=(3)
        # of q-points
    nq : integer
        # of irreducible q-points
    """
    ofs = open(OFILE, "w")
    ofs.write("# ------------------------------------------\n")
    ofs.write("# Phonon scattering rate due to an impurity\n")
    ofs.write("# calculated by phtmat\n")
    ofs.write("# ------------------------------------------\n")
    ofs.write("#\n")
    ofs.write("# Units\n")
    ofs.write("# frequency, scattering rate [THz]\n")
    ofs.write("#\n")
    ofs.write("# q-mesh : ")
    for j in range(3):
        ofs.write("{:d} ".format(mesh[j]))
    ofs.write("\n")
    ofs.write("# irreducible q-points : {:d}\n".format(nq))
    ofs.write("#\n")
    ofs.close()

def dump_scattering_rate(iq, qpoint, im, f2, rscat, OFILE=OUTFILE):
    """Make file of scattering rates
    Parameters
    ------------
    qs : ndarray, float, shape=(nq, 3)
        q-points
    f2s : ndarray, float, shape=(nq, nmodes)
        Squared frequencies
    rscat : ndarray, float, shape=(nq, nmodes)
        Scattering rate
    OFILE : straing
        Output file name
    """
    freq = np.sqrt(abs(f2)) * np.sign(f2)
    ofs = open(OFILE, "a")
    if im == 0:
        if iq != 0:
            ofs.write("\n")
        ofs.write("# {:2d} qpoint : ".format(iq))
        for j in range(3):
            ofs.write("{:13.7e} ".format(qpoint[j]))
        ofs.write("\n")
    
    ofs.write("{:2d} {:13.3f} ".format(im, freq))
    if rscat is None:
        ofs.write("None\n")
    else:
        ofs.write("{:18.5e}\n".format(rscat))
    ofs.close()

def output_scattering_rates(qs, f2s, rscat, OFILE=OUTFILE):
    """Make file of scattering rates
    Parameters
    ------------
    qs : ndarray, float, shape=(nq, 3)
        q-points
    f2s : ndarray, float, shape=(nq, nmodes)
        Squared frequencies
    rscat : ndarray, float, shape=(nq, nmodes)
        Scattering rate
    OFILE : straing
        Output file name
    """
    nq = len(qs)
    nmodes = len(f2s[0])
    freqs = np.sqrt(abs(f2s)) * np.sign(f2s)
    ofs = open(OFILE, "w")
    ofs.write("# ------------------------------------------\n")
    ofs.write("# Phonon scattering rate due to an impurity\n")
    ofs.write("# calculated by phtmat\n")
    ofs.write("# ------------------------------------------\n")
    ofs.write("#\n")
    for iq, qq in enumerate(qs):
        ofs.write("# {:2d} qpoint : {:8.5f} {:8.5f} {:8.5f}\n".format(
            iq, qq[0], qq[1], qq[2]))
        for im, freq in enumerate(freqs[iq]):
            ofs.write("{:2d} {:13.3f} ".format(im, freq))
            if rscat[iq,im] is None:
                ofs.write("None\n")
            else:
                ofs.write("{:18.5e}\n".format(rscat[iq,im]))
        ofs.write("\n")
    ofs.close()
    print(" Output: ", OFILE)


