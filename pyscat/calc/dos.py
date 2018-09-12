import numpy as np

def get_dos_iws(iws):
    """Calculate DOS with IWs.
    Parameter
    ----------
    iws : ndarray, float, shape=(((nqpoints, nfreq, nmodes)))
    nqpoints : integer
        # of q-points
    nfreq : integer
        # of frequencies
    nmodes : integer
        # of modes

    Return
    ----------
    dos : array, float, shape=(nfreq,)

    """
    nfreq = len(iws[0])
    dos = np.zeros(nfreq)
    for ifreq in range(nfreq):
        dos[ifreq] = np.sum(iws[:,ifreq,:])
    return dos

def cal_dos_f2(thm, weights, freqs=None, f2s=None):
    """Calculator of DOS from integration weights
    Parameters
    ---------------
    freqs : array, float, shape=(nfreq)
    thm : TetrahedronMethod
        See descriptions in Phonopy in detail
    weights : array, float, shape=(nqpoints,)
        Weight of each q-points
        Type is float, but in fact, integer.
    
    nfreq : integer
        = len(freqs)
    nqpoints : integer
        # of q-points
    
    Return
    ---------
    dos : array, float, shape=(ndiv,)
        Calculated DOS
    """
    if f2s is None:
        if freqs is not None:
            f2s = freqs**2
        else:
            print("Error: input a list of frequencies.")
            import sys
            sys.exit()
    thm.set(value='I', frequency_points=f2s)
    dos = np.zeros_like(f2s)
    for iq, iw in enumerate(thm):
        dos += np.sum(iw * weights[iq], axis=1)
    
    return dos
    
def get_dos_green(g0, multiplicity=None):
    """Calculate phonon DOS using Green's function at a given frequency
    Parameters
    ------------
    g0 : ndarray, complex, shape=((3*natoms, 3*natoms))
        Green's function of the pure crystal
    #nat_prim : integer
    #    # of atoms in the primitive cell
    """
    if multiplicity is None:
        multi_long = np.ones(len(g0))
    else:
        if 3*len(multiplicity) != len(g0):
            print("Error {:d} != {:d}".format(
                    3*len(multiplicity),
                    len(g0)))
            import sys
            sys.exit()
        
        multi_long = np.zeros(len(g0))
        for i in range(len(multiplicity)):
            multi_long[3*i:3*(i+1)] = multiplicity[i] * np.ones(3)
    dos = (np.sum(np.diag(np.imag(g0))/multi_long) / np.pi)
    return dos


