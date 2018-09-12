import numpy as np
from pyscat.tmat import Tmat
from pyscat.crystal.crys import Grid
from pyscat.crystal.vesta import mkvesta_IFCs
from pyscat.latdynam.phonon import get_grid
from pyscat.calc.dos import get_dos_green
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# --------------------------------
# input parameters
# --------------------------------
nc_scell = [[1,0,0], [0,1,0], [0,0,1]]
#primat = [[0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]]
primat = [[0, 0.25, 0.25], [0.25, 0, 0.25], [0.25, 0.25, 0]]

POS_PURE, FORCE_PURE = "POSCAR222.pure", "FORCE_pure.txt"
POS_IMP,  FORCE_IMP  = "POSCAR222.imp", "FORCE_imp.txt"

def draw(freqs, d1, d2):
    #fig, ax1 = plt.subplots()
    plt.plot(freqs, d1, '-', label="DOS(pure)")
    plt.plot(freqs, d2, '.', c='orange', label="DOS(imp)")
    plt.legend()
    #plt.show()
    
    FNAME = "dos_imp.png"
    plt.savefig(FNAME)
    print("Output:", FNAME)


def main():
    
    mesh = [3,3,3]
    
    ndiv_integral = 200
    
    tmat = Tmat( 
            impurity_site = 10,
            mesh          = mesh, 
            primat_pure   = primat,
            ncells_pure   = nc_scell,
            fposcar_pure  = POS_PURE,
            fforce_pure   = FORCE_PURE,
            fborn_pure    = "BORN",
            fposcar_imp   = POS_IMP,
            fforce_imp    = FORCE_IMP,
            fborn_imp     = "BORN" 
            )
    
    tmat.set_ifcs_and_phonon()
    tmat.check_structures()
    tmat.set_newcelles_around_impurity()
    tmat.set_ifcs4newcells()
    #tmat.visualize_IFCs_diff()
    
    # --- output     
    Nfreq = 100
    frequencies = np.linspace(0., 7., Nfreq)
    dos_green = np.zeros_like(frequencies)
    dos_tmat = np.zeros_like(frequencies)
    #for i in range(int(Nfreq/2), int(Nfreq/2)+1):
    for i, freq in enumerate(frequencies):
        
        tmat.set_green_pure(frequency=freq, ndiv_integral=200)
        if tmat.flag_calc:
            #dos_green[i] = get_dos_green(
            #        tmat.g0, len(tmat.ph_pure.get_primitive().masses))
            dos_green[i] = get_dos_green(tmat.g0)
        else:
            dos_green[i] = 0.0
        
        tmat.set_Tmatrix()
        dos_tmat[i] = tmat.get_dos_tmat() * 2.0
        print("{:10.4f} {:15.10f} {:15.10f} {:10.7f}".format(
            freq, dos_green[i], dos_tmat[i], dos_green[i]/dos_tmat[i]))
    
    draw(frequencies, dos_green, dos_tmat)

if __name__ == "__main__":
    main()

