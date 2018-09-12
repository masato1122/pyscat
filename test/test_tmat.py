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
    
    mesh = [9, 9, 9]
    
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
    
    # --- ver.1 : auto loop
    qs, freqs, rscat, dos_imp = tmat.autoloop4tmat(ndiv_integral=200)
    
    #draw(frequencies, dos_green, dos_tmat)

if __name__ == "__main__":
    main()

