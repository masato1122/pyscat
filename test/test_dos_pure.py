import numpy as np
from pyscat.green import Green
from pyscat.tmat import Tmat
from pyscat.crystal.crys import Grid
from pyscat.crystal.vesta import mkvesta_IFCs
from pyscat.latdynam.phonon import get_grid
from pyscat.calc.dos import cal_dos_f2

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

def draw(f0, d0, f1, d1, d2):
    #fig, ax1 = plt.subplots()
    plt.plot(f0, d0, '-', c='black', label="DOS(normal)")
    plt.plot(f1, d1, '.', c='blue', label="DOS(Green)")
    plt.plot(f1, d2, '_', c='orange', label="DOS(Tmat)")
    plt.legend()
    #plt.show()
    
    FNAME = "dos_pure.png"
    plt.savefig(FNAME)
    print("Output:", FNAME)

def main():
    
    mesh = [5,5,5]
    ndiv_integral = 200
    
    #--- use Green class
    green = Green(mesh, ndiv=ndiv_integral,
            fposcar=POS_PURE, primat=primat, ncells=nc_scell,
            fforce=FORCE_PURE)
    green.set_qmesh()

    #--- use Tmat class
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
    
    Nfreq = 3
    frequencies = np.linspace(2., 6., Nfreq)
    dos_G = np.zeros_like(frequencies)
    dos_T = np.zeros_like(frequencies)
    for i, freq in enumerate(frequencies):
        
        green.cal_green_function(frequency=freq)
        dos_G[i] = green.dos_green
        
        tmat.set_green_pure(frequency=freq, ndiv_integral=200)
        dos_T[i] = tmat.get_dos_pure() * 3*2
        print("{:10.4f} {:15.10f} {:15.10f}".format(freq, dos_G[i], dos_T[i]))
    
    #--- normal method
    freqs = np.linspace(0., 7., 200)
    weights = np.ones(len(green.grid.qs))
    dos_normal = cal_dos_f2(green.grid.thm, weights, f2s=freqs**2)
        
    #nat_prim = len(tmat.ph_pure.get_primitive().masses)
    
    dos_T *= float(len(tmat.g0)/3)
    draw(freqs, dos_normal, frequencies, dos_G, dos_T)

if __name__ == "__main__":
    main()

