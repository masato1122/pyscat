import numpy as np
from pyscat.tmat import Tmat
from pyscat.crystal.crys import Grid
from pyscat.crystal.vesta import mkvesta_IFCs
from pyscat.latdynam.phonon import get_grid
import matplotlib.pyplot as plt

# --------------------------------
# input parameters
# --------------------------------
nc_scell = [[1,0,0], [0,1,0], [0,0,1]]
#primat = [[0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]]
primat = [[0, 0.25, 0.25], [0.25, 0, 0.25], [0.25, 0.25, 0]]

POS_PURE, FORCE_PURE = "POSCAR222.pure", "FORCE_pure.txt"
POS_IMP,  FORCE_IMP  = "POSCAR222.imp", "FORCE_imp.txt"

def draw(f1, d1, f2, d2):
    fig, ax1 = plt.subplots()
    ax1.plot(f1, d1, '-', label="normal")
    ax1.plot(f2, d2, '.', c='orange', label="G_0")
    ax1.legend()
    plt.show()

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
    
    # --- ver.1 : auto loop
    #tmat.autoloop4tmat(ndiv_integral=200)
    
    # --- ver.2 : manual
    qs, weights, f2s, evecs = tmat.get_grid4summation() 
    nq = len(qs)
    nmodes = len(f2s)
    dos_tmat = np.zeros_like(qs)
    for iq, qpoint in enumerate(qs):
        for im, f2 in enumerate(f2s[iq]):
               
            # --- calculate Green's function
            tmat.set_green_pure(f2=f2, ndiv_integral=200)
            
            #nat = tmat.nc_pure.nat
            #mkvesta_IFCs("green.vesta",
            #            np.imag(tmat.g0).reshape(nat, nat, 3, 3),
            #            tmat.nc_pure.scaled_positions,
            #            tmat.nc_pure.imp_site, 0, 1,
            #            TYPE="log")
            
            # --- calculate T-matrix
            tmat.set_Tmatrix()
            dos_tmat[iq,im] = tmat.get_dos_tmat()
            rscat = (
                    tmat.get_scattering_rate(
                        qpoint=qpoint,
                        f2=f2,
                        evec=evecs[iq,:,im],
                        iq=iq, imode=im
                        )
                    )
     
    
    
    exit()
    
    Nfreq = 5
    frequencies = np.linspace(0., 7., Nfreq)
    dos_green = np.zeros_like(frequencies)
    #for i in range(int(Nfreq/2), int(Nfreq/2)+1):
    for i in range(Nfreq):
        
        green.cal_green_function( frequency=frequencies[i] )
        
        dos_green[i] = green.dos_green
        print("{:10.4f} {:15.10f}".format(frequencies[i], green.dos_green))
    
    
    draw( np.sqrt(np.sqrt( green.f2s_integral**2)) * np.sign(green.f2s_integral) , \
            green.dos_normal, \
            frequencies, dos_green)

if __name__ == "__main__":
    main()

