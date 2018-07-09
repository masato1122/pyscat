import numpy as np
from phonopy import Phonopy
from phtmat.tmat import Tmat
from phtmat.crystal.crys import Grid
from phtmat.crystal.vesta import mkvesta_IFCs
from phtmat.latdynam.phonon import get_grid
import matplotlib.pyplot as plt

# --------------------------------
# input parameters
# --------------------------------
nc_scell = [[1,0,0], [0,1,0], [0,0,1]]
primat = [[0, 0.25, 0.25], [0.25, 0, 0.25], [0.25, 0.25, 0]]

POS_PURE, FORCE_PURE = "SPOSCAR_pure", "FORCE_SETS_pure"
POS_IMP,  FORCE_IMP  = "SPOSCAR_imp", "FORCE_SETS_imp"

def draw(freqs, rscat):

    nq = len(freqs)
    nmodes = len(freqs[0])

    fig, ax1 = plt.subplots()
    for iq in range(nq):
        ax1.plot(freqs[iq,:], rscat[:], '.', label="normal")
    ax1.legend()
    plt.show()

def main():
    
    mesh = [5, 5, 5]
    
    ndiv_integral = 200
     
    tmat = Tmat( 
            impurity_site = 0,
            mesh          = mesh, 
            primat_pure   = primat,
            ncells_pure   = nc_scell,
            fposcar_pure  = POS_PURE,
            fforce_pure   = FORCE_PURE,
            fposcar_imp   = POS_IMP,
            fforce_imp    = FORCE_IMP,
            )
    
    tmat.set_ifcs_and_phonon()
    tmat.check_structures(cell=False)
    tmat.set_newcelles_around_impurity()
    tmat.set_ifcs4newcells()
    
    #mkvesta_IFCs("ifcs_pure_00.vesta",
    #        tmat.nc_pure.IFCs,
    #        tmat.nc_pure.scaled_positions,
    #        tmat.nc_pure.imp_site, 0, 0,
    #        cell=tmat.nc_pure.cell,
    #        IFC_min=-0.7, IFC_max=0.7)
    #mkvesta_IFCs("ifcs_imp_00.vesta",
    #        tmat.nc_imp.IFCs,
    #        tmat.nc_imp.scaled_positions,
    #        tmat.nc_imp.imp_site, 0, 0,
    #        cell=tmat.nc_imp.cell,
    #        IFC_min=-0.7, IFC_max=0.7)
    #mkvesta_IFCs("ifcs_diff_00.vesta",
    #        tmat.nc_imp.IFCs - tmat.nc_pure.IFCs,
    #        tmat.nc_pure.scaled_positions,
    #        tmat.nc_pure.imp_site, 0, 0,
    #        cell=tmat.nc_pure.cell,
    #        IFC_min=-0.15, IFC_max=0.15)

    # --- ver.1 : auto loop
    qs, freqs, rscat = tmat.autoloop4tmat(ndiv_integral=200)

    # --- ver.2 : manual
    #qs, weights, f2s, evecs = tmat.get_grid4summation() 
    #nq = len(qs)
    #nmodes = len(f2s)
    #for iq, qpoint in enumerate(qs):
    #    for im, f2 in enumerate(f2s[iq]):
    #           
    #        # --- calculate Green's function
    #        tmat.set_green_pure(f2=f2, ndiv_integral=200)
    #        
    #        #nat = tmat.nc_pure.nat
    #        #mkvesta_IFCs("green.vesta",
    #        #            np.imag(tmat.g0).reshape(nat, nat, 3, 3),
    #        #            tmat.nc_pure.scaled_positions,
    #        #            tmat.nc_pure.imp_site, 0, 1,
    #        #            TYPE="log")
    #        
    #        # --- calculate T-matrix
    #        tmat.set_tmatrix()
    #         
    #        rscat, _ = (
    #                tmat.get_scattering_rate(
    #                    qpoint=qpoint,
    #                    f2=f2,
    #                    evec=evecs[iq,:,im],
    #                    iq=iq, imode=im
    #                    )
    #                )
    # 
    
     
    draw(freqs, rscat)
    
if __name__ == "__main__":
    main()

