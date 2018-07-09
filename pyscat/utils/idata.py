import numpy as np
from phonopy import Phonopy

def _check_true_none(variable, label):
    if variable is None:
        print("Error: \"{:s}\" is not given.".format(label))
        import sys
        sys.exit()

def _check_file(file):
    import os
    if os.path.exists(file) == False:
        print("Cannot find", file)
        import sys
        sys.exit()

class Idata():
    def __init__( self,
                  mesh     = None,
                  fposcar  = None,
                  fforce   = None,
                  fborn    = None,
                  primat = None,
                  ncells   = None):
        """
        Parameters:
            impurity_site: atom site for the impurity
                           Currently, only one substitution is supported.
            fposcar     : POSCAR file name (unit cell)
            fforce      : file of force & displacement sets
            fborn       : File of Born effective charge
            nmesh       : # of meshes
            primat      : primivice vectors scaled by unitcell (POSCAR file)
            ncells      : # of unit cells in the structure used to make the force sets
        
            nat_prim    : # of atoms in a primitive cell
            nat_super   : # of atoms in a supercell
            nmod        : # of modes

            scID        : order of atom indices
        """
        self.nmesh    = mesh
        self.fposcar  = fposcar
        self.fforce   = fforce
        self.fborn    = fborn
        self.primat   = primat
        self.ncells   = ncells

        self.nat_prim  = None
        self.nat_super = None
        self.nmod      = None

        self.scID      = None
        
        _check_true_none( self.primat, "primat" )
        _check_true_none( ncells, "ncells" )
        _check_file( fposcar )
        _check_file( fforce )
        
        if fborn is not None:
            _check_file( fborn )

    def set_natoms( self, phonon ):
        self.nat_prim  = len( phonon.get_primitive().masses )
        self.nat_super = len( phonon.get_supercell().masses )
        self.nmod      = 3 * self.nat_prim
        


