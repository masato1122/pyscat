# -*- coding: utf-8 -*-
import numpy as np
from phonopy import Phonopy
from phonopy.structure.cells import get_smallest_vectors

from idata import Idata
from green import Green
from phonon import get_phonon
from utils import set_frequency2
import check

"""
  
"""
class Tmat(object):
    def __init__(self,
                 impurity_site=None,
                 mesh_pure=None, primat_pure=None, ncells_pure=None,
                 fposcar_pure = "POSCAR_pure", 
                 fforce_pure  = "FORCE_SET_pure",
                 fborn_pure   = None,
                 mesh_imp     = None, 
                 primat_imp   = np.eye(3),
                 ncells_imp   = np.eye(3, dtype=int),
                 fposcar_imp  = "POSCAR_imp",
                 fforce_imp   = "FORCE_SET_imp",
                 fborn_imp    = None ):
        """
        Parameters:
            impurity_site: atom site for the impurity 
                           Currently, only one substitution is supported.
            ***_pure    : for pure crystal
            ***_imp     : for the structure with an impurity
            fposcar_*, fforce_*, fborn_*
                        : files for POSCAR, force sets, and Born effective charge
            mesh_*, primat_*, ncells_*
                        : # of meshes, primitive matrix, and # of cells
            
            See also the description of Idata class.

        Variables:
            ph_***          : phonon information
            _nat_prim_***   : # of atoms in the primitive cell
            _nat_super_***  : # of atoms in the super cell
            _nmod_***       : # of modes
        
            _flag_cell      : result of the check of cell size 
            _flag_natom     : result of the check of number of atoms
            _flag_positions : result of the check of atomic positions
          
        Note:
            The order of atoms in the two POSCAR file and the structure used for
            the force sets should correspond to each other.
        """
        # ---------------------------
        # input variables
        # ---------------------------
        # -- for pure crystal
        self.idat_pure = Idata( mesh     = mesh_pure,
                                fposcar  = fposcar_pure,
                                fforce   = fforce_pure,
                                fborn    = fborn_pure,
                                primat   = primat_pure,
                                ncells   = ncells_pure)
        
        # -- for the structure w/ an impurity
        if mesh_imp is None:
            mesh_imp = mesh_pure
        
        self.idat_imp = Idata( mesh     = mesh_imp,
                               fposcar  = fposcar_imp,
                               fforce   = fforce_imp,
                               fborn    = fborn_imp,
                               primat   = primat_imp,
                               ncells   = ncells_imp)
        
        # -- impurity site
        if impurity_site is None:
            print("Error: impurity site is not given.")
            import sys
            sys.exit()
        self.imp_site = impurity_site
        
        # --------------------------------------------------
        # Get "Phonopy", force constants, masses, etc. and
        # set q-meshes
        # --------------------------------------------------
        # -- for the pure crystal
        self.ph_pure = get_phonon( self.idat_pure )
        
        self.ph_pure.set_mesh( self.idat_pure.nmesh, is_mesh_symmetry=False, is_eigenvectors=True)
        self._nat_prim_pure  = len( self.ph_pure.get_primitive().masses )
        self._nat_super_pure = len( self.ph_pure.get_supercell().masses )
        self._nmod_pure      = self._nat_prim_pure * 3
        
        # -- for the structure w/ an impurity
        self.ph_imp = get_phonon( self.idat_imp )
        
        self.ph_imp.set_mesh( self.idat_imp.nmesh, is_mesh_symmetry=False, is_eigenvectors=True)
        self._nat_prim_imp  = len( self.ph_imp.get_primitive().masses )
        self._nat_super_imp = len( self.ph_imp.get_supercell().masses )
        self._nmod_imp      = self._nat_prim_imp * 3
        
        # -----------------------------------
        # variables that will be used
        # -----------------------------------
        self.f2_target = None
        
        self._flag_cell      = None
        self._flag_natom     = None
        self._flag_positions = None
        
     
    def check_structures( self ):
        """ Check the structures.
        Note:
            What is checked may be changed depending on what kind of impurities 
            can be calculated with this code.
        """
        sc_pure = self.ph_pure.get_supercell()
        sc_imp  = self.ph_imp.get_supercell()
        
        self._flag_cell  = check.check_cells( sc_pure.cell, sc_imp.cell, tolerance=1e-5)
        
        self._flag_natom = check.check_number_of_atoms( self._nat_super_pure, self._nat_super_imp )
        
        self._flag_positions = check.check_positions(
                np.dot( sc_pure.scaled_positions, sc_pure.cell ),
                np.dot( sc_imp.scaled_positions, sc_imp.cell),
                tolerance=1.0 )
        
        print(" Structures are compatible.")
        exit()
    
    #def set_atoms_around_impurity( self ):

    def cal_tmatrix( self, 
                     frequency=None, f2=None, 
                     ph_pure=None, ph_imp=None):
        
        self.f2_target = set_frequency2( frequency=frequency, f2=f2 )
         
        # -----------------------------------------------------------
        # set q-mesh using the same mesh points as G0 (green.mesh).
        # -----------------------------------------------------------
        ph_pure.set_mesh( green.mesh, is_mesh_symmetry=False, is_eigenvectors=True )
        ph_imp.set_mesh( green.mesh, is_mesh_symmetry=False, is_eigenvectors=True )
         
        # ------------------------------------------------
        # get tetrahedron meshes for the both system
        # ------------------------------------------------
        thm_pure, qpoints, f2P_all, _ = get_tetrahedron_mesh_f2( ph_pure, green.mesh )
        thm_imp, _, f2I_all, _ = get_tetrahedron_mesh_f2( ph_pure, green.mesh )
        
        # -------------------------------------
        # get supercells and check their size
        # -------------------------------------
        prim_pure = ph_pure.get_primitive()
        prim_imp  = ph_imp.get_primitive()
        sc_pure = ph_pure.get_supercell()
        sc_imp  = ph_imp.get_supercell()
        self._check_cell( sc_pure.get_cell(), sc_imp.get_cell() )
        
        # -----------------------------------------------
        # get atom positions with the unit of angstrom
        # and check positions of the pure system
        # -----------------------------------------------
        self.pos_pure = np.dot(sc_pure.scaled_positions, sc_pure.cell)
        self.pos_imp = np.dot(sc_imp.scaled_positions, sc_imp.cell)
        self._check_number_atoms( len(self.pos_pure), len(self.pos_imp) )
        
        self._check_positions( green.positions, self.pos_pure )
        
        # ------------------------------------------
        # get interatomic force constants (IFCs)
        # ------------------------------------------
        
        IFCs = ph_pure.get_force_constants()
        print(IFCs.shape)
        exit()


