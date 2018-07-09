# -*- coding: utf-8 -*-
import numpy as np

class Grid:
    def __init__(
            self, thm=None, qs=None, 
            f2s=None, evecs=None
            ):
        """Class to put q-grid information
        
        thm : TetrahedronMethod
            Containing info. of Tetrahedron grid
            Integration weights can be obtained with this.
        qs : ndarray, float, shape=(nq, 3)
            q-points
        f2s : array, float, shape=(nq, nmodes)
            Frequencies**2 on every q-points
        evecs : ndarray, float, shape=(3*nat_prim, nmodes)
            Eigenvectors on every q-points.
            Note that the eivenvector of the i-th mode is
            evecs[:,i], but NOT evecs[i,:].
        """
        self.thm = thm
        self.qs = qs
        self.f2s = f2s
        self.evecs = evecs
    
    def set_grid(
            self, thm, qs,
            f2s, evecs
            ):
        self.thm = thm
        self.qs = qs
        self.f2s = f2s
        self.evecs = evecs

class Newcell:
    """Class to treat the cell used for the calculation of T-matrix 
    """
    def __init__( self,
                  natoms      = None, 
                  cell        = None,
                  ids_orig    = None,
                  scaled_displacements  = None,
                  scaled_positions  = None,
                  chem_syms   = None):
        """
        Note
        -------
        This is a class for a new cell created to calculate T-matrix. 
        To calculate T-matrix, it is necesarry to put the impurity on 
        the cetner of the cell. Therefore, the impurity site
        is on the center of the new cell. In addition, when Phonopy treats
        a super cell, it considers multiplicity of atoms; equilibrium 
        atoms from a certain atom. However, during the calculation of
        T-matrix, multipilicity is always "1" and every atoms considered 
        appear explicitly. Consequently, the new cell contains every 
        equilibrium atoms and the number of atoms in the new cell is 
        generally larger than that of the original super cell.
        Here, "original cell" denotes the original super cell and 
        "new cell" the cell for the T-matrix calculation.
        
        Variables
        ------------
        nat : integer
            # of atoms
        
        cell : ndarray, float, shape=(3,3)
            Cell vectors
        
        ids_orig : array, integer, shape=(nat,)
            List of atom indices of atoms in the original supercell 
            corresponding to those in the new cell
        
        n2s_map : array, integer, shape=(nat,)
            # of equilibrium atoms
        
        scaled_displacements: narray, float, shape=(nat,3)
            Vectors from positions in the original cell to those in this cell.
            While dtype is "float", this should be the array of "integer".
        
        scaled_positions: narray, float, shape=(nat,3)
            Atomic positions in the new cell scaled by the cell vectors
        
        chem_syms : string, shape=(nat,)
            Chemical symbols such as H, He, ...
        
        IFCs : narray, shape=(nat, nat)
            Interatomic force constants
        
        masses : array, float, shape=(nat,)
            Mass of atoms

        # -- set later
        imp_site : integer
            Impurity site in the new cell

        """
        self.ids_orig   = None
        self.nequivs    = None
        self.scaled_displacements = None
        self.scaled_displacements = None
        self.imp_site   = None
        
        self.nat       = None
        self.cell      = None
        self.chem_syms = None
        self.masses    = None

        self.IFCs      = None
         
        if natoms is not None:
            self.nat = natoms
        if cell is not None:
            self.cell = cell
        if ids_orig is not None:
            self.ids_orig = ids_orig
        if scaled_positions is not None:
            self.scaled_positions = scaled_positions
        if scaled_displacements is not None:
            self.dips_from_orig = scaled_displacements
        if chem_syms is not None:
            self.chem_syms = chem_syms
        
    def set_crystal(
            self, ids_orig, nequivs, 
            scaled_displacements, r_imp, sc_orig):
        """Set variables of the class.
        See the above descriptions for parameters.
        """
        self.nat        = len(ids_orig)
        self.cell       = sc_orig.cell
        self.ids_orig   = ids_orig
        self.nequivs    = nequivs
        self.scaled_displacements = scaled_displacements
        self.scaled_positions = (
                sc_orig.get_scaled_positions()[ids_orig] + 
                scaled_displacements)
        
        chems = sc_orig.get_chemical_symbols()
        self.chem_syms = []
        self.masses = np.zeros(self.nat)
        for ia in range(self.nat):
            self.chem_syms.append(chems[ids_orig[ia]])
            self.masses[ia] = sc_orig.masses[ids_orig[ia]]
        
        self.imp_site = self._get_impurity_site(
                self.cell, r_imp, 
                self.scaled_positions)
        
    def _get_impurity_site( self, cell, r_imp, scaled_positions, eta=1e-5 ):
        """Get impurity site
        Parameters:
            cell       : cell vectors
            r_imp      : scaled position of the impurity
            scaled_positions : scaled positions of the atoms in the cell
        Return:
            atom index of the impurity
        """
        for ia in range( len(scaled_positions) ):
            length = np.linalg.norm(np.dot(scaled_positions[ia] - r_imp, cell))
            if length < eta:
                return ia
        return None

    def mkxyzfile( self, OFILE ):
        """
        Make xyz file. Use for the check.
        """
        coords = np.dot( self.scaled_positions, self.cell )
        ofs = open(OFILE, "w")
        ofs.write("{}\n".format(self.nat))
        for j in range(3):
            ofs.write("0.0 {0:15.8f} ".format(self.cell[j,j]))
        ofs.write("\n")
        for ia in range(self.nat):
            if self.imp_site is not None and ia == self.imp_site:
                ofs.write("Ar ")
            else:
                ofs.write("{} ".format(self.chem_syms[ia]))
            for j in range(3):
                ofs.write("{0:15.8f} ".format( coords[ia,j] ))
            ofs.write("\n")
        ofs.close()
        print("Output:", OFILE)

