# -*- coding: utf-8 -*-
import numpy as np
from phonopy import Phonopy

from pyscat.utils.idata import Idata
from pyscat.utils.freq  import get_frequency2
from pyscat.green import Green
from pyscat.latdynam.phonon import (
        get_phonopy, 
        set_phonopy_structure,
        get_grid)
from pyscat.latdynam.tetrahedron import get_tetrahedron_mesh_f2
from pyscat.crystal.map import (
        get_corresp_sites,
        offset_cell_position,
        get_IFCs)
from pyscat.crystal.crys import Newcell, Grid
from pyscat.crystal.vesta import mkvesta_IFCs
from pyscat.utils.comp import (
        check_cells, check_number_of_atoms, check_positions)
from pyscat.calc.tips import (
        cal_green_function, get_dphi, get_UTU_braket,
        conv_tmat_L2s, get_dos_imp, conv_n2ss2weights)
from pyscat.calc.dos import get_dos_green
from pyscat.utils.files import (
        initialize_scattering_file,
        dump_scattering_rate,
        output_scattering_rates)

def _check_input_parameter(para, label):
    if para is None:
        print("Error.", label, "is not given.")
        import sys
        sys.exit()

"""  
"""
class Tmat(object):
    def __init__(self,
                 mesh = None,
                 impurity_site = None,
                 primat_pure   = None,
                 ncells_pure   = None,
                 fposcar_pure = "POSCAR_pure", 
                 fforce_pure  = "FORCE_SET_pure",
                 fborn_pure   = None,
                 primat_imp   = np.eye(3),
                 ncells_imp   = np.eye(3, dtype=int),
                 fposcar_imp  = "POSCAR_imp",
                 fforce_imp   = "FORCE_SET_imp",
                 fborn_imp    = None,
                 tol_disp     = 1.0
                 ):
        """
        Parameters
        -----------
        impurity_site : integer
            Atom site for the impurity 
            Currently, only one substitution is supported.
        ***_pure     : for pure crystal
        ***_imp      : for the structure with an impurity
        fposcar_*, fforce_*, fborn_* : string
            File names for POSCAR, force sets, and Born effective charge
        mesh, primat_*, ncells_* :
            # of meshes, primitive matrix, and # of cells
            Mesh size is used for the pure crystal.
        See also the description of Idata class.

        tol_disp : float, unit[A]
            tolerant displacement for the structure with an impurity

        Variables
        ------------
        ph_***          : phonon information
        _nat_prim_***   : # of atoms in the primitive cell
        _nat_super_***  : # of atoms in the super cell
        _nmod_***       : # of modes
    
        _flag_cell      : result of the check of cell size 
        _flag_natom     : result of the check of number of atoms
        _flag_positions : result of the check of atomic positions
          
        Note
        ------------
        The order of atoms in the two POSCAR file and the structure used for
        the force sets should correspond to each other.
        """
        _initialize_pyscat()
        _check_input_parameter( mesh, "q-mesh (mesh)" )
        _check_input_parameter( primat_pure, 
                "primitive matrix of the pure crystal (primat_pure)" )
        _check_input_parameter( ncells_pure, 
                "# of unit cells in the pure crystal (ncells_pure)" )
        
        # ---------------------------
        # Set input parameters
        # ---------------------------
        # -- for pure crystal
        self.idat_pure = Idata(
                mesh     = mesh,
                fposcar  = fposcar_pure,
                fforce   = fforce_pure,
                fborn    = fborn_pure,
                primat   = primat_pure,
                ncells   = ncells_pure)
        
        # -- for the structure w/ an impurity
        self.idat_imp = Idata( 
                fposcar  = fposcar_imp,
                fforce   = fforce_imp,
                fborn    = fborn_imp,
                primat   = primat_imp,
                ncells   = ncells_imp)
        
        # -- impurity site
        self.imp_site = impurity_site
        
        # Tolerance of displacement of atoms in the pure crystal
        # and in the structure with an impurity
        self.tol_disp = tol_disp
        
        # -------------------------------------------
        # Phonopy class for the pure crystal and
        # the structure with an impurity
        # -------------------------------------------
        self.ph_pure = None
        self.ph_imp  = None
        
        # -----------------------------------
        # variables that will be used
        # -----------------------------------
        self.f2_target = None
        self.g0 = None
        
        self._flag_structure = None
        self.flag_calc = None
        
        # ----------------------------------------------------
        # new cell around an impurity in which the impurity
        # is located on the center.
        # ----------------------------------------------------
        self.nc_pure = Newcell()
        self.nc_imp  = Newcell()
        self.grid_pure = Grid()
        self.grid_imp = Grid()
    
    def set_ifcs_and_phonon(self):
        
        # -- for the pure crystal
        self.ph_pure = get_phonopy(self.idat_pure)
        self.ph_pure.set_mesh(
                self.idat_pure.nmesh, 
                is_mesh_symmetry=False, 
                is_eigenvectors=True) 
        
        # -- for the structure w/ an impurity
        self.ph_imp = get_phonopy( self.idat_imp )
          
        self.idat_pure.set_natoms(self.ph_pure)
        self.idat_imp.set_natoms(self.ph_imp)
         
        
    def check_structures(
            self, natom=True, cell=True, position=True):
        """ Check the structures.
        Note
        ------
        What is checked may be changed depending on what kind of impurities
        can be calculated with this code.
        
        Parameters
        ------------
        unit1 and unit2: Supercell and Primitive
            See descriptions in Phonopy
        """
        unit1 = self.ph_pure.get_supercell()
        unit2 = self.ph_imp.get_supercell()
        
        flag_cell = check_cells(
                unit1.cell, unit2.cell, flag=cell, tolerance=1e-5)
        
        flag_natom = check_number_of_atoms(
                len(unit1.masses), len(unit2.masses), flag=natom)
        
        flag_positions = check_positions(
            np.dot(unit1.scaled_positions, unit1.cell),
            np.dot(unit2.scaled_positions, unit2.cell),
            cell=unit1.cell,
            tolerance=self.tol_disp,
            flag=position
            )
        
        if flag_positions is True:
            self.idat_pure.scID = np.arange(len(unit1.masses))
            self.idat_imp.scID  = np.arange(len(unit2.masses))
        else:
            print("Error. This is not supported yet.")
            print("Please make POSCAR files with the same order of atoms.")
            import sys
            sys.exit()
        
        print(" Given structures are compatible.")
    
    def set_newcelles_around_impurity(self):
        """Identify atoms around the impurity and set new cells around it.
        Note
        ----------
        This identifies in which supercell the atoms are located around the
        impurity, which roughly means that putting the impurity on the center of
        the supercell. After this procedure, atom indices for G0 should be also
        modified. Here, "new cell" denotes the cell created around the impurity,
        which is almost the same as the supercell, but may contain more atoms
        that the supercell because there may exist equilibrium atoms.
        
        Return
        ----------
        ids_orig : array, integer, shape=(natoms_new,)
            List of atom IDs of the supercell corresponding to 
            atoms in the new cell
        scaled_displacements : narray, float, shape=(natoms_new,3)
            Scaled displacement vectors from the original points,
            which is float, but should be integer
        natoms_new : integer
            # of atoms around the impurity,
            should be >= 3*nat_super,
        nat_super : integer
            # of atoms in a supercell
        
        For example, if id_list[inew] = iorig, the iorig-th atom in the 
        original cell is displaced with scaled_displacements[iorig] 
        to make the new cell.
        """
        sc_imp = self.ph_imp.get_supercell()
        poses_imp = sc_imp.get_scaled_positions()
        
        # -- Position of the impurity, which can be changed to 
        # -- an arbitrary position.
        R_imp = poses_imp[self.imp_site]
        
        ids_orig, nequivs, disps_from_orig = (
                offset_cell_position(sc_imp.cell, poses_imp, R_imp))
        
        self.nc_imp.set_crystal(
                ids_orig, nequivs, disps_from_orig, R_imp, sc_imp)
        self.nc_pure.set_crystal(
                ids_orig,
                nequivs,
                disps_from_orig,
                R_imp, 
                self.ph_pure.get_supercell()
                )
        
        #self.nc_imp.mkxyzfile( "imp.xyz" )
        #self.nc_pure.mkxyzfile( "pure.xyz" )
        
    def set_ifcs4newcells( self ):
        """Set force constants for the new cells,
        the pure structure and the structure with an impurity.
        
        Return
        -----------
        IFCs : ndarray, float, shape=(nat, nat, 3, 3)
            Force constants.
            Unit : eV/A^2, same as Phonopy.
        """
        sc_pure = self.ph_pure.get_supercell()
        sc_imp = self.ph_imp.get_supercell()
        
        self.nc_pure.IFCs = get_IFCs(
                self.ph_pure, 
                self.nc_pure.ids_orig, 
                self.nc_pure.scaled_displacements
                )
        self.nc_imp.IFCs = get_IFCs(
                self.ph_imp, 
                self.nc_imp.ids_orig, 
                self.nc_imp.scaled_displacements
                )
        
        print(" Force constants are set.")
        
    def visualize_IFCs_diff(
            self, IFCs=None, idir1=0, idir2=0, IFC_max=0.1, IFC_min=-0.1
            ): 
        """Visualization of the difference of IFCs
        Parameters:
            idir{i-th}  : integer
                direction for i-th atom (0.x, 1.y, 2.z)
            IFC_{min, max}: float
                min and max IFC for the visualization
        """
        if IFCs is None:
            IFCs_view = self.nc_imp.IFCs - self.nc_pure.IFCs
        else:
            IFCs_new = IFCs
        OFILE = "ifcs_diff_{:d}{:d}.vesta".format(idir1, idir2)
        mkvesta_IFCs( OFILE,
                      IFCs_view,
                      self.nc_pure.scaled_positions,
                      self.nc_pure.imp_site, idir1, idir2, 
                      IFC_max = IFC_max, IFC_min = IFC_min)
        
    def set_green_pure(
            self, frequency=None, f2=None, 
            ndiv_integral=200, feta=1e-3):
        """Calculator of Green's function of the pure crystal
        Parameters
        ------------
        frequency, f2 : float
            Frequency [THz] or f2[TH**2] at which Green's function
            will be calculated.
        
        ndiv_integral : integer
            # of frequency points for the integral
        """
        self.f2_target = get_frequency2(frequency=frequency, f2=f2)
        if np.sqrt(abs(self.f2_target)) < feta:
            self.flag_calc = False
            return False
        else:
            self.flag_calc = True
        
        # -- get tetrahedron meshes for the both system
        ph_pure = self.ph_pure
        mesh    = self.idat_pure.nmesh
        
        # -- get phonons at every q-points
        thm, qs, f2s, evecs = get_tetrahedron_mesh_f2(ph_pure, mesh)   
        self.grid_pure.set_grid(thm, qs, f2s, evecs)
         
        # -- extract Primitive cell
        pc = self.ph_pure.get_primitive()
        
        n2pp = _get_IDmap4newcell(pc, self.nc_pure)
        ndiv = ndiv_integral
        
        # -- calculate Green's function!!
        self.g0 = cal_green_function(
                self.f2_target, self.grid_pure,
                pc.cell, self.nc_pure.cell,
                self.nc_pure.scaled_positions,
                n2pp, ndiv_integral = ndiv)
        #print(" Green's function is calculated.")
    
    def set_Tmatrix(self):
        """Calculator of T-matrix
        """
        if self.flag_calc is not True:
            return None
        dphi = self._get_dphi()
        Imat = np.eye(len(self.g0))
        self.tmat = np.dot(
                np.linalg.inv(Imat - np.dot(dphi, self.g0)), 
                dphi
                )
        
    
    def autoloop4tmat(self, ndiv_integral=200):
        """Calculate scattering rate at every irreducible q-point and mode.
        Parameters
        -----------
        qs : ndarray, float, shape=(nq, 3)
            q-points
        weights : array, float, shape=(nq,)
            weights of each q-point
        f2s : ndarray, float, shape=(nq, nmodes)
            frequencies
        evecs : ndarray, complex, shape=(nq, 3*nsites, nmodes)
            eigenvectors
        dos_imp : ndarray, float, shape=(nq, nmodes)
            DOS of the system with an impurity
        nq, nsites, and nmodes : integer
            # of q-points, sites, and modes.
            nmodes == 3*nsites
        """
        qs, weights, f2s, evecs = self.get_grid4summation()
        nq = len(qs)
        nmodes = len(f2s)
        rscat = np.zeros_like((f2s), dtype=float)
        dos_imp = np.zeros_like((f2s), dtype=float)
        ddump = np.zeros_like((f2s), dtype=np.complex)
        
        for iq, qpoint in enumerate(qs):
            for im, f2 in enumerate(f2s[iq]):
                self.set_green_pure(f2=f2, ndiv_integral=ndiv_integral)
                self.set_Tmatrix()
                dos_imp[iq,im] = self.get_dos_tmat()
                rscat[iq,im], ddump[iq,im] = (
                        self.get_scattering_rate(
                            qpoint=qpoint,
                            f2=f2,
                            evec=evecs[iq,:,im],
                            iq=iq, imode=im
                            )
                        )
                dump_scattering_rate(iq, qpoint, im, f2, rscat[iq,im])
        _end_of_simulation()
        #output_scattering_rates(qs, f2s, rscat)
        
        frequencies = np.sqrt(abs(f2s)) * np.sign(f2s)
        return qs, frequencies, rscat, dos_tmat

    def get_grid4summation(self):
        """Get phonon properteis on grid
        """
        qs, weights, f2s, evecs = get_grid(self.idat_pure)
        
        mesh = self.idat_pure.nmesh
        print("")
        print(" q-mesh : {:d} x {:d} x {:d}".format(mesh[0], mesh[1], mesh[2]))
        print(" # of irreducible q-points : {:d}".format(len(qs)))
        print("")
        initialize_scattering_file(mesh, len(qs))

        return qs, weights, f2s, evecs
        
    def get_scattering_rate(
            self, qpoint=None, f2=None, evec=None, iq=None, imode=None):
        """Calculator of the scattering rate
        Parameters
        -----------
        qpoint : float
            scaled wavevector
        f2 : float
            Squared frequency
        evec : array, float, shape=(nmodes)
            Eigenvector
        iq : integer (optional)
            q-point index
        imode : integer (optional)
            mode index

        nmodes : integer
            # of modes
        
        Return
        ---------
        rscat : float
            Scattering rate
        """
        if self.flag_calc is not True:
            rscat = 0.0
            _print_rscat(iq, qpoint, imode, f2, rscat)
            return rscat, None
        if qpoint is None:
            print("Error: qpoint should be given.")
            import sys
            sys.exit()
        
        #====== UNIT !! =====
        evecs = self.grid_pure.evecs
        pc = self.ph_pure.get_primitive()
        nc = self.nc_pure
        n2pp = _get_IDmap4newcell(pc, self.nc_pure)
        nat_super = len(n2pp)
        
        UTU = get_UTU_braket(
                self.tmat, qpoint=qpoint, f2=f2, evec=evec,
                pcell=pc.cell,
                scell=nc.cell,
                scaled_positions=nc.scaled_positions,
                s2pp=n2pp)
        
        # --- calculate the scattering rate
        freq = np.sqrt(abs(f2)) * np.sign(f2)
        real_imag = UTU / float(freq * nat_super)  
        self.rscat = np.imag(real_imag)
         
        _print_rscat(iq, qpoint, imode, f2, self.rscat)
        return self.rscat, real_imag

    def _get_dphi(self):
        """Calculator of dphi, phi is force constants
        Return
        -------
        dphi : ndarray, float, shape=(3*nat, 3*nat)
            matrix representing the difference of IFCs
        nat : integer
            # of atoms
        """
        dphi = get_dphi(
                self.f2_target, 
                self.nc_pure.masses, 
                self.nc_imp.masses, 
                self.nc_pure.IFCs, 
                self.nc_imp.IFCs
                )
        return dphi
    
    def get_dos_tmat(self):
        """Calculator of DOS with T-matrix
        """
        if self.flag_calc is not True:
            dos = 0.0
            return dos
        stmat = conv_tmat_L2s(self.tmat, self.nc_imp.masses)
        n2ss = _get_IDmap4newcell(self.ph_imp.get_primitive(), self.nc_imp)
        multi = conv_n2ss2weights(n2ss)
        dos = get_dos_imp(stmat, self.g0, multipilicity=multi)
        return dos

    def get_dos_pure(self):
        """Calculate DOS of the pure system
        """
        if self.flag_calc is not True:
            dos = 0.0
            return dos
        primitive = self.ph_pure.get_primitive()
        n2pp = _get_IDmap4newcell(primitive, self.nc_pure)
        multi = conv_n2ss2weights(n2pp)
        nat_prim = len(primitive.masses)

        print(n2pp)
        print(multi)
        dos = get_dos_green(self.g0, multiplicity=multi)
        return dos
        
def _initialize_pyscat():
    print("")
    print(" -----------------------------------------------------------------")
    print(" pyscat : Calculator of phonon scattering rate due to an impurity")
    print(" -----------------------------------------------------------------")
    print("")

def _get_IDmap4newcell(pc, nc):
    """Get atomic ID mapping for the new cell.
    Parameters
    -----------
    pc : Primitive
    nc : Newcell
    
    Return
    --------
    s2pp : array, integer, nat_new
        Atomic indices of the primitive cell corresponging to 
        atoms in the new cell
    """
    s2pp = get_corresp_sites(pc)
    n2pp = np.array([s2pp[i] for i in nc.ids_orig], dtype=int)
    return n2pp

def _print_rscat(iq, qpoint, imode, f2, rscat):
    """Print a calculated scattering ratio
    iq, imode : integer
        q-point and mode index
    qpoint : array, float, shape=(3)
        q-point
    f2 : float
        Squared frequency
    rscat : float
        Scattering rate
    """
    if iq is None or imode is None:
        return False
    if imode == 0:
        print("")
        print(" {:2d}-th q-point : {:6.2f} {:6.2f} {:6.2f}".
                format(iq, qpoint[0], qpoint[1], qpoint[2]))
    
    freq = np.sqrt(abs(f2)) * np.sign(f2)
    print(" {:2d} {:6.2f} THz : ".format(imode, freq), end="")
    if rscat is not None:
        #print("{:13.3e}".format(rscat))
        print("{:13.8f}".format(rscat))
    else:
        print(" None")

def _end_of_simulation():
    print("")
    print(" All done!")
    print("----------------------------------------")
    print("")

