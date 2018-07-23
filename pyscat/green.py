# -*- coding: utf-8 -*-
import numpy as np
from phonopy import Phonopy
from pyscat.crystal.crys import Grid
from pyscat.utils.idata import Idata
from pyscat.utils.freq  import get_frequency2
from pyscat.latdynam.tetrahedron import get_tetrahedron_mesh_f2
from pyscat.latdynam.phonon import get_phonopy
from pyscat.crystal.map import get_corresp_sites
from pyscat.calc.tips import (get_f2s4integral,
                              cal_green_function)
from pyscat.calc.dos import get_dos_green

"""
Note:
    Calculation of Green's function of a pure crystal
"""
class Green():
    def __init__(
            self, mesh, 
            fposcar="POSCAR", primat=None, ncells=None,
            fforce="FORCE_SET", fborn=None,
            ndiv=200
            ):
        """
        Parameters
        -----------
        mesh : array, integer, shape=(3,)
            # of q-meshes, (n1, n2, n3)
        fposcar : string
            POSCAR file name, will be named "unitcell"
        fforce : string
            File name of force and displacement sets.
            File format is the same as that used for Phonopy
        fborn : string
            File for Born effective charge, if exists
        primat : ndarray, float, shape=(3,3)
            primitive vector scaled by the unit cell
            ex) [[0.0,0.5,0.5], [0.5,0.0,0.5], [0.5,0.5,0.0]]
        ncells : array, integer, shape=(3,)
            # of unit cells in the supercell used for the force sets
        ndiv : integer
            # of frequencies for the integration of the real 
            part of the term excluding the wavefunction
        
        Variables
        -----------
        phonon : Phonopy
        grid : Grid
            Containing thm, qs, f2s, and evecs.
            See descriptions of Grid

        nq : integer
             = len(q_mesh)
        nat_prim : integer
            # of atoms in a primitive cell.
        nmodes : integer
            # of modes, = 3*nat_prim.
        
        g0 : ndarray, complex, shape=(3*natoms, 3*natoms)
            Green's function!
            "natoms" depends on the structure size
            at least equal to or more than # of atoms
            in a supercell, nat_super.
        f2_target : float
            frequency**2 at which g0 is calculated 
        f2s_integral : array, float, shape=(ndiv,)
            frequency**2 list used for an integral,
            ndiv = len(f2s_integral)
        dos_green : float
            DOS at f2_target calculated with Green's function
          
        """
        #--- input files
        self.idat = Idata(
                mesh    = mesh,
                fposcar = fposcar,
                fforce  = fforce,
                fborn   = fborn,
                primat  = primat,
                ncells  = ncells
                )
        
        self.ndiv = ndiv
         
        #--- for q-meshes
        self.phonon    = None
        self.grid      = Grid()
         
        #--- for Green's function
        self.g0           = None
        self.f2_target    = None
        self.f2s_integral = None
        self.dos_green    = None
        
    def set_qmesh(self, mesh=None):
        """
        Calculator of phonon properties at every q-points
        """
        if mesh is not None:
            self.idat.mesh = mesh
        self.phonon = get_phonopy(self.idat)
         
        self.phonon.set_mesh(
                self.idat.nmesh, 
                is_mesh_symmetry=False, 
                is_eigenvectors=True
                )
        
        thm, qs, f2s, evecs = get_tetrahedron_mesh_f2(
                self.phonon, self.idat.nmesh
                )
        self.grid = Grid(thm=thm, qs=qs, f2s=f2s, evecs=evecs)

        # -- Read number of atoms in the primitive and super cells
        self.idat.set_natoms(self.phonon)
    
    def cal_green_function(self, frequency=None, f2=None):
        """
        Note
        -----
        Calculation of Green's function of a pure crystal 
        at a given frequency (f2_tar).
        Note that input is a squared frequency
        
        Parameters
        -----------
        frequency, f2 : float
            frequency [THz] and squared frequency [THz**2]

        Returns:
        ----------
        g0 : ndarray, complex, shape=((3*natom, 3*natom)),
             where natom is # of atoms in a supercell
         
        """
        self.f2_target = get_frequency2(frequency=frequency, f2=f2)
        
        pc = self.phonon.get_primitive()
        sc = self.phonon.get_supercell()
        s2pp = get_corresp_sites( pc )
        
        self.g0 = cal_green_function(
                self.f2_target,
                self.grid,
                pc.cell,
                sc.cell,
                sc.get_scaled_positions(),
                s2pp,
                ndiv_integral = self.ndiv)
        
        self.dos_green = get_dos_green(self.g0) * 3*self.idat.nat_prim
        

