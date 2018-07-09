# -*- coding: utf-8 -*-
import numpy as np
from phonopy import Phonopy

#---- developped codes
#from phtmat.tetrahedron import get_tetrahedron_mesh_f2
from phtmat.utils.idata import Idata
from phtmat.utils.freq  import set_frequency2
from phtmat.latdynam.tetrahedron import get_tetrahedron_mesh_f2
from phtmat.latdynam.phonon import get_phonon

"""
Note:
    Calculation of Green's function of a pure crystal
"""
class Green():
    def __init__(self, mesh, 
                 fposcar="POSCAR", primat=None, ncells=None,
                 fforce="FORCE_SET", fborn=None,
                 ndiv=None):
        """
        Parameters:
            mesh    : # of q-mesh. shape=(3,)
            fposcar : POSCAR file (unit cell)
            fforce  : file of force and displacement sets
            fborn   : Born effective charge, if exists
            primat  : primitive vector scaled by the unit cell
            ncells  : # of unit cell in the supercell used for the force sets
            ndiv    : number of frequency points for the integration of the real part of
                      the term excluding the wavefunction
        
        Variables:
            phonon    : Phonopy
            thm       : containing info of Tetrahedron grids
                        Integration weights can be obtained.
            q_mesh    : q-meshes
            f2_mesh   : frequencies**2 on every q-points
            evec_mesh : eigenvectors on every q-points
            
            g0          : Green's function!
            f2_target   : frequency**2 at which g0 is calculated 
            f2s_integral: frequency**2 list used for an integral
            dos_normal  : DOS calculated with integration weights. shape=(ndiv,)
            dos_green   : DOS calculated with Green's function. (scalar)
            
            _df2        : interval of f2s_integral

        """
        #--- input files
        self.idat = Idata( mesh    = mesh,
                           fposcar = fposcar,
                           fforce  = fforce,
                           fborn   = fborn,
                           primat  = primat,
                           ncells  = ncells)
         
        if ndiv is not None:
            self.ndiv = ndiv
        else:
            self.ndiv = 200
         
        #--- for q-meshes
        self.phonon    = None
        self.thm       = None
        self.q_mesh    = None
        self.f2_mesh   = None
        self.evec_mesh = None
        
        #--- for Green's function
        self.g0           = None
        self.f2_target    = None
        self.f2s_integral = None
        self.dos_normal   = None
        self.dos_green    = None
        
        self._df2 = None
        
    def set_phonon_properties(self):
        """
        Calculate phonon properties at every q-points
        """    
        self.phonon = get_phonon( self.idat )
         
        self.phonon.set_mesh( self.idat.nmesh, is_mesh_symmetry=False, is_eigenvectors=True)
        
        self.thm, self.q_mesh, self.f2_mesh, self.evec_mesh \
                = get_tetrahedron_mesh_f2( self.phonon, self.idat.nmesh )
        
        # -------------------------------------------------------
        # Read number of atoms in the primitive and super cells
        # -------------------------------------------------------
        self.idat.set_natoms( self.phonon )
     
    def cal_green_function( self, frequency=None, f2=None ):
        """
        Note:
            Calculation of Green's function of a pure crystal at a given frequency (f2_tar).
            Note that input is a squared frequency
        
        Parameters:
            f2_tar: frequency**2
            
        """
        self.f2_target = set_frequency2( frequency=frequency, f2=f2 )
        
        # ------------------------------------
        # calculate integration weights (IW)
        # ------------------------------------
        iw4real, iw4imag = self._get_integration_weights_for_real_imag( self.thm )
         
        # ------------------------------------
        # calculate Green's function!!
        # ------------------------------------
        self.g0 = np.zeros((3*self.idat.nat_super, 3*self.idat.nat_super), dtype=np.complex_)
        for iq in range( len(self.q_mesh) ):
            self.g0 += self._cal_green_at_q( self.phonon,
                                             self.q_mesh[iq], self.evec_mesh[iq], 
                                             iw4real[iq], iw4imag[iq] )
        
        self.g0 /= float(3*self.idat.nat_super)
        
        self.get_dos_green()
    
    def _get_integration_weights_for_real_imag( self, thm ):
        """
        Calculate integration weights for the real and imaginary parts.
        DOS is also calculated at the same time.
        """
                  
        nqpoints = len(self.q_mesh)
        
        #- for real part: iw4imag[nqpoints x ndiv x nmod]
        self._get_f2points4integral( np.min(self.f2_mesh), np.max(self.f2_mesh) )
        thm.set(value='I', frequency_points=self.f2s_integral)
        iw4real = np.zeros(((nqpoints, self.ndiv, self.idat.nmod)))
        self.dos_normal = np.zeros(self.ndiv)
        for iq, iw in enumerate( thm ):
            iw4real[iq,:,:] = iw[:,:]
            self.dos_normal += np.sum(iw, axis=1)
        
        #- for imaginary part: iw4imag[nqpoints x nmod]
        thm.set(value='I', frequency_points=np.array([self.f2_target]))
        iw4imag = np.zeros((nqpoints, self.idat.nmod))
        for iq, iw in enumerate( thm ):
            iw4imag[iq,:] = iw[0,:]
        
        return iw4real, iw4imag
    
    def _cal_green_at_q( self, phonon, qpoint, evecs, iw4real, iw4imag ):
        """Calculte Green's function at q
        Parameters:
            phonon (Phonopy)
            qpoint: (q1,q2,q3)
            evec:
                eigenvector at q and a certain band.
                shape=(3*nsites, nmodes)
                Note that evecs[:,i] is the i-th eigenvector.
            iw4real, iw4imag:
                integration weights for the real and imaginary part of
                the term excluding the wavefunction
        
        """
        # --------------------------------------------------
        # Calculate real and imaginary parts of the term
        # excluding wave function
        # --------------------------------------------------
        Re_part, Im_part = self._cal_dos_part( iw4real, iw4imag )
        
        nmat = 3 * self.idat.nat_super
        g0q = np.zeros(( nmat, nmat ), dtype=np.complex_)
        
        for imod in range(self.idat.nmod):
            
            g0q += self._cal_green_at_q_j( qpoint, evecs[:,imod] ) \
                   * ( Re_part[imod] + 1j * Im_part[imod] )
        
        return g0q
    
    def _cal_green_at_q_j( self, qpoint, evec ):
        """Caluate Green's function at q and j-th band
        
        Parameters:
            evec:
                eigenvector at q and j-th band.
                shape=(3*nsites)
            others:
                see descriptions for "cal_green_at_q"
        
        """
        primitive = self.phonon.get_primitive()
        supercell = self.phonon.get_supercell()
        
        # -----------------------------------------------------
        # Extract smallest vectors.
        # p2p_map[Ns], p2s_map[Np], s2p_map[Ns]
        # svecs  : smallest vectors, shape=(Ns, Np, 27, 3).
        # multi  : multiplicitly,    shape=(Ns, Np)
        # Np, Ns: # of atoms in the primitive and super cells
        # -----------------------------------------------------
        p2p_map = primitive.get_primitive_to_primitive_map()
        #p2s_map = primitive.get_primitive_to_supercell_map()
        s2p_map = primitive.get_supercell_to_primitive_map()
         
        #------------------------------------------------------------------
        # Transform matrix from the real space to the primitive cell base
        #------------------------------------------------------------------
        M_r2p = np.linalg.inv( primitive.cell )
        
        # ---------------------------------------------
        # Get wavefunction expanded in the super cell:
        # [3*self._na_super] array
        # ---------------------------------------------
        nmat = 3 * self.idat.nat_super
        wf_expand = np.zeros(nmat, dtype=np.complex_)
        for ia in range(self.idat.nat_super):
             
            # convert a coordinate from real scale to that scaled by the primitive cell
            ri = np.dot( supercell.scaled_positions[ia], supercell.cell )
            ri = np.dot( ri, M_r2p )
             
            # get the site (in a primitive cell) of the ia-th atom (in a supercell)
            isite = p2p_map[ s2p_map[ia] ]
            
            wf_expand[ ia*3: (ia+1)*3 ] = \
                    np.exp( 1j * 2.*np.pi*np.dot( qpoint, ri )) \
                    * evec[3*isite:3*(isite+1)]
        
        g_qj = np.tile( wf_expand, nmat ).reshape( nmat, nmat ) \
                 * np.tile( np.conj(wf_expand), ( nmat, 1 ) )
        
        return g_qj
    
    def get_dos_green( self ):
        """Calculate phonon DOS using Green's function at a given frequency
        
        """
        self.dos_green = np.sum( np.diag( np.imag(self.g0) ) ) \
                * 3 * self.idat.nat_prim / np.pi
        
        return self.dos_green
    
    def _cal_dos_part( self, iw4real, iw4imag ):
        """
        Note:
            Calculate the real and imaginary parts
        
        Parameters:
            iw4real[nqpoints x nmod], iw4imag[nmod]:
                integration weights for real and imaginary parts
                nqpoints and nmod are number of q-points and modes at each q
        
        Return:
            Re_part and Im_part
        """
        #--- real part
        Re_part = self._cal_real_part( iw4real )
        
        #--- imaginary part: [nmod]
        Im_part = np.pi * iw4imag
         
        return Re_part, Im_part
    
    def _cal_real_part(self, iw4real):
        """
        Note:
            Calculate the real part of the term excluding wave function.
        
        Parameters:
            iw4real[self._ndiv x nmod]:
                delta function at w'^2 - w(q)^2
        
        Return:
            calculated value
        """
        self._df2 = self.f2s_integral[1] - self.f2s_integral[0]
        Re_part = np.zeros( self.idat.nmod )
        for i2 in range(self.ndiv):
            Re_part += self._df2 * iw4real[i2] / (- self.f2_target + self.f2s_integral[i2])
        
        return Re_part
     
    def _get_f2points4integral( self, f2min, f2max ):
        """
        Note:
            Get frequencies for the integration.
            Adjust the frequency points to make self.f2_target locate on 
            the middle of adjacent two f^2 points.
        
        Parameters:
            f2min, f2max: minimum and maximum frequencies
        
        Return:
            f2s_integral: frequencies used for the integration
        """
        self.f2s_integral = np.linspace(f2min, f2max, self.ndiv)
        if f2min < self.f2_target and self.f2_target < f2max:
            for i in range(self.ndiv-1):
                if self.f2s_integral[i] <= self.f2_target and self.f2_target <= self.f2s_integral[i+1]:
                    self.f2s_integral -= \
                            (self.f2s_integral[i] + self.f2s_integral[i+1]) * 0.5 - self.f2_target
        

