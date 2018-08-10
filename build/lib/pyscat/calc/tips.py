#
# Functions that are mainly used for a calculation of 
# Green's function and T-matrix
#
# -*- coding: utf-8 -*-
import numpy as np
import phonopy.units

AMU = 1.6605402e-27 # [kg]
EV = 1.60217733e-19 # [J]
Angstrom = 1.0e-10  # [m]
VaspToTHz = np.sqrt(EV/AMU)/Angstrom/(2*np.pi)/1e12 # [THz] 15.633302


def cal_green_function(
        f2_target, grid,
        pcell, scell,
        scaled_positions,
        s2pp, ndiv_integral = 200):
    """Calculator of Green's function
    Note
    -----
    Calculation of Green's function of a pure crystal 
    at a given frequency (frequency or f2_tar).
    
    Parameters
    -----------
    f2_target : float
        frequency [THz**2] at which Green's function
        is being calculated.
    grid : Grid
        thm (ThetrahedronMethod), qs (q-points), f2s
        (frequencies), and evecs (eigenvectors) are
        contained.
    pcell, scell : ndarray, float, shape=(3, 3)
        Cell vectors of the primitive and super cells
    scaled_positions : ndarray, float, shape=(nat_super, 3)
        Atom positions in the supercell scaled by scell
    s2pp : array, integer, shape=(nat_super,)
        List of atom indices in the primitive cell
        corresponding to those in the super cell
    ndiv_integral : integer
        Frequency points fore the integral
     
    nat_super : integer
        # of atoms in the supercell
        = len(scaled_positions)
     
    Returns
    ----------
    g0 : ndarray, complex, shape=(3*nat_super, 3*nat_super)
        Green's function of the pure crystal
    
    """
    nqpoints  = len(grid.qs)
    nmodes    = len(grid.evecs[0,0])
    nat_prim  = int(nmodes / 3)
    nat_super = len(scaled_positions)
    
    # --- get frequencies
    f2s_integral = (
            get_f2s4integral(
                f2_target,
                np.min(grid.f2s),
                np.max(grid.f2s),
                ndiv_integral
            )
        )
    
    # -- calculate integration weights (IW)
    iw4real, iw4imag = get_iws4both(
            grid.thm, f2s_integral, 
            f2_target, nqpoints, nmodes
            )
    
    # -- calculate DOS using IWs
    #dos_normal = get_dos_iws(iw4real)
    
    # -- calculate Green's function!!
    g0 = np.zeros((3*nat_super, 3*nat_super), dtype=np.complex)
    for iq in range(len(grid.qs)):
        g0 += cal_green_at_q(
                grid.qs[iq], grid.evecs[iq],
                iw4real[iq], iw4imag[iq],
                f2_target, f2s_integral,
                pcell, scell, scaled_positions, s2pp
                )
    
    g0 /= float(3. * nat_super)  #==== this coeff. is to be checked!!

    return g0

def get_iws4both(thm, freq1s, freq2, nqpoints, nmodes):
    """Get integration weights (IWs) for two sets of frequencies.
    Obtain IWs at freq1s and freq2.

    Parameters
    -------------
    thm : TetrahedronMethod
        See descriptions in Phonopy for details.
    freq1s : array, float, shape=(len(freq1s),)
        List of frequencies
    freq2  : float
        A frequency
    nqpoints : integer
        # of q-points
    nmodes : integer
        # of modes

    Returns
    -----------
    iws1 : ndarray, float, shape=(((nqpoints, len(freq1s), nmodes)))
    iws2 : ndarray, float, shape=((nqpoints, nmodes))
        IWs at freq1s and freq2, respectively.
        
    """
    # --- for real part
    iws1 = extract_iws(thm, freq1s, nqpoints, nmodes)
    
    # --- for imaginary part
    iw_temp = extract_iws(thm, np.array([freq2]), nqpoints, nmodes)
    iws2 = np.zeros((nqpoints, nmodes))
    iws2[:,:] = iw_temp[:,0,:]
    
    return iws1, iws2

def extract_iws(thm, fpoints, nqpoints, nmodes):
    """Calculate integration weights at given frequencies, fpoints.
    
    Parameters
    -----------
    thm : TetrahedronMesh
    fpoints : ndarray, float, shape=(nfreq,)
        List of frequencies at which IW will be calculated
    nqpoints : integer
        # of q-points
    nmodes : integer
        # of modes
    
    Return
    -----------
    iw : ndarray, float, shape=(((nqpoints, ndiv, nmodes)))
        Integration weights
    
    """
    nfreq = len(fpoints)
    thm.set(value='I', frequency_points=fpoints)
    iw_dump = np.zeros(((nqpoints, nfreq, nmodes)))
    for iq, iw in enumerate(thm):
        iw_dump[iq,:,:] = iw[:,:]
    
    return iw_dump

def cal_green_at_q( qpoint, evecs,
                    iw4real, iw4imag,
                    f2_target, f2s_integral, 
                    pcell, scell, scaled_positions, s2pp ):
    """Calculte Green's function at q
    Parameters
    -----------
    qpoint : ndarray, float, shape=(3,)
        wavevector, q = (q1, q2,q3 )
    evec : ndarray, shape=(3*nsites, nmodes)=(nmodes, nmodes)
        eigenvector at q and a certain band.
        Note that evecs[:,i] is the i-th eigenvector.
    iw4real, iw4imag : ndarray, float, shape=((nqpoints, nmodes))
        integration weights for the real and imaginary part of
        the term excluding the wavefunction
    
    nqpoints : integer
        # of q-points
    nsites : integer
        # of atoms in a primitive cell
    nmodes : integer
        # of modes, which should be equal to 3*nsites
    
    Returns
    ---------
    g0q : ndarray, complex, shape=(nmat, nmat)
        Green function at q
    """
    # --------------------------------------------------
    # Calculate real and imaginary parts of the term
    # excluding wave function
    # --------------------------------------------------
    Re_part, Im_part = cal_dos_part(iw4real, iw4imag, f2_target, f2s_integral)
    
    nat_super = len(scaled_positions)
    nmat = 3 * nat_super
    g0q = np.zeros((nmat, nmat), dtype=np.complex)
    
    nmodes = len(evecs[0])
    for imod in range(nmodes):
        
        g0q += (cal_green_at_qj(qpoint, evecs[:,imod], 
                                pcell, scell,
                                scaled_positions, s2pp) 
                * (Re_part[imod] + 1j * Im_part[imod]))
    
    return g0q

def cal_green_at_qj( qpoint, evec, pcell, scell, scaled_positions, s2pp ):
    """Caluation of Green's function at q and j-th band
    
    Parameters
    -----------
    evec : ndarray, complex, shape=(3*nsites)
        eigenvector of j-th band at q.
    {p,s}cell: ndarray, shape=((3,3))
        cell vectors of the primitive and super cell
    scaled_positions: ndarray, shape=((nat_super,3))
        Atomic positions in the supercell
    nat_super:
        number of atoms in the supercell

    """  
    wf_expand = get_expanded_wf(
            qpoint, evec, pcell, scell, scaled_positions, s2pp)
    nmat = 3 * len(scaled_positions)
    g_qj = (np.tile(wf_expand, nmat).reshape(nmat, nmat) *
            np.tile(np.conj(wf_expand), (nmat, 1)))  
    return g_qj

def get_expanded_wf(
        qpoint, evec,
        pcell, scell, scaled_positions, s2pp):
    """Make wavefunction expanded in the supercell
    Parameters
    ----------
    {p,s}cell : ndarray, float, shape=(3,3)
        Cell vectors of the primitive and super cells
    qpoint : array, float, shape=(3,)
        q-point scaled by the primitive BZ (pcell^-1)
    evec : arrya, complex, shape=(3*nat_prim)
        eigenvector of a certain band at q
    scaled_positions : ndarray, float, shape=(nat_super, 3)
        Atomic positions scaled by scell
    s2pp : array, integer, shape=(nat_super)
        List of atomic indices in a primitive cell
        corresponding to those in a super cell
    
    nat_{prim, super} : integer
        # of atoms in a primitive and super cells

    Return
    -------
    wf_expand : ndarray, complex, shape=(3*nat_super, 3)
    """
    nat_prim = int(len(evec)/3)
    nat_super = len(scaled_positions)
    nmat = 3 * nat_super
    
    # -- Transform matrix from the real space to the primitive cell base
    M_r2p = np.linalg.inv(pcell)
    
    # -- Get wavefunction expanded in the super cell
    wf_expand = np.zeros(3*nat_super, dtype=np.complex)
    for ia in range( nat_super ):
         
        # convert a coordinate from real scale to 
        # that scaled by the primitive cell
        ri = np.dot(np.dot(scaled_positions[ia], scell), M_r2p)
         
        # get the site (in a primitive cell) of 
        # the ia-th atom (in a supercell)
        isite = s2pp[ia]
        
        wf_expand[3*ia:3*(ia+1)] = (
                np.exp( 1j * 2.*np.pi*np.dot( qpoint, ri )) *
                evec[3*isite:3*(isite+1)]
                )
    return wf_expand 

def cal_dos_part(iw4real, iw4imag, f2_target, f2s_integral):
    """
    Note
    -----------
    Calculate the real and imaginary parts
    
    Parameters
    -----------
    iw4real : ndarray, double, shape=((nqpoints, nmodes))
    iw4imag : ndarray, double, shape=(nmodes,)
        Integration weights for real and imaginary parts
    
    nqpoints : integer
        # of q-points
    nmodes : integer
        # of modes
    
    Return
    -----------
    Re_part and Im_part
    """
    #--- real part
    Re_part = cal_real_part( iw4real, f2_target, f2s_integral )
    
    #--- imaginary part: [nmod]
    Im_part = np.pi * iw4imag
      
    return Re_part, Im_part

def cal_real_part(iw4real, f2_target, f2s_integral):
    """Calculator of the real part of the term excluding 
    wave function.
    
    Parameters
    -----------
    iw4real : ndarray, shape=(ndiv, nmodes)
        Delta function at wi^2 - w(q)^2,
        where wi = f2s_integral[i] and w(q) is a frequency at q.
    f2_target : float
        Frequency at which Green's function is being calculated
    f2s_integral : array, float, shape=(ndiv,)
        Frequency points for the integration
    
    ndiv : integer
        # of points for the integration
    nmodes : integer
        # of modes

    Return
    ----------
    Re_part : array, float, shape=(nmodes,)
        Calculated value of each mode
    """
    _df2 = f2s_integral[1] - f2s_integral[0]
    
    ndiv = len(iw4real)
    nmodes = len(iw4real[0])
    
    Re_part = np.zeros(nmodes)
    for i2 in range(ndiv):
        Re_part += (
                _df2 * iw4real[i2] / 
                (- f2_target + f2s_integral[i2])
                )
        if abs(f2s_integral[i2] - f2_target) < 1e-3:
            print("ERROR")
            exit()
    
    return Re_part

def get_f2s4integral(f2_center, f2min, f2max, ndiv):
    """
    Note
    -------
    Get frequencies for the integration.
    This adjusts the frequency points to make f2_center locate on 
    the middle of adjacent two f^2 points.
    
    Parameters
    -----------
    f2min, f2max : float
        Minimum and maximum frequencies
    
    Return
    ---------
    f2s_integral : array, float, shape=(ndiv,)
        frequencies used for the integration
    """
    f2s_list = np.linspace(f2min, f2max, ndiv)
    
    if f2min <= f2_center and f2_center <= f2max:
        for i in range(ndiv - 1):
            if (f2s_list[i] <= f2_center and 
                    f2_center <= f2s_list[i+1]):
                f2s_list -= (
                        (f2s_list[i] + f2s_list[i+1]) * 0.5 - f2_center
                        )
     
    return f2s_list

def get_dphi(f2_target, Mmat0, Mmat1, Cmat0, Cmat1):
    """Calculator of the term representing the difference of
    force constants that appears for the calculation of T-matrix
    as T(w^2) = [1 - dphi * G0]^-1 * dphi.
    
    Note
    --------
    Unit of frequency and IFCs is changed to (rad/s) and 
    
    Parameters
    ------------
    f2_target : float
        Squared frequency [THz**2]
    Mmat0, Mmat1 : array, float, shape=(natoms,), unit=[AMU]
        Atomic masses of the pure and impurity systems
    Cmat0, Cmat1 : ndarray, float, shape=(natoms,natoms,3,3)
        Force constants of the pure crystal and 
        the impurity system. unit=[eV/A^-2] (for VASP)
    
    natoms : integer
        # of atoms

    Return
    ---------
    dphi : ndarray, float, shape=(3*natoms, 3*natoms)
    """
    natoms = len(Mmat0)
    dphi = np.zeros((3*natoms, 3*natoms))
    dp_ij = np.zeros((3, 3))
    for i1 in range(natoms):
        for i2 in range(i1,natoms):
            dp_ij[:,:] = (
                    (Cmat0[i1,i2,:,:] - Cmat1[i1,i2,:,:]) /
                    np.sqrt(Mmat0[i1] * Mmat0[i2])
                    ) * VaspToTHz**2    # Unit is changed to THz^-2
            if i1 == i2:
                dp_ij += (
                        np.eye(3) * f2_target * 
                        (Mmat1[i1] / Mmat0[i1] - 1.)
                        )
            
            dphi[3*i1:3*(i1+1), 3*i2:3*(i2+1)] = dp_ij
            if i1 != i2:
                dphi[3*i2:3*(i2+1), 3*i1:3*(i1+1)] = dp_ij.T
    
    return dphi

def get_UTU_braket(
        tmat, qpoint=None, f2=None, evec=None,
        pcell=None,
        scell=None,
        scaled_positions=None,
        s2pp=None):
    """Calculator of the scattering rate
    Parameters
    ------------
    tmat : ndarray, complex, shape=(3*nat_super, 3*nat_super)
        T-matrix
    qpoint : array, float, shape=(3,)
        Scaled wave vector
    f2 : float
        Squared frequency [THz**2]
    evec : array, float, shape=(3*nat_prim,)
    
    nat_{prim, super} : integer
        # of atoms in a primitive and super cells

    Return
    ---------
    rscat : float
        Scattering rate
    """
    nat_super = len(scaled_positions)
    wf_expand = get_expanded_wf(
            qpoint, evec, pcell, scell, scaled_positions, s2pp)
    rscat = np.dot(
            np.dot(np.conjugate(wf_expand), tmat), 
            np.matrix(wf_expand).T)
    
    return rscat[0,0]

def conv_tmat_L2s(LTmat, masses):
    """Conversion from t-matrix to T-matrix
    Note
    ------
    stmat = (1 - dv * g0)^-1 * dv
    LTmat = (1 - dphi * g0)^-1 * dphi
    
    Parameters
    -----------
    LTmat : complex, shape=(3*nat, 3*nat)
    masses : float, shape=(nat,)
        mass of each atom in the system with impurities
    n2s : integer, shape=(nat,)
        
    nat : integer
        # of atoms
    """
    nmat = len(LTmat)
    nat = len(masses)
    if nmat != 3*nat:
        print("Error. matrix shape is incompatible.")
        print(nmat, "!=", nat*3)
        import sys
        sys.exit()
    stmat = np.zeros_like(LTmat, dtype=np.complex)
    mlong = np.zeros(nmat)
    for iat in range(nat):
        mlong[3*iat:3*(iat+1)] = masses[iat] * np.ones(3)
    for iat in range(nat):
        stmat[3*iat:3*(iat+1)] = LTmat[3*iat:3*(iat+1)] / mlong
    return stmat

def get_dos_imp(stmat, g0, multiplicity=None):
    """Calculator of DOS of structures with an impurity
    Parameters
    ----------
    tmat : complex, shape=(3*nat, 3*nat)
        t-matrix, t = (1 - dv * g0)^-1
    g0 : complex, shape=(3*nat, 3*nat)
        Green's function of pure crystal
    n2ss : array, integer, shape=(natom,)
        atomic indices of supercell corresponding to the new cell
    """
    #if multiplicity is None:
    #    multiplicity = np.ones(len(g0))
    #nat = len(n2ss)
    #mul_long = np.zeros(3*nat)
    #for i in range(nat):
    #    mul_long[3*i:3*(i+1)] = multiplicity[i] * np.ones(3)
    
    # -- calculate Green's function of impurity system
    g1 = np.dot(np.eye(len(stmat)) + np.dot(g0, stmat), g0)
    dos_imp = get_dos_green(g1, multiplicity=multiplicity)
    return dos_imp

def conv_n2ss2weights(n2ss):
    """Get # of equivalent atoms
    Definitions
    ------------
    new cell : cell made for the calculation of T-matrix
    supercell : cell used for the calculation of IFCs

    Parameters
    ----------
    n2ss : array, integer, shape=(natom,)
        atomic indices of supercell corresponding to the new cell
    natom : integer
        # of atoms

    Return
    --------
    multi : array, float (but integer), shape=(natom,)
        # of equivalent atoms for each atom in the new cell
    """
    multi = np.zeros_like(n2ss)
    for i1 in range(len(n2ss)):
        for i2 in range(len(n2ss)):
            if n2ss[i2] == n2ss[i1]:
                multi[i1] += 1.
    return multi

