# -*- coding: utf-8 -*-
import numpy as np
from phonopy import Phonopy
from phonopy.interface.vasp import read_vasp_from_strings, read_vasp
from phonopy.file_IO import (parse_FORCE_SETS, parse_BORN)
from pyscat.utils.idata import Idata


def get_grid(
        idat,
        is_mesh_symmetry=True,
        is_eigenvectors=True):
    """Get phonon properties on grid
    Parameters
    ------------
    idat : Idata

    Return
    -------
    qs : ndarray, float, shape=(nqpoints, 3)
        q-points
    weights : ndarray, float, shepe=(nqpoints)
        Multiplicity flowing the symmetry.
        Type is float, but values are integer.
    f2s : ndarray, float, shape=(nqpoints, nmodes)
        Squared frequencies [THz**2]
    evecs : ndarray, float, shape=(nqpoints, 3*nsites, nmodes) 
        Eigenvectors.
    """
    phonon = get_phonopy(idat)
    phonon.set_mesh(
            idat.nmesh,
            is_mesh_symmetry=is_mesh_symmetry,
            is_eigenvectors=is_eigenvectors)
    qs, weights, frequencies, evecs = phonon.get_mesh()
    f2s = frequencies**2 * np.sign(frequencies)
    return qs, weights, f2s, evecs


def get_phonopy(idat):
    """Make and return a Phonpy class
    containing infomration of structure and IFCs

    Parameters
    ------------
    idat : Idata
        Containing basic file names and vectors 
        See the description of Idata class for details.
    
    Return
    --------
    phonon : Phonpy
    """
    phonon = set_phonopy_structure(idat)
    add_ifcs2phonopy(phonon, idat)
    return phonon
 
def set_phonopy_structure(idat):
    """
    Parameters
    -------------
    idat : Idata
        Containing basic file names and vectors 
        See the description of Idata class for details.
    
    Return
    --------
    phonon : Phonopy
        Structure of POSCAR file
    """
    unitcell = read_vasp(idat.fposcar)
    phonon = Phonopy(unitcell, idat.ncells, primitive_matrix=idat.primat)
    return phonon

def add_ifcs2phonopy(phonon, idat):
    """Addition of IFCs data to Phonopy calss
    
    Parameters
    -----------
    phonon : Phonpy
    idat : Idata

    Return
    ---------
    phonon : Phonopy
        There is a return, but it is not necessary because 
        the contents of phonon is added.
    """
    # -- Read force sets and calculate IFCs
    force_sets = parse_FORCE_SETS(filename=idat.fforce)
    phonon.set_displacement_dataset(force_sets)
    phonon.produce_force_constants()
    phonon.symmetrize_force_constants()
    
    # -- Set Born effective charge if provided.
    if idat.fborn is not None:
        primitive = phonon.get_primitive()
        nac_params = parse_BORN(primitive, filename=idat.fborn)
        phonon.set_nac_params(nac_params)
    return phonon


