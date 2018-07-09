# -*- coding: utf-8 -*-
import numpy as np
from phonopy import Phonopy
from phonopy.phonon.tetrahedron_mesh import TetrahedronMesh

"""
Note:
    Get Tetrahedron meshes for frequency**2

Parameters:
    phonon (Phonopy)
    mesh[3]: (nq1, nq2, nq3)
"""
def get_tetrahedron_mesh_f2( phonon, mesh ):

    qpoints, weights, frequencies, eigenvecs = phonon.get_mesh()

    grid_address, ir_grid_points, grid_mapping_table = phonon.get_mesh_grid_info()

    primitive = phonon.get_primitive()

    f2s = frequencies**2 * np.sign(frequencies)
    thm = TetrahedronMesh(
            primitive, f2s, mesh,
            grid_address, grid_mapping_table,
            ir_grid_points)
    
    return thm, qpoints, f2s, eigenvecs
 
