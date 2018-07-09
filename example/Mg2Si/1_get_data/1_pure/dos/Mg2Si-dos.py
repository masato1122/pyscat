import numpy as np
from phonopy import Phonopy
from phonopy.interface.vasp import read_vasp_from_strings, read_vasp
from phonopy.file_IO import (parse_FORCE_SETS, parse_BORN)
from phonopy.phonon.tetrahedron_mesh import TetrahedronMesh

import matplotlib.pyplot as plt

def get_phonon(mesh):
    """Initialize phonopy for phonos on mesh grid sampling

    For initializing phonopy, more details are found in NaCl.py.
    
    """
     
    unitcell = read_vasp("POSCAR")
    
    phonon = Phonopy(unitcell,
                     [[2, 0, 0],
                      [0, 2, 0],
                      [0, 0, 2]],
                     primitive_matrix=[[0, 0.5, 0.5],
                                       [0.5, 0, 0.5],
                                       [0.5, 0.5, 0]])
    
    force_sets = parse_FORCE_SETS()      #--- read "FORCE_SETS"
    
    phonon.set_displacement_dataset( force_sets )
    
    phonon.produce_force_constants()
    
    phonon.symmetrize_force_constants() 
    
    primitive = phonon.get_primitive()
    
    #nac_params = parse_BORN(primitive)
    #phonon.set_nac_params(nac_params)
    
    #--- calculation of eigenvalues & eigenvectors!!
    phonon.set_mesh(mesh, is_eigenvectors=True)
    
    return phonon

def get_tetrahedron_mesh(mesh):
    
    phonon = get_phonon(mesh)

    #--- 1.qpoints, 2.weights, 3.frequencies, and 4.eigenvectors (= None)
    #--- eigenvecs[iq, iat, iw] ??
    qpoints, weights, frequencies, eigenvecs = phonon.get_mesh()
    
    (grid_address, 
            ir_grid_points, 
            grid_mapping_table) = phonon.get_mesh_grid_info()
    
    primitive = phonon.get_primitive()
    
    #--- 
    thm = TetrahedronMesh(primitive,
                          frequencies,
                          mesh,
                          grid_address,
                          grid_mapping_table,
                          ir_grid_points)
    
    fmin, fmax = frequencies.min(), frequencies.max()
    frequency_points = np.linspace(fmin, fmax, 200)
    thm.set(value='I', frequency_points=frequency_points)
    return thm, weights, frequencies, eigenvecs

def get_dos(mesh):
    
    thm, weights, frequencies, eigenvecs = get_tetrahedron_mesh(mesh)
    dos = np.zeros_like(thm.get_frequency_points())
    
    #--- calculation of DOS
    for i, iw in enumerate(thm):
        dos += np.sum(iw * weights[i], axis=1)

    return thm.get_frequency_points(), dos

def draw(f, d):
    for i in range(len(f)):
        print("{:10.4f} {:10.4f}".format(f[i], d[i]))
    plt.style.use('seaborn-darkgrid')
    plt.plot(f, d, '-')
    plt.legend()
    plt.xlabel("Frequency (THz)")
    plt.ylabel("DOS (a.u.)")
    #plt.show()
    plt.savefig("dos.png")
    
def main():
    mesh = [11, 11, 11]
    #mesh = [5, 5, 5]
    frequency_points, dos = get_dos(mesh)
    draw(frequency_points, dos)

if __name__ == "__main__":
    main()
