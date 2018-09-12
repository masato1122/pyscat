import numpy as np
from phonopy.structure.cells import get_smallest_vectors

def get_corresp_sites( prim ):
    """
    Input:
        prim: Primitive (see Phonopy)
            primitive cell (PC)
        p2p_map[Ns], s2p_map[Ns]:
            See descriptions in Phonopy
        Ns : number of atoms
    
    Return:
        coressponding sites in PC to the atoms in SC 
    
    """
    p2p_map = prim.get_primitive_to_primitive_map()
    s2p_map = prim.get_supercell_to_primitive_map()
    s2pp = np.array([p2p_map[i] for i in s2p_map], dtype=int)
    return s2pp

def offset_cell_position( cell, pos_sc, base_posi ):
    """Set atoms around an arbitrary position

    Parameters
    --------------
    cell        : ndarray((3,3,))
        cell vectors
    pos_sc      : ndarray((natom,3))
        atom positions scaled by the cell vectors
    base_posi   : array(3)
        the arbitrary position around which atoms are aranged
    
    Variables
    -------------
    svecs : ndarray, complex, shape=(nat_sc, 1, 27, 3)
        smallest vectors
    multi : multiplicity, shape(Na,1)

    Returns
    -----------
    ids_orig : array, integer, shape=(nat_nc)
        List of original atom IDs
    nequives : array, integer, shape=(nat_nc)
        # of equilibrium atoms
    scaled_disps : ndarray, float, shape=(nat_nc, 3)
        scaled displacement that should be integer like [1,0,0]
        
    """
    svecs, multi = get_smallest_vectors(cell, pos_sc, np.array([base_posi]))
    
    nat_sc = len(svecs)
    nat_nc = np.sum(multi)
    ids_orig = []
    scaled_disps = np.zeros((nat_nc,3))
    pos_nc = np.zeros((nat_nc,3))
    count = 0
    for ia_sc in range(nat_sc):
        nmulti = multi[ia_sc,0]
        for im in range(nmulti):
            posi = base_posi + svecs[ia_sc, 0, im]
            ids_orig.append(ia_sc)
            scaled_disps[count,:] = posi - pos_sc[ia_sc,:]
            count += 1
    
    nequivs = []
    for ia in range(count):
        nequives = multi[ids_orig[ia]]
    
    #pos_nc = pos_sc[ids_orig] + scaled_disps
    #__mkxyzfile( "check.xyz", cell, ids_orig, pos_sc, scaled_disps )
    
    return ids_orig, nequives, scaled_disps

def get_IFCs(ph_orig, ids_orig, scaled_disps):
    """Get interatomic force constants of an arbitrary atomic pair
    in the new cell.
    
    Note
    ---------
    Vectors will be scaled by supercell.

    Paraemters
    -----------
    ph_orig : Phonopy
        See descriptions in Phonopy for details
    ids_orig : array, integer, shape=(nat_nc,)
        List of atomic IDs in the original cell corresponding
        to the atoms in the new cell.
    scaled_disps : narray, float, shape=(nat_nc,3)
        Displacement of atoms from the original position
        to a position in the new cell.
        Unit is scaled by supercell
    nat_nc : # of atoms in a new cell
        = len(ids_orig)
    
    Return
    -------
    IFCs : ndarray, float, shape=(nat_nc, nat_nc, 3, 3)
        force constants
    
    """
    nat_nc = len(ids_orig)
    IFCs = np.zeros((((nat_nc, nat_nc, 3, 3))))
    
    pcell = ph_orig.get_primitive()
    scell = ph_orig.get_supercell()
    s2pp = get_corresp_sites(pcell)

    # -- transform matrix
    M_p2s = np.dot( pcell.cell, np.linalg.inv(scell.cell) )
    
    # -- svecs.shape = (nat_sc, nat_pc, 27, 3)
    # -- multi.shape = (nat_sc, nat_pc)
    svecs, multi = pcell.get_smallest_vectors()
    
    # -- convert to vectors scaled by supercell
    svecs = np.dot( svecs, M_p2s )
    
    # -----------------------------------
    # atomic positions in each cell
    # **_pos.shape = ((nat_**,3))
    # pc, sc, and nc are
    # primitive, super, and new cells
    # nat_nc >= nat_sc >= nat_pc
    # -----------------------------------
    pos_pc = np.dot(pcell.get_scaled_positions(), M_p2s)
    pos_sc = scell.get_scaled_positions()
    pos_nc = pos_sc[ids_orig] + scaled_disps

    for i1_nc in range(nat_nc):
        
        # ---------------------------------------------------
        # i1_{sc,pc}:
        # atom index in the SC or PC corresponding to i1_nc            
        # ---------------------------------------------------
        i1_sc = ids_orig[i1_nc]
        i1_pc = s2pp[i1_sc]
        
        for i2_nc in range(i1_nc, nat_nc):
            
            R12 = pos_nc[i2_nc] - pos_nc[i1_nc]
            R2_search = pos_pc[i1_pc] + R12
            
            i2_sc_trans, nmulti_i2 = _get_atom_in_SC(
                    R2_search, svecs[:,i1_pc,:,:], multi[:,i1_pc])
            
            if i2_sc_trans is None or nmulti_i2 is None:
                IFCs[i1_nc,i2_nc,:,:] = np.zeros((3,3))
                if i1_nc != i2_nc:
                    IFCs[i2_nc,i1_nc,:,:] = np.zeros((3,3))
            else:
                IFCs[i1_nc,i2_nc,:,:] = (
                        ph_orig.force_constants[i1_pc,i2_sc_trans,:,:] /
                        float(nmulti_i2))
                
                # -- transpose a matrix for the same pair
                if i1_nc != i2_nc:
                    IFCs[i2_nc,i1_nc] = IFCs[i1_nc,i2_nc].T
    
    return IFCs

def _get_atom_in_SC( R2, R2_cands, multi, eta=1e-5 ):
    """ Get the atom index in the supercell
        corresponding to "isite1-position + R12"
    Note:
        Find "i2" that satisfies "R2_cands[i2, j] == R2",
        (0 <= j <= multi[i2])
    Parameters:
        R2       : array, shape=(3,)
        R2_cands : narray, shape=(nat_super, 27, 3)
            in which R2 is searched.
        multi : multiplicity, shape=(nat_super,),
            where
            nat_prim and nat_super are # of atom in PC and SC.
    """
    I2_GET, NMUL_GET = None, None
    for i2 in range(len(multi)):
        nmulti = multi[i2]
        for im in range(nmulti):
            delta = R2 - R2_cands[i2,im]
            if np.linalg.norm(delta) < eta:
                I2_GET = i2
                NMUL_GET = nmulti
                return I2_GET, NMUL_GET
    return I2_GET, NMUL_GET


def _check_existence( pos, pos_list, eta=1e-5 ):
    for ia in range(len(pos_list)):
        if np.linalg.norm( pos_list[ia] - pos ) < eta:
            return True
    return None

def __mkxyzfile( OFILE, cell, id_list, scaled_posi, scaled_disps ):
    """Make xyz file to check the structure
    """
    new_posi = scaled_posi[id_list] + scaled_disps
    
    ofs = open(OFILE, "w")
    ofs.write("{:2d}\n".format(len(id_list)))
    ofs.write("\n")
    for ia in range(len(id_list)):
        ofs.write("H ")
        ri = np.dot(new_posi[ia], cell)
        for j in range(3):
            ofs.write("{:10.7f} ".format( ri[j] ))
        ofs.write("\n")
    print("Output:", OFILE)

def _check_IFCs_symmetry( IFCs, pos_nc, eta=0.01 ):
    nat = len(IFCs)
    for i1 in range(1,nat):
        for j1 in range(3):
            for i2 in range(i1+1,nat):
                L12 = np.linalg.norm( pos_nc[i1] - pos_nc[i2] )
                for j2 in range(3):
                    if abs(IFCs[i1,i2,j1,j2] - IFCs[i2,i1,j2,j1]) > eta:

                        print("L12 = ", L12)
                        print("Error", i1, j1, i2, j2)
                        for i in range(3):
                            for j in range(3):
                                print("{:15.7f} ".format(IFCs[i1,i2,i,j]), end="")
                            print("")
                        print("")
                        for i in range(3):
                            for j in range(3):
                                print("{:15.7f} ".format(IFCs[i2,i1,i,j]), end="")
                            print("")
                        print("")

                        import sys
                        sys.exit()
    return True

