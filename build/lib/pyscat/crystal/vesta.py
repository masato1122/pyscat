# -*- coding: utf-8 -*-
import numpy as np

 
def mkvesta_IFCs(
        OFILE, 
        IFCs,
        scaled_positions,
        ibase, idir1, idir2, 
        IFC_max = None, IFC_min = None,
        TYPE="linear", cell=np.eye(3)
        ):
    """Make vesta file in whic atoms are colored 
    with following the magnitude of IFCs.
    
    Note:
        Visualize absolute values of IFC_{ibase, idir1; iat2, idir2},
        where ibase and iat2 (0 <= iat2 <= nat) are atom indices,
        idir1 and idir2 are direction, and nat is the number of
        atoms.
        It is preferable to use cubic cell.
        To visualize the output file, you may need "VESTA".
        http://jp-minerals.org/vesta/en/
    
    Parameters:
        OFILE   : output file name
        IFCs    : ndarray, shape=(nat, nat, 3, 3)
            Force constants
        scaled_positions: ndarray, shape=(nat, 3)
            Positions scaled by cell vectors,
            ranging from 0 to 1.
        ibase   : integer
            base atom index
        idir*   : integer
            direction (0.x, 1.y, 2.z)
    """
    ifcs_used = np.zeros_like(IFCs[ibase, :, idir1, idir2])
    ifcs_used[:] = IFCs[ibase, :, idir1, idir2]
    coords = np.dot(scaled_positions, cell)
    if np.sum(abs(cell - np.eye(3))) < 1e-5:
        size = 0.1
    else:
        size = 1.2
    
    if TYPE == "log":
        ifcs_used = np.log10(ifcs_used)
    
    # --- get the maximum and minimum IFCs
    if IFC_max is not None:
        ifcmax = IFC_max
    else:
        ifcmax = np.max(ifcs_used)

    if IFC_min is not None:
        ifcmin = IFC_min
    else:
        ifcmin = np.min(ifcs_used)
    
    ifcbase = (ifcmin + ifcmax) * 0.5
    
    nat = len(coords)
    ofs = open(OFILE, "w")
    ofs.write("#VESTA_FORMAT_VERSION 3.1.9\n")
    ofs.write("\n")
    ofs.write("MOLECULE\n")
    ofs.write("\n")
    ofs.write("TITLE\n")
    ofs.write("{:s}\n".format(OFILE))
    ofs.write("\n")
    ofs.write("STRUC\n")
    for ia in range(nat):
        ofs.write("{:d} H     A{:d} 1.00  ".format(ia+1, ia+1))
        for j in range(3):
            ofs.write("{:12.5f} ".format(coords[ia,j]))
        ofs.write("{:d}\n".format(ia+1))
        ofs.write("                 0.000000   0.000000   0.000000  0.00\n")
    ofs.write("  0 0 0 0 0 0 0\n")
    ofs.write("SITET\n")
    # -- {red, **, blue}
    col = np.zeros(3, dtype=int)
    for ia in range(nat):
        if ia == ibase:
            col[0] = 0
            col[1] = 255
            col[2] = 0
        else:
            if ifcs_used[ia] >= ifcbase:
                col[0] = 255
                ratio = (ifcs_used[ia] - ifcbase) / (ifcmax - ifcbase)
                if ratio > 1.0:
                    ratio = 1.0
                col[1] = col[2] = int(255 * (1. - ratio))
            else:
                col[2] = 255
                ratio = (ifcs_used[ia] - ifcbase) / (ifcmin - ifcbase)
                if ratio > 1.0:
                    ratio = 1.0
                col[0] = col[1] = int(255 * (1. - ratio))
        
        ofs.write("{:2d}     A{:d}  {:f} {:3d} {:3d} {:3d}"
                .format(ia+1, ia+1, size, col[0], col[1], col[2]))
        ofs.write("      255 204 204 204  0\n")
    ofs.write("  0 0 0 0 0 0 0\n")
    ofs.close()
    print(" Output:", OFILE, " Max:{:10.2e}, min:{:10.2e}"
            .format(ifcmax, ifcmin))


