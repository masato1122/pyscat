# -*- coding: utf-8 -*-
import numpy as np

def check_cells(cell1, cell2, flag=True, tolerance=1e-5):
    """Check of cell size
    cell{1,2}: ndmap, float, shape=(3,3)
        cell vectors
    flag : True or False
    """
    if flag is not True:
        print(" Cell size check is OFF.")
        return False
    
    diff = cell1 - cell2
    FLAG = 0
    for i in range(3):
        for j in range(3):
            if abs(diff[i][j]) > tolerance:
                print("Error: the difference of the cell sizes is too large.")
                print(" Pure system: ", end="")
                print(cell1)
                print(" w/ impurity: ", end="")
                print(cell2)
                print("")
                import sys
                sys.exit()
    return True

def check_number_of_atoms(nat1, nat2, flag=True):
    """Check number of atoms
    Parameters
    -----------
    nat{1,2} : integer
        # of atoms in the cell 1 and 2
    """
    if flag is not True:
        print(" natom check is OFF.")
        return False
    
    if nat1 != nat2:
        print("Error: number of atoms is different in the pure system and "
                "in the system with impurities.")
        import sys
        sys.exit()
    return True

def check_positions(pos1, pos2, cell=np.eye(3), flag=True, tolerance=1.0):
    """Check positions
    Parameters
    -----------
    pos1, pos2 : ndarray, float, shape=(nat,3)
        scaled positions
    cell : ndarray, float, shape=(3,3)
        scaled positions
    """
    if flag is not True:
        print(" Position check is OFF.")
        return False
    diff = pos1 - pos2
    
    for ia in range(len(pos1)):
        print(pos1[ia], pos2[ia])

    exit()
    diff = np.linalg.norm(np.dot(diff, cell), axis=1)
    print(diff)
    exit()
    FLAG = 0
    for ia in range(len(diff)):
        for j in range(3):
            if abs(diff[ia,j]) > tolerance:
                print(" Error: difference of an atomic position of the pure "
                        "crystal is too large.")
                print(" {:d}-th atom has been moved {:f} A.".format(ia,
                    diff[ia,j]))
                import sys
                sys.exit()
    if FLAG == 0:
        return True
    else:
        return None

