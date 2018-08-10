# -*- coding: utf-8 -*-
import numpy as np

def get_frequency2( frequency=None, f2=None ):
    """
    Note:
        Set the frequency at which Green's function is calculated
        from a given frequency with the unit of [THz] or [THz**2].
        \"frequency\" or \"f2\" have to be given.

    """
    if f2 is not None:
        f2_target = f2
    else:
        if frequency is not None:
            f2_target = frequency**2 * np.sign(frequency)
        else:
            import sys
            print("Error: \"frequency\" or \"f2\", is not given.")
            sys.exit()
    
    return f2_target

