import numpy as np
from numba import jit

@jit(nopython=True)
def fast_cross(v1, v2):
    return np.array([v1[1]*v2[2] - v1[2]*v2[1],
                     v1[2]*v2[0] - v1[0]*v2[2],
                     v1[0]*v2[1] - v1[1]*v2[0]])
