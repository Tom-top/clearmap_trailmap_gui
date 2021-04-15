"""
Convolution on a subset of points ina big array
"""

cimport cython
from cython.parallel import prange, parallel
import numpy as np
cimport numpy as np

@cython.boundscheck(False)
@cython.wraparound(False)
def cy_convolve3D(long[:, :, :] data, long[:, :, :] kernel, Py_ssize_t[:, :]  points, int processes = 4):
    """Convolves binary data with a specified kernel at specific points only
    
    Arguments:
        data (array): 3d binary array to convolve
        kernel (array): kernel to convolve
        points (array): list of points to convolve
        
    Returns:
        array: list of results of convolution at specified points
    """
    shape = np.asarray(data).shape
    cdef Py_ssize_t i, j, k, y, z, x, n
    cdef Py_ssize_t ks = kernel.shape[0]
    cdef Py_ssize_t ks2 = ks / 2;
    cdef Py_ssize_t npoints = points.shape[0]
    cdef long[:] result = np.zeros(npoints, dtype = int);
    with nogil, parallel(num_threads = processes):    
      for n in prange(npoints, schedule='static'):
        x = points[n, 0]; 
        y = points[n, 1]; 
        z = points[n, 2];
        #print n,x,y,z, data[x,y,z], ks2
        for k in range(ks):
          for i in range(ks):
            for j in range(ks):
              result[n] += data[x + i - ks2, y + j - ks2, z + k - ks2] * kernel[i, j, k]
    
    return np.array(result);

  