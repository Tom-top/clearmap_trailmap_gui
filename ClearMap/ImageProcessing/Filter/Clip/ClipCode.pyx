#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False

cimport cython
from cython.parallel import prange, parallel

cimport numpy as np

ctypedef fused source_t:
  np.int32_t
  np.int64_t
  np.uint8_t
  np.uint16_t
  np.uint32_t
  np.uint64_t
  np.float32_t
  np.double_t
  
ctypedef fused sink_t:
  np.int32_t
  np.int64_t
  np.uint8_t
  np.uint16_t
  np.uint32_t
  np.uint64_t
  np.float32_t
  np.double_t
  
ctypedef Py_ssize_t index_t


cdef extern from "stdio.h":
    int printf(char *format, ...) nogil

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef void clip(source_t[:,:,:] source, sink_t[:, :, :] sink,
                double clipmin, double clipmax, double norm, int processes) nogil:
    """Clip image and normalize"""

    # array sizes
    cdef index_t rows = source.shape[0]
    cdef index_t cols = source.shape[1]
    cdef index_t plns = source.shape[2]   
  
    # local variable types
    cdef index_t r,c,p
    cdef double temp
    
    with nogil, parallel(num_threads = processes):   
      for r in prange(rows, schedule='guided'):
        for c in range(cols):
          for p in range(plns):
            temp = <double>source[r,c,p];
            if temp < clipmin:
              temp = clipmin;
            if temp > clipmax:
              temp = clipmax;
            sink[r, c, p] = <sink_t>(norm * (temp - clipmin)/(clipmax-clipmin));
            