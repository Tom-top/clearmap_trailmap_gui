#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False

#cimport numpy as cnp

from libc.math cimport sqrt

from .HessianCoreCode cimport dtype_t_source, dtype_t_sink, core

cdef inline void kernel_tubeness(dtype_t_sink* sink,
                                 double e1s, double e2s, double e3s, 
                                 double par) nogil:
  
  if e2s < 0 and e3s < 0:
    sink[0] = <dtype_t_sink>sqrt(e2s * e3s);
  else:
    sink[0] = 0;


def tubeness(dtype_t_source[:, :, ::1] source,
             dtype_t_sink[:, :, :, ::1] sink,
             double par):
  
  core(kernel_tubeness[dtype_t_sink], source, sink, par);



cdef inline void kernel_tubeness_threshold(dtype_t_sink* sink,
                                           double e1s, double e2s, double e3s, 
                                           double par) nogil:
  if e2s < 0 and e3s < 0:
    if sqrt(e2s * e3s) > par:
      sink[0] = 1;
    else:
      sink[0] = 0;
  else:
    sink[0] = 0;
    
    
def tubeness_threshold(dtype_t_source[:, :, ::1] source,
                       dtype_t_sink[:, :, :, ::1] sink,
                       double par):
  
  core(kernel_tubeness_threshold[dtype_t_sink], source, sink, par);
