#!python
#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False
"""
Cython Code for the LargeData module
"""
__author__    = 'Christoph Kirst <ckirst@rockefeller.edu>'
__license__   = 'MIT License <http://www.opensource.org/licenses/mit-license.php>'
__copyright__ = 'Copyright (c) 2017 by Christoph Kirst, The Rockefeller University, New York City'

import numpy as np
cimport numpy as np

from numpy cimport uint8_t, uint16_t, float32_t, double_t

cimport cython
from cython.parallel import parallel, prange, threadid

from libc.stdlib cimport abort, malloc, free

ctypedef fused index_t:
  Py_ssize_t

cdef extern from "stdio.h":
    int printf(char *format, ...) nogil

@cython.boundscheck(False)
@cython.wraparound(False)
cdef index_t neighbour26(index_t* indices, index_t n, index_t* offsets, index_t max_offset,
                         index_t center, index_t other) nogil:
  #center is index of position in indices, other is also an index in indices!                                                 
  #cdef Py_ssize_t n = indices.shape[0];
  #cdef Py_ssize_t no = offsets.shape[0];
  cdef Py_ssize_t no = 26;
 
  cdef index_t plo = max(-1, center - max_offset - 1);
  cdef index_t phi = min(n , center + max_offset + 1);
  
  cdef index_t c = indices[center];
  cdef index_t ilo = min(0, c - max_offset);
  cdef index_t ihi = c + max_offset;
  
  cdef index_t i, o, ic;
  
  #printf('nb: cc=%d (%d) plo=%d, phi=%d, ilo=%d, ihi=%d, n=%d, no=%d\n', center, c, plo, phi, ilo, ihi, n, no)
  #TODO: no need to check all offsets -> could check only the ones with are small/large enough
  for p in range(center - 1, plo, -1):
    i = indices[p];
    if i < ilo:  #no more neighbours to test
      break;
    ic = i - c;
    #printf('nd d: p=%d, i=%d, ic=%d\n', p,i,ic);
    for o in range(no):
      if (ic == offsets[o]):
        if (p != other):
          return p;
  
  for p in range(center + 1, phi,  1):
    i = indices[p];
    #printf('nd u: p=%d, i=%d\n', p,i);
    if i > ihi:  #no more neighbours to test
      break;
    ic = i - c;
    #printf('nd u: p=%d, i=%d, ic=%d\n', p,i,ic);
    for o in range(no):
      if (ic == offsets[o]):
        if (p != other):
          return p;
          
  return n+1;

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef np.ndarray[uint8_t, ndim=1] cleanBranchesIndex(index_t[:] points, Py_ssize_t[:] strides,
                                                     index_t length, 
                                                     index_t[:] startpoints, np.ndarray[uint8_t, ndim=1] stoppoints, int processes):
  """Remove branches from start points to end points if below certain length (index version)
  """

  # array sizes
  cdef Py_ssize_t ni = points.shape[0];
  cdef Py_ssize_t ns = startpoints.shape[0];
  
  #offsets of neighbours
  cdef Py_ssize_t[:] offsets = np.empty(26, dtype = np.int);
  cdef Py_ssize_t max_offsets = 0;
  cdef int dx,dy,dz;
  cdef Py_ssize_t i = 0;
  for dx in range(-1, 2):
    for dy in range(-1, 2):
      for dz in range(-1, 2):
        if (dx == 0 and dy == 0 and dz == 0):
          continue;
        offsets[i] = dx * strides[0] + dy * strides[1] + dz * strides[2];
        max_offsets = max(max_offsets, offsets[i]);
        #printf('o = %d %d %d = %d (max=%d)\n', dx,dy,dz,offsets[i], max_offsets);
        i += 1;  
  
  cdef Py_ssize_t l, ll, nb, pid, k;
  cdef index_t e, pre;
  cdef index_t* tmp;
  cdef index_t* local_points;
  cdef Py_ssize_t* local_offsets;
  
  cdef np.ndarray[uint8_t, ndim=1] keep = np.ones(ni, dtype = 'uint8');
  

  with nogil, parallel(num_threads = processes):
    
    tmp = <index_t *> malloc(sizeof(index_t) * length)
    if tmp == NULL:
        abort()
    
    local_points = <index_t *> malloc(sizeof(index_t) * ni);
    if local_points == NULL:
        abort()
    for k in range(ni):
      local_points[k] = points[k];
      
    local_offsets = <index_t *> malloc(sizeof(Py_ssize_t) * 26);
    if local_offsets == NULL:
        abort()
    for k in range(26):
      local_offsets[k] = offsets[k];
        
    
    for i in prange(ns, schedule = 'guided'):
      pid = threadid();
      e = startpoints[i];
      tmp[0] = ni + 1;
      #if pid == 5:
      if i % 100 == 0:
        printf('-----\ni=%d/%d pid=%d e=%d\n', i, ns, pid, e);
      for l in range(length):
        nb = neighbour26(local_points, ni, local_offsets, max_offsets, e, tmp[l]);
        #printf('nb=%d l=%d\n', nb, l);
        if nb == ni + 1: # no further neighbours -> skip
          break;
        tmp[l] = e;
        
        if stoppoints[nb] > 0: # found a branch to delete
          for ll in range(l+1):
            keep[tmp[ll]] = 0;
          break;
        
        # go to next point
        e = nb;
      
    free(tmp);
    free(local_points);
    free(local_offsets);
  
  return keep;
