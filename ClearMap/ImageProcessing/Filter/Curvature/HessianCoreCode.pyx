#!python
#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False

import numpy as np

cimport numpy as cnp

cimport HessianCoreCode

#from libc.stdlib cimport malloc, free

from libc.math cimport sqrt, atan2, cos, sin, pow

cdef extern from "stdio.h":
    int printf(char *format, ...) nogil

cdef inline double _abs(double a) nogil:
    return a if a >= 0 else -a;


cdef void core(void kernel(dtype_t_sink*, double, double, double, double) nogil,
               dtype_t_source[:, :, ::1] source,
               dtype_t_sink[:, :, :, ::1] sink,
               double par) except *:
    """Compute Hessian eigenvalues for each pixel and apply a measure defined by the kernel
    """

    # array sizes
    cdef Py_ssize_t rows = source.shape[0]
    cdef Py_ssize_t cols = source.shape[1]
    cdef Py_ssize_t plns = source.shape[2]   
    
    cdef Py_ssize_t max_img = rows * cols * plns
    
    #cdef Py_ssize_t odepth = sink.shape[3]
    
    # local variable types
    cdef Py_ssize_t r,c,p,rm,rp,pm,pp,cm,cp

    # Hessian matrix is represented as A = h[0,0], D = h[1,1], F = h[2,2], B = h[0,1]=h[1,0], C = h[0,2]=h[2,0], E = h[1,2]=h[2,1]
    #cdef double* hessian = <double*>malloc(6 * sizeof(double))

    cdef double A, B, C, D, E, F;
    cdef double temp, a, b, cc, d, third, rootThree, q, rr, discriminant
    cdef double innerSize, innerAngle, stSize, sPlusT, sAngle, firstPart, lastPart, bOver3A
    cdef double e1, e2, e3
    cdef double e1a, e2a, e3a
    cdef double e1s, e2s, e3s
    
    a = -1.0;
    third = 1.0/3;
    rootThree = 1.7320508075688772935;
    
    with nogil:
        for r in range(rows):
       
          if r > 0:
            rm = r - 1;
          else:
            rm = r;
          if r < rows - 1:
            rp = r + 1;
          else:
            rp = rows - 1;
          
          for c in range(cols):
            
            if c > 0:
              cm = c - 1;
            else:
              cm = c;
            if c < cols - 1:
              cp = c + 1;
            else:
              cp = cols - 1;
            
            for p in range(plns):
              if p > 0:
                pm = p - 1;
              else:
                pm = p;
              if p < plns - 1:
                pp = p + 1;
              else:
                pp = plns - 1;
            

              # create hessian
              temp = 2.0 * <double>source[r,c,p];
              
              A = source[rm,c, p ] - temp + source[rp,c, p ];
              D = source[r, cm,p ] - temp + source[r, cp,p ];
              F = source[r, c, pm] - temp + source[r, c ,pp];

              B = (<double>source[rp,cp,p ] - source[rm,cp,p ] - source[rp,cm,p ] + source[rm,cm,p ]) / 4.0;
              C = (<double>source[rp,c ,pp] - source[rm,c ,pp] - source[rp,c ,pm] + source[rm,c ,pm]) / 4.0;
              E = (<double>source[r ,cp,pp] - source[r ,cm,pp] - source[r ,cp,pm] + source[r ,cm,pm]) / 4.0;
 
              # calculate eigenvalues
              b = A + D + F;
              cc = B * B + C * C + E * E - A * D - A * F - D * F;
              d = A * D * F - A * E * E - B * B * F + 2 * B * C * E - C * C * D;
              
              q = (3*a*cc - b*b) / (9*a*a);
              rr = (9*a*b*cc - 27*a*a*d - 2*b*b*b) / (54*a*a*a);
              
              discriminant = q*q*q + rr*rr;

              if discriminant < 0:
                innerSize  = sqrt(rr*rr - discriminant);
                innerAngle = atan2(sqrt(-discriminant) , rr);
                stSize = pow(innerSize, third);
                sAngle = innerAngle / 3;
                sPlusT = 2 * stSize * cos(sAngle);
                firstPart = - (sPlusT / 2) - (b / 3*a);
                lastPart = - rootThree * stSize * sin(sAngle);
      
                e1 = ( sPlusT - (b / (3*a)) );
                e2 = ( firstPart + lastPart );
                e3 = ( firstPart - lastPart );
                
              elif discriminant == 0:
                
                if r >= 0:
                  sPlusT = 2 * pow(rr,third);
                else:                  
                  sPlusT = -2 * pow(-rr,third);
                bOver3A = b / (3 * a);
      
                e1 = ( sPlusT   - bOver3A );
                e2 = (-sPlusT/2 - bOver3A );
                e3 = e2;
              
              else:
                # discriminant > 0 should not happen for a symmetric matrix
                e1 = e2 = e3 = 0.0;
                
              # sort eigenvalues
              e1a = _abs(e1); e2a = _abs(e2); e3a = _abs(e3);
              
              if e1a <= e2a:
                if e2a <= e3a :
                  e1s = e1; e2s = e2; e3s = e3;
                else:
                  if e1a <= e3a:
                    e1s = e1; e2s = e3; e3s = e2;
                  else:
                    e1s = e3; e2s = e1; e3s = e2;
              else: #e2a < e1a
                if e1a <= e3a:
                  e1s = e2; e2s = e1; e3s = e3;
                else:
                  if e2a <= e3a:
                    e1s = e2; e2s = e3; e3s = e1;
                  else:
                    e1s = e3; e2s = e2; e3s = e1;
              
              # apply kernel
              kernel(&sink[r, c, p, 0], e1s, e2s, e3s, par);
      
