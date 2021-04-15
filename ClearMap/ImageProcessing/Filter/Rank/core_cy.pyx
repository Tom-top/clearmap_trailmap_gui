#!python
#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False

import numpy as np

cimport numpy as cnp
from libc.stdlib cimport malloc, free

cimport core_cy

cdef extern from "stdio.h":
    int printf(char *format, ...) nogil


cdef inline dtype_t _max(dtype_t a, dtype_t b) nogil:
    return a if a >= b else b


cdef inline dtype_t _min(dtype_t a, dtype_t b) nogil:
    return a if a <= b else b


cdef inline void histogram_increment(Py_ssize_t* histo, double* pop,
                                     dtype_t value) nogil:
    histo[value] += 1
    pop[0] += 1


cdef inline void histogram_decrement(Py_ssize_t* histo, double* pop,
                                     dtype_t value) nogil:
    histo[value] -= 1
    pop[0] -= 1


cdef inline char is_in_mask(Py_ssize_t rows, Py_ssize_t cols, Py_ssize_t plns, 
                            Py_ssize_t r, Py_ssize_t c, Py_ssize_t p,
                            char* mask) nogil:
    """Check whether given coordinate is within image and mask is true."""
    if r < 0 or r > rows - 1 or c < 0 or c > cols - 1 or p < 0 or p > plns - 1:
        return 0
    else:
        if mask[p + c * plns + r * cols * plns]:
            return 1
        else:
            return 0

cdef inline char is_in_image(Py_ssize_t rows, Py_ssize_t cols, Py_ssize_t plns, 
                            Py_ssize_t r, Py_ssize_t c, Py_ssize_t p) nogil:
    """Check whether given coordinate is within image and mask is true."""
    if r < 0 or r > rows - 1 or c < 0 or c > cols - 1 or p < 0 or p > plns - 1:
      return 0
    else:
      return 1



cdef void _core(void kernel(dtype_t_out*, Py_ssize_t, Py_ssize_t*, double,
                            dtype_t, Py_ssize_t, Py_ssize_t, double,
                            double, Py_ssize_t, Py_ssize_t) nogil,
                dtype_t[:, :, ::1] image,
                char[:, :, ::1] selem,
                char[:, :, ::1] mask,
                dtype_t_out[:, :, :, ::1] out,
                signed char shift_x, signed char shift_y, signed char shift_z,
                double p0, double p1,
                Py_ssize_t s0, Py_ssize_t s1,
                Py_ssize_t max_bin) except *:
    """Compute histogram for each pixel neighborhood, apply kernel function and
    use kernel function return value for output image.
    """

    cdef Py_ssize_t rows = image.shape[0]
    cdef Py_ssize_t cols = image.shape[1]
    cdef Py_ssize_t plns = image.shape[2]   
    
    cdef Py_ssize_t max_img = rows * cols * plns
    
    cdef Py_ssize_t srows = selem.shape[0]
    cdef Py_ssize_t scols = selem.shape[1]
    cdef Py_ssize_t splns = selem.shape[2]
    
    cdef Py_ssize_t max_se  = srows * scols * splns
    
    cdef Py_ssize_t odepth = out.shape[3]

    cdef Py_ssize_t centre_r = <Py_ssize_t>(selem.shape[0] / 2) + shift_y
    cdef Py_ssize_t centre_c = <Py_ssize_t>(selem.shape[1] / 2) + shift_x
    cdef Py_ssize_t centre_p = <Py_ssize_t>(selem.shape[2] / 2) + shift_z

    # check that structuring element center is inside the element bounding box
    assert centre_r >= 0
    assert centre_c >= 0
    assert centre_p >= 0
    
    assert centre_r < srows
    assert centre_c < scols
    assert centre_p < splns
    
    #print centre_r, centre_c, centre_p
    #print max_se


    # add 1 to ensure maximum value is included in histogram -> range(max_bin)
    max_bin += 1

    cdef Py_ssize_t mid_bin = max_bin / 2

    # define pointers to the data
    #cdef char* mask_data = &mask[0, 0, 0]

    # define local variable types
    cdef Py_ssize_t r, c, p, rr, cc, pp, s, value, local_max, i, dir_r, dir_c

    # number of pixels actually inside the neighborhood (double)
    cdef double pop = 0

    # the current local histogram distribution
    cdef Py_ssize_t* histo = <Py_ssize_t*>malloc(max_bin * sizeof(Py_ssize_t))
    for i in range(max_bin):
        histo[i] = 0

    # these lists contain the relative pixel row and column for each of the 4
    # attack borders east, west, north and south e.g. se_e_r lists the rows of
    # the east structuring element border

    # number of element in each attack border
    cdef Py_ssize_t num_se_n, num_se_s, num_se_e, num_se_w, num_se_d, num_se_u
    num_se_n = num_se_s = num_se_e = num_se_w = num_se_d = num_se_u = 0

    cdef Py_ssize_t* se_e_r = <Py_ssize_t*>malloc(max_se * sizeof(Py_ssize_t))
    cdef Py_ssize_t* se_e_c = <Py_ssize_t*>malloc(max_se * sizeof(Py_ssize_t))
    cdef Py_ssize_t* se_e_p = <Py_ssize_t*>malloc(max_se * sizeof(Py_ssize_t))
    
    cdef Py_ssize_t* se_w_r = <Py_ssize_t*>malloc(max_se * sizeof(Py_ssize_t))
    cdef Py_ssize_t* se_w_c = <Py_ssize_t*>malloc(max_se * sizeof(Py_ssize_t))
    cdef Py_ssize_t* se_w_p = <Py_ssize_t*>malloc(max_se * sizeof(Py_ssize_t)) 

    
    cdef Py_ssize_t* se_n_r = <Py_ssize_t*>malloc(max_se * sizeof(Py_ssize_t))
    cdef Py_ssize_t* se_n_c = <Py_ssize_t*>malloc(max_se * sizeof(Py_ssize_t))
    cdef Py_ssize_t* se_n_p = <Py_ssize_t*>malloc(max_se * sizeof(Py_ssize_t))
    
    cdef Py_ssize_t* se_s_r = <Py_ssize_t*>malloc(max_se * sizeof(Py_ssize_t))
    cdef Py_ssize_t* se_s_c = <Py_ssize_t*>malloc(max_se * sizeof(Py_ssize_t))
    cdef Py_ssize_t* se_s_p = <Py_ssize_t*>malloc(max_se * sizeof(Py_ssize_t))
    
    
    cdef Py_ssize_t* se_u_r = <Py_ssize_t*>malloc(max_se * sizeof(Py_ssize_t))
    cdef Py_ssize_t* se_u_c = <Py_ssize_t*>malloc(max_se * sizeof(Py_ssize_t))
    cdef Py_ssize_t* se_u_p = <Py_ssize_t*>malloc(max_se * sizeof(Py_ssize_t))
    
    cdef Py_ssize_t* se_d_r = <Py_ssize_t*>malloc(max_se * sizeof(Py_ssize_t))
    cdef Py_ssize_t* se_d_c = <Py_ssize_t*>malloc(max_se * sizeof(Py_ssize_t))
    cdef Py_ssize_t* se_d_p = <Py_ssize_t*>malloc(max_se * sizeof(Py_ssize_t))

    # build attack and release borders by using difference along axis
    t = np.concatenate((selem, np.zeros((1, selem.shape[1], selem.shape[2]))), axis = 0)
    cdef unsigned char[:, :, :] t_e = (np.diff(t, axis=0) < 0).view(np.uint8)

    t = np.concatenate((np.zeros((1, selem.shape[1], selem.shape[2])), selem), axis = 0)
    cdef unsigned char[:, :, :] t_w = (np.diff(t, axis=0) > 0).view(np.uint8)
    
    
    t = np.concatenate((selem, np.zeros((selem.shape[0], 1, selem.shape[2]))), axis = 1)
    cdef unsigned char[:, :, :] t_n = (np.diff(t, axis=1) < 0).view(np.uint8)

    t = np.concatenate((np.zeros((selem.shape[0], 1, selem.shape[2])), selem), axis = 1)
    cdef unsigned char[:, :, :] t_s = (np.diff(t, axis=1) > 0).view(np.uint8)
    
    
    t = np.concatenate((selem, np.zeros((selem.shape[0], selem.shape[1], 1))), axis = 2)
    cdef unsigned char[:, :, :] t_u = (np.diff(t, axis=2) < 0).view(np.uint8)

    t = np.concatenate((np.zeros((selem.shape[0], selem.shape[1], 1)), selem), axis = 2)
    cdef unsigned char[:, :, :] t_d = (np.diff(t, axis=2) > 0).view(np.uint8)
  
    
    cdef Py_ssize_t mode;
    
    with nogil:

        #define attack borders
        for p in range(splns):
          for c in range(scols):
            for r in range(srows):
              if t_e[r, c, p]:
                  se_e_r[num_se_e] = r - centre_r
                  se_e_c[num_se_e] = c - centre_c
                  se_e_p[num_se_e] = p - centre_p         
                  num_se_e += 1
              if t_w[r, c, p]:
                  se_w_r[num_se_w] = r - centre_r
                  se_w_c[num_se_w] = c - centre_c
                  se_w_p[num_se_w] = p - centre_p
                  num_se_w += 1
              
              if t_n[r, c, p]:
                  se_n_r[num_se_n] = r - centre_r
                  se_n_c[num_se_n] = c - centre_c
                  se_n_p[num_se_n] = p - centre_p  
                  num_se_n += 1
              if t_s[r, c, p]:
                  se_s_r[num_se_s] = r - centre_r
                  se_s_c[num_se_s] = c - centre_c
                  se_s_p[num_se_s] = p - centre_p 
                  num_se_s += 1
                  
              if t_u[r, c, p]:
                  se_u_r[num_se_u] = r - centre_r
                  se_u_c[num_se_u] = c - centre_c
                  se_u_p[num_se_u] = p - centre_p  
                  num_se_u += 1 
              if t_d[r, c, p]:
                  se_d_r[num_se_d] = r - centre_r
                  se_d_c[num_se_d] = c - centre_c
                  se_d_p[num_se_d] = p - centre_p  
                  num_se_d += 1

                  

        # create historgram
        for r in range(srows):
          for c in range(scols):
            for p in range(splns):
              rr = r - centre_r
              cc = c - centre_c
              pp = p - centre_p
              if selem[r, c, p]:
                  #if is_in_mask(rows, cols, plns, rr, cc, pp, mask_data):
                  if is_in_image(rows, cols, plns, rr, cc, pp):
                      histogram_increment(histo, &pop, image[rr, cc, pp])
                      #printf('hist adding (%d,%d,%d) [%d]\n', rr, cc, pp, image[rr,cc,pp]);
        r = 0
        c = 0
        p = 0;
        kernel(&out[r, c, p, 0], odepth, histo, pop, image[r, c, p], max_bin, mid_bin,
               p0, p1, s0, s1)

        # main loop
        dir_r = 1;
        dir_c = 1;
        mode = 0; # modes 0 r+=1, 1 r-=1, 2 c+=1, 3 c-=1, 4 p+=1
        
        while True:
          
          #printf('-----\n');          
          
          if dir_r == 1 and r < rows - 1:
            r += 1;
            mode = 0;
          elif dir_r == -1 and r > 0:
            r -= 1;
            mode = 1;
          else:
            if dir_c == 1 and c < cols - 1:
              c += 1;
              dir_r *= -1;
              mode = 2;
            elif dir_c == -1 and c > 0:
              c -= 1;
              dir_r *= -1;
              mode = 3;
            else:
              p += 1
              if p == plns:
                break;
              dir_r *= -1;
              dir_c *= -1;
              mode = 4;

          if mode == 0:
            for s in range(num_se_e):
                  rr = r + se_e_r[s]
                  cc = c + se_e_c[s]
                  pp = p + se_e_p[s]
                  #if is_in_mask(rows, cols, plns, rr, cc, pp, mask_data):
                  if is_in_image(rows, cols, plns, rr, cc, pp):
                      histogram_increment(histo, &pop, image[rr, cc, pp])
                      #printf('mode %d adding (%d,%d,%d) [%d]\n', mode, rr, cc, pp, image[rr,cc,pp]);
    
            for s in range(num_se_w):
                  rr = r + se_w_r[s] - 1
                  cc = c + se_w_c[s] 
                  pp = p + se_w_p[s]
                  #if is_in_mask(rows, cols, plns, rr, cc, pp, mask_data):
                  if is_in_image(rows, cols, plns, rr, cc, pp):
                      histogram_decrement(histo, &pop, image[rr, cc, pp])
                      #printf('mode %d removing (%d,%d,%d) [%d]\n', mode, rr, cc, pp, image[rr,cc,pp]);
    
            kernel(&out[r, c, p, 0], odepth, histo, pop, image[r, c, p], max_bin,
                   mid_bin, p0, p1, s0, s1)
                     
          
          elif mode == 1:
              for s in range(num_se_w):
                  rr = r + se_w_r[s]
                  cc = c + se_w_c[s]
                  pp = p + se_w_p[s] 
                  #if is_in_mask(rows, cols, plns, rr, cc, pp, mask_data):
                  if is_in_image(rows, cols, plns, rr, cc, pp):
                      histogram_increment(histo, &pop, image[rr, cc, pp])
                      #printf('mode %d adding (%d,%d,%d) [%d]\n', mode, rr, cc, pp, image[rr,cc,pp]);
    
              for s in range(num_se_e):
                  rr = r + se_e_r[s] + 1
                  cc = c + se_e_c[s]
                  pp = p + se_e_p[s] 
                  #if is_in_mask(rows, cols, plns, rr, cc, pp, mask_data):
                  if is_in_image(rows, cols, plns, rr, cc, pp):
                      histogram_decrement(histo, &pop, image[rr, cc, pp])
                      #printf('mode %d removing (%d,%d,%d) [%d]\n', mode, rr, cc, pp, image[rr,cc,pp]);
    
              kernel(&out[r, c, p, 0], odepth, histo, pop, image[r, c, p], max_bin,
                     mid_bin, p0, p1, s0, s1)

          elif mode == 2:
              for s in range(num_se_n):
                  rr = r + se_n_r[s]
                  cc = c + se_n_c[s]
                  pp = p + se_n_p[s] 
                  #if is_in_mask(rows, cols, plns, rr, cc, pp, mask_data):
                  if is_in_image(rows, cols, plns, rr, cc, pp):
                      histogram_increment(histo, &pop, image[rr, cc, pp])
                      #printf('mode %d adding (%d,%d,%d) [%d]\n', mode, rr, cc, pp, image[rr,cc,pp]);
    
              for s in range(num_se_s):
                  rr = r + se_s_r[s] 
                  cc = c + se_s_c[s] - 1
                  pp = p + se_s_p[s] 
                  #if is_in_mask(rows, cols, plns, rr, cc, pp, mask_data):
                  if is_in_image(rows, cols, plns, rr, cc, pp):
                      histogram_decrement(histo, &pop, image[rr, cc, pp])
                      #printf('mode %d removing (%d,%d,%d) [%d]\n', mode, rr, cc, pp, image[rr,cc,pp]);
    
              kernel(&out[r, c, p, 0], odepth, histo, pop, image[r, c, p], max_bin,
                     mid_bin, p0, p1, s0, s1)
     
          elif mode == 3:
              for s in range(num_se_s):
                  rr = r + se_s_r[s]
                  cc = c + se_s_c[s]
                  pp = p + se_s_p[s] 
                  #if is_in_mask(rows, cols, plns, rr, cc, pp, mask_data):
                  if is_in_image(rows, cols, plns, rr, cc, pp):
                      histogram_increment(histo, &pop, image[rr, cc, pp])
                      #printf('mode %d adding (%d,%d,%d) [%d]\n', mode, rr, cc, pp, image[rr,cc,pp]);
    
              for s in range(num_se_n):
                  rr = r + se_n_r[s] 
                  cc = c + se_n_c[s] + 1
                  pp = p + se_n_p[s] 
                  #if is_in_mask(rows, cols, plns, rr, cc, pp, mask_data):
                  if is_in_image(rows, cols, plns, rr, cc, pp):
                      histogram_decrement(histo, &pop, image[rr, cc, pp])
                      #printf('mode %d removing (%d,%d,%d) [%d]\n', mode, rr, cc, pp, image[rr,cc,pp]);
    
              kernel(&out[r, c, p, 0], odepth, histo, pop, image[r, c, p], max_bin,
                     mid_bin, p0, p1, s0, s1)
          
                
          elif mode == 4:
              for s in range(num_se_u):
                  rr = r + se_u_r[s]
                  cc = c + se_u_c[s]
                  pp = p + se_u_p[s] 
                  #if is_in_mask(rows, cols, plns, rr, cc, pp, mask_data):
                  if is_in_image(rows, cols, plns, rr, cc, pp):
                      histogram_increment(histo, &pop, image[rr, cc, pp])
                      #printf('mode %d adding (%d,%d,%d) [%d]\n', mode, rr, cc, pp, image[rr,cc,pp]);
    
              for s in range(num_se_d):
                  rr = r + se_d_r[s]
                  cc = c + se_d_c[s]
                  pp = p + se_d_p[s] - 1
                  #if is_in_mask(rows, cols, plns, rr, cc, pp, mask_data):
                  if is_in_image(rows, cols, plns, rr, cc, pp):
                      histogram_decrement(histo, &pop, image[rr, cc, pp])
                      #printf('mode %d removing (%d,%d,%d) [%d]\n', mode, rr, cc, pp, image[rr,cc,pp]);
    
              kernel(&out[r, c, p, 0], odepth, histo, pop, image[r, c, p], max_bin,
                     mid_bin, p0, p1, s0, s1)

        
        # release memory allocated by malloc
        free(se_e_r)
        free(se_e_c)
        free(se_e_p)
        
        free(se_w_r)
        free(se_w_c)
        free(se_w_p)
        
        free(se_n_r)
        free(se_n_c)
        free(se_n_p)
        
        free(se_s_r)
        free(se_s_c)
        free(se_s_p)
        
        free(se_u_r)
        free(se_u_c)
        free(se_u_p)        
        
        free(se_d_r)
        free(se_d_c)
        free(se_d_p)
        
        free(histo)
