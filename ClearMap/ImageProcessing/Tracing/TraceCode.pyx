#distutils: language = c++
#distutils: sources = trace.cpp
#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False

cimport numpy as np
import numpy as np

from numpy cimport uint8_t, uint16_t, float32_t, double_t

from libcpp cimport bool as bool_t, long as long_t, int as int_t

ctypedef fused dtype_t_source:
    #uint8_t
    #uint16_t
    #float32_t
    double_t

#ctypedef fused dtype_t_tubness:
#    #uint8_t
#    #uint16_t
#    #float32_t
#    double_t

ctypedef fused dtype_t_mask:
     uint8_t
     uint16_t
     float32_t
     double_t

ctypedef fused index_t:
  Py_ssize_t


cdef extern from "trace.h":
    cdef cppclass Tracer[D, I]:
        Tracer()
        
        int run(D* source_, I shape_x_, I shape_y_, I shape_z_,
                             I stride_x_,I stride_y_,I sstride_z_,
                D* tubeness_,
                I start_x_, I start_y_, I start_z_,
                I goal_x_,  I goal_y_,  I goal_z_)
          
        int getPathSize()
        
        void getPath(I* path_array)
        
        double getPathQuality()
        
        double_t cost_per_distance        
        double_t minimum_cost_per_distance
        
        double_t tubeness_multiplier
        
        double_t minimal_tubeness
        
        long max_step

        bool_t verbose

cdef extern from "trace.h":
    cdef cppclass TracerToMask[D, I, M]:
        TracerToMask()
        
        int run(D* source_, I shape_x_, I shape_y_, I shape_z_,
                            I stride_x_,I stride_y_,I sstride_z_,
                D* tubeness_,
                I start_x_, I start_y_, I start_z_,
                M* mask_)
          
        int getPathSize()
        
        void getPath(I* path_array)
        
        double getPathQuality()
        
        double_t cost_per_distance        
        double_t minimum_cost_per_distance
        
        double_t tubeness_multiplier
        double_t minimal_tubeness
        
        long max_step

        bool_t verbose



def trace(dtype_t_source[:,:,:] source, 
          dtype_t_source[:,:,:] tubeness, 
          index_t[:] start, index_t[:] goal,
          double_t cost_per_distance, double_t minimum_cost_per_distance,
          double_t tubeness_multiplier, double_t minimal_tubeness, 
          bool_t quality,
          long_t max_step,
          bool_t verbose):
    
    #source  = np.ascontiguousarray(source);
    #tubness = np.ascontiguousarray(tubeness);
    
    cdef Tracer[dtype_t_source, index_t] tracer = Tracer[dtype_t_source, index_t]();
    tracer.cost_per_distance = cost_per_distance;
    tracer.minimum_cost_per_distance = minimum_cost_per_distance;    
    tracer.tubeness_multiplier = tubeness_multiplier;
    tracer.minimal_tubeness = minimal_tubeness;
    tracer.max_step = max_step;
    tracer.verbose = verbose;
    
    strides = np.array(source.strides) / source.itemsize;
    res = tracer.run(&source[0,0,0], source.shape[0], source.shape[1], source.shape[2],
                                     strides[0], strides[1], strides[2],
                     &tubeness[0,0,0],
                     start[0], start[1], start[2],
                     goal[0],  goal[1],  goal[2]);
   

                  
    if res < 0:
      if quality:
        return  np.zeros((0,3), dtype = int), 0;
      else:
        return np.zeros((0,3), dtype = int);
  
    n = tracer.getPathSize();
    cdef np.ndarray[index_t, ndim=2] path = np.zeros((n,3), dtype = int);
    tracer.getPath(&path[0,0]); 
    if quality:
      return path, tracer.getPathQuality();
    else:
      return path;
    
    
def traceToMask(dtype_t_source[:,:,:] source, 
                dtype_t_source[:,:,:] tubeness, 
                index_t[:] start, dtype_t_mask[:,:,:] mask,
                double_t cost_per_distance, double_t minimum_cost_per_distance,
                double_t tubeness_multiplier, double_t minimal_tubeness,
                bool_t quality,
                long_t max_step,
                bool_t verbose):
    
    #source  = np.ascontiguousarray(source);
    #tubness = np.ascontiguousarray(tubeness);
    #mask = np.ascontiguousarray(mask);
    
    cdef TracerToMask[dtype_t_source, index_t, dtype_t_mask] tracer = TracerToMask[dtype_t_source, index_t, dtype_t_mask]();
    tracer.cost_per_distance = cost_per_distance;
    tracer.minimum_cost_per_distance = minimum_cost_per_distance;    
    tracer.tubeness_multiplier = tubeness_multiplier;
    tracer.minimal_tubeness = minimal_tubeness;
    tracer.max_step = max_step;
    tracer.verbose = verbose;
    
    strides = np.array(source.strides) / source.itemsize;
    #print('strides: ', strides)
    res = tracer.run(&source[0,0,0], source.shape[0], source.shape[1], source.shape[2],
                                     strides[0], strides[1], strides[2],
                     &tubeness[0,0,0],
                     start[0], start[1], start[2],
                     &mask[0,0,0]);
                 
    if res < 0:
      if quality:
        return  np.zeros((0,3), dtype = int), 0;
      else:
        return np.zeros((0,3), dtype = int);
  
    n = tracer.getPathSize();
    cdef np.ndarray[index_t, ndim=2] path = np.zeros((n,3), dtype = int);
    tracer.getPath(&path[0,0]);    
    if quality:
      return path, tracer.getPathQuality();
    else:
      return path;