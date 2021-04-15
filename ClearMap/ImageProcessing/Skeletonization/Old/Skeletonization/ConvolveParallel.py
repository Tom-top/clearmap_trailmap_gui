# -*- coding: utf-8 -*-
"""
Convolution on a subset of points ina big array
"""

import numpy as np;
import pyximport; 
pyximport.install(setup_args={"include_dirs":np.get_include()}, reload_support=True)
 
from ConvolveParallelCy import cy_convolve3D;

def convolve3D(data, kernel, points, processes = 4):
    """Convolves binary data with a specified kernel at specific points only
    
    Arguments:
        data (array): 3d binary array to convolve
        kernel (array): kernel to convolve
        points (array): list of points to convolve
        
    Returns:
        array: list of results of convolution at specified points
    """
    
    return cy_convolve3D(data, kernel, points, processes);


def test():
  import numpy as np
  from ConvolveParallel import convolve3D
  data = np.array(np.random.rand(1500,500,500) > 0.75, dtype = int);
  data[[0,-1],:,:] = 0;
  data[:,[0,-1],:] = 0;
  data[:,:,[0,-1]] = 0;
  pts = np.where(data);
  pts = np.array(pts, dtype = int).T;
  
  from Topology import n6
  n6 = np.array(n6, dtype = int)
  
  import time;
  t0 = time.time();
  result = convolve3D(data, n6, pts)
  t1 = time.time();
  print '%f secs' % ((t1-t0));
  
  import scipy.ndimage
  t0 = time.time();
  good = scipy.ndimage.filters.convolve(np.array(data, dtype = float), np.array(n6, dtype = float), mode = 'constant')
  t1 = time.time();
  print '%f secs' % ((t1-t0));
  
  x,y,z = pts.T
  np.any(good[x,y,z] -result)
