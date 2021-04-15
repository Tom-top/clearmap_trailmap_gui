# -*- coding: utf-8 -*-
"""
Distance Transform module
"""
import numpy as np
import time

import ClearMap.ImageProcessing.Skeletonization.Topology3D as t3d
import ClearMap.DataProcessing.LargeData as ld
import ClearMap.DataProcessing.ConvolvePointList as cpl

def distanceTransform(data, points = None, steps = None, verbose = True):
  """Distrance tranform of a binary 3d array 
  
  Arguments:
  ----------
    data : array
      binary image
    points : array or None
      foreground points in data, if None determined from data
    verbose : bool
      if True print progress info
    
  Returns:
    array: distances as 1d array corresponding to the point indices
  """
    
  if verbose:    
    print('#############################################################'); 
    print('Distance Transform [index]');
  tstart = time.time();  
  
  data_flat = data.reshape(-1, order = 'A');
  
  if points is None:
    points = ld.where(data_flat);
    
  if verbose:
    print('Foreground points: %d, time %0.2f s' % (points.shape[0], time.time() - tstart));

  step = 1;
  while True:
    if verbose:
      print('#############################################################');
      print('Iteration %d' % step);
    titer = time.time();
  
     border = cpl.convolve3DIndexCondition(data, top.n6, points, 6);
     borderpoints = points[border];
 

    
    if verbose:
      print('-------------------------------------------------------------');
      print('Iteration time: %0.2f s' % (time.time() - titer));
      
      
      step += 1;
      if step >= steps:
        break