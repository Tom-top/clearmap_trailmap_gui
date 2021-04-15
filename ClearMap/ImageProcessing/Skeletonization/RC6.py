# -*- coding: utf-8 -*-
"""
Raynal Couprie 6-Subiteration Istmus based thinning

Reference:
  B. Raynal and M. Couprie, Istmus-Based 6-Directional Parallel Thinning Algorithms
"""
import os
import numpy as np
from multiprocessing import Pool, cpu_count

from Topology import rotate, cubeBool, cubeBase, rotations6

def match(cube):
  """Match one of the masks in the algorithm 
  
  Arguments:
    cube (3x3x3 bool array): the local binary image
  
  Returns:
    bool: True if one of the masks matches
  
  Note:
    cf. Fig 4 Raynal & Coupire
  """
  
  #M1'
  M1 = (cube[1,1,0] & cube[1,1,1] & 
        (not cube[0,0,2]) & (not cube[1,0,2]) & (not cube[2,0,2]) &
        (not cube[0,1,2]) & (not cube[1,1,2]) & (not cube[2,1,2]) &
        (not cube[0,2,2]) & (not cube[1,2,2]) & (not cube[2,2,2]));
  if M1:
    return True;
  
  # gerate rotations around z/vertical axis
  cuberots = [rotate(cube, axis = 2, steps = rot) for rot in range(4)];
  #print('Cube rotations:');
  #[printCube(c) for c in cuberots]  
  
  # M2' and all rotations
  for curo in cuberots:
    M2 = (curo[1,1,0] & curo[1,1,1] & curo[1,2,1] &
          (not curo[0,0,2]) & (not curo[1,0,2]) & (not curo[2,0,2]) &
          (not curo[0,1,2]) & (not curo[1,1,2]) & (not curo[2,1,2]));
    if M2:
      return True;
      
  # M3' and all rotations
  for curo in cuberots:
    M3 = (curo[1,1,0] & curo[1,1,1] & curo[1,2,1] & curo[2,1,1] &
          (not curo[0,0,2]) & (not curo[1,0,2]) &
          (not curo[0,1,2]) & (not curo[1,1,2]));
    if M3:
      return True;
      
  # M4' and all rotations
  for curo in cuberots:
    M4 = (curo[1,1,0] & curo[1,1,1] & curo[2,2,1] & curo[2,2,2] &
          (not curo[0,0,2]) & (not curo[1,0,2]) & (not curo[2,0,2]) &
          (not curo[0,1,2]) & (not curo[1,1,2]) & (not curo[2,1,2]) &
          (not curo[0,2,2]) & (not curo[1,2,2]));
    if M4:
      return True;
  
  # M5' and all rotations
  for curo in cuberots:
    M5 = (curo[1,2,0] & curo[1,1,1] & 
          (not curo[0,0,0]) & (not curo[1,0,0]) & (not curo[2,0,0]) &
          (not curo[1,1,0]) &
          (not curo[0,0,1]) & (not curo[1,0,1]) & (not curo[2,0,1]) &
          (not curo[0,0,2]) & (not curo[1,0,2]) & (not curo[2,0,2]) &
          (not curo[0,1,2]) & (not curo[1,1,2]) & (not curo[2,1,2]) &
          (not curo[0,2,2]) & (not curo[1,2,2]) & (not curo[2,2,2]));
    if M5:
      return True;
      
  # M6' and all rotations
  for curo in cuberots:
    M6 = (curo[2,1,0] & curo[1,2,0] & curo[1,1,1] &
          (not curo[0,0,0]) & (not curo[1,0,0]) &
          (not curo[0,1,0]) & (not curo[1,1,0]) &
          (not curo[0,0,1]) & (not curo[1,0,1]) &
          (not curo[0,1,1]) &
          (not curo[0,0,2]) & (not curo[1,0,2]) & (not curo[2,0,2]) &
          (not curo[0,1,2]) & (not curo[1,1,2]) & (not curo[2,1,2]) &
          (not curo[0,2,2]) & (not curo[1,2,2]) & (not curo[2,2,2]));
    if M6:
      return True;
      
  # M7' and all rotations
  for curo in cuberots:
    M7 = (curo[2,2,0] & curo[1,1,1] &
          (not curo[0,0,0]) & (not curo[1,0,0]) & (not curo[2,0,0]) &
          (not curo[0,1,0]) & (not curo[1,1,0]) & (not curo[2,1,0]) &
          (not curo[0,2,0]) & (not curo[1,2,0]) &
          (not curo[0,0,1]) & (not curo[1,0,1]) & (not curo[2,0,1]) &
          (not curo[0,1,1]) & (not curo[2,1,1]) &
          (not curo[0,2,1]) & (not curo[1,2,1]) & (not curo[2,2,1]) &
          (not curo[0,0,2]) & (not curo[1,0,2]) & (not curo[2,0,2]) &
          (not curo[0,1,2]) & (not curo[1,1,2]) & (not curo[2,1,2]) &
          (not curo[0,2,2]) & (not curo[1,2,2]) & (not curo[2,2,2]));
    if M7:
      return True;
    
    return False;
 

def matchIndex(index, verbose = True):
  if verbose and index % 2**14 == 0:
    print('RC6 LUT: %d / %d' % (index, 2**26));
  cube = cubeBool(index = index, center = True);
  return match(cube);


def generateLookupTable(base = None, verbose = True):
  """Generates lookup table for templates"""
   
  pool = Pool(cpu_count());
  lut = pool.map(matchIndex, range(2**26),chunksize=2**26/8/cpu_count());
  
  return np.array(lut, dtype = bool);


filename = "RC6.npy";
"""Filename for the Lookup table"""


def initializeMatches():
  """Initialize the module"""
  
  fn = os.path.join(os.path.dirname(os.path.abspath(__file__)), filename);
  if os.path.exists(fn):
    return np.load(fn);
  else:
    lut = generateLookupTable();
    np.save(fn, lut);
    return lut;


base = cubeBase(center = False);
"""Base kernel to multiply with cube to obtain index of cube"""

delete = initializeMatches();
"""Lookup table mapping cube index to its deleteability"""

keep = np.logical_not(delete);
"""Lookup table mapping cube index to its non-deleteability"""

rotations = rotations6(base);
"""Rotations of the base cube for the sub-iterations"""



### 1d Isthmuses 
  
def labelNeighbours26(data, label, x0,y0,z0, index):
  """Recursively label all 26 neighbours of (x0,y0,z0) with index"""
  shape = label.shape;
  for xp in range(max(0,-1+x0),min(2+x0, shape[0])):
    for yp in range(max(0,-1+y0),min(2+y0, shape[1])):
      for zp in range(max(0,-1+z0),min(2+z0, shape[2])):
        if data[xp,yp,zp] and label[xp,yp,zp] == 0:
          label[xp,yp,zp] = index;
          label = labelNeighbours26(data, label, xp,yp,zp, index);
  return label;
        

def labelComponents26(cube):
  """Counts number of 26-connected components in a binary array"""
  x,y,z = np.where(cube);
  label = np.zeros(cube.shape, dtype = 'uint8');
  ncomp = 0;
  for xp,yp,zp in zip(x,y,z):
    if label[xp,yp,zp] == 0:
      ncomp += 1;
      label = labelNeighbours26(cube, label, xp,yp,zp, ncomp);
  return ncomp, label
      

def countComponents26(cube):
  """Counts number of 26-connected components in a binary array"""
  n,l = labelComponents26(cube);
  return n;


def isthmus1D(cube):
  """Detect if the the center of the cube is a 1d istmus
  
  Arguments:
    cube (3x3x3 bool array): 26 neighbourhood of the point
  
  Returns:
    bool: True if the center point is a 1d istmus
    
  Note:
    The center point is a 1d istmuss if removal of center point gives rise to two 26-connected components
  """
  
  return countComponents26(cube) >= 2;

  
filename1DIsthmus = 'RC61DI.npy';


def isthmus1DFromIndex(index, verbose = True):
  if verbose and index % 2**14 == 0:
    print('RC6 LUT: %d / %d' % (index, 2**26));
  cube = cubeBool(index = index, center = False);
  return isthmus1D(cube);


def generateIstmus1D():
  """Generates lookup table for isthmuses"""
   
  pool = Pool(cpu_count());
  i1d = pool.map(isthmus1DFromIndex, range(2**26), chunksize=2**26/8/cpu_count());
  
  return np.array(i1d, dtype = bool);


def initializeIsthmus1D():
  """Initialize the 1-d istmus lookup table"""
  
  fn = os.path.join(os.path.dirname(os.path.abspath(__file__)), filename1DIsthmus);
  if os.path.exists(fn):
    return np.load(fn);
  else:
    i1d = generateIstmus1D();
    np.save(fn, i1d);
    return i1d;


isthmus = initializeIsthmus1D();
"""Lookup table for isthmus detection"""


### Skeletonization

import time
from ConvolvePointList import convolve_3d
#from Topology import n6

def skeletonize3D(data, steps = None, verbose = True):
    """Skeletonize a binary 3d array using RC6 istmus based algorithm
    
    Arguments:
      data (array): binary image
      steps (int or None): number of maximal iteration steps
      verbose (bool): if True print progrss info
      
    Returns:
      array: the skeleton
      array: the point coordinates of the skeleton nx3
    """
    
    if verbose:    
      print('#############################################################'); 
      print('Skeletonization RC6 [convolution]');
    tstart = time.time();
    
    # detect points
    points = np.array(np.nonzero(data)).T;
    if verbose:
      print('Foreground points: %d' % points.shape[0]);    
    
    #if info is not None:
    #  #birth = np.zeros(data.shape, dtype = 'uint16');
    #  death = np.zeros(data.shape, dtype = 'uint16');
    
    K = np.zeros(len(points.shape[0]), dtype = bool);
    
    # iterate
    if steps is None:
      steps = -1;

    step = 1;
    removed = 0;
    while True:
      if verbose:
        print('#############################################################');
        print('Iteration %d' % step);
      titer = time.time();
    
      # 1d istmusses on remaining points
      notK = np.logical_not(K);
      constrained[notK] = isthmus1D[convolve_3d(data, base, points[notK])];
            
      
      if verbose:
        print('-------------------------------------------------------------');
        print('Constrained %d' % constrained.sum());
      
      #if info is not None:
      #  b = birth[borderpoints[:,0], borderpoints[:,1], borderpoints[:,2]];
      #  bids = b == 0;
      #  birth[borderpoints[bids,0], borderpoints[bids,1], borderpoints[bids,2]] = step;
        
      # sub iterations over 6 directions
      remiter = 0;
      for i in range(6):
        if verbose:
          print('-------------------------------------------------------------');
          print('Sub-Iteration %d' % i);
        tsubiter = time.time();
        
        remborder = delete[convolve3D(data, rotations[i], free_points)];
        constrained[not_constrained] = rempoints = borderpoints[remborder];
        if verbose:
          print('Matched points: %d' % len(rempoints));
        

        data[rempoints[:,0], rempoints[:,1], rempoints[:,2]] = 0;
        keep[borderids[remborder]] = False;
        rem = len(rempoints);
        remiter += rem;
        removed += rem;
        if verbose:
          print('Deleted points: %d' % (rem));
          print('Sub-Iteration %d time: %0.2f s' % (i, time.time() - tsubiter));

      #update foreground
      points = points[keep];
      if verbose:
        print('Foreground points: %d' % points.shape[0]);  
        
      #death times
      #if info is not None:
      #  #remo = np.logical_not(keep);
      #  death[points[:,0], points[:,1], points[:,2]] = step;
      
      if verbose:
        print('-------------------------------------------------------------');
        print('Iteration time: %0.2f s' % (time.time() - titer));
      
      step += 1;
      if steps >= 0 and step >= steps:
        break
      if remiter == 0:
        break
  
    if verbose:
      print('#############################################################');
      print('Skeletonization time %0.2f s' % (time.time()-tstart));
      print('Total removed:   %d' % (removed));
      print('Total remaining: %d' % (len(points)));
  
    #if info is not None:
    #  #return data, points, birth, death
    #  return data, points, death
    #else:
    return data, points;





if __name__ == "__main__":
  import numpy as np;
  import RC6;
  reload(RC6);
  lut = RC6.generateLookupTable();
  np.save(RC6.filename, lut);
  
  reload(RC6);
  i1d = RC6.generateIstmus1D();
  np.save(RC6.filename1DIsthmus, i1d);
  
  import numpy as np;
  import Topology as top;
  data = np.load('test_bin.npy');
  data = top.deleteBorder(data);
  
  #data = np.zeros((5,5,5), dtype = bool);
  #data[1:4,1:4,1:4] = 1
  
  reload(RC6);
  skel, pts = RC6.skeletonize3D(data.copy(), verbose = True);
  
  
  #plot the result
  import os;
  cmpath = '/home/ckirst/Science/Projects/WholeBrainClearing/Vasculature/Analysis/ClearVessel'
  skpath = '/home/ckirst/Science/Projects/WholeBrainClearing/Vasculature/Analysis/ClearVessel/ClearMap/ImageProcessing/Skeletonization'
  os.chdir(cmpath);
  import ClearMap.GUI.DataViewer as dv;
  dv.DataViewer([data, skel])
  os.chdir(skpath)
  
  from mayavi import mlab
  mlab.figure(bgcolor=(1,1,1));
  mlab.contour3d(np.array(data, dtype = int), color = (0.5, 0.5, 0.5), opacity = 0.2);
  mlab.contour3d(np.array(skel, dtype = int), color = (1,0,0), opacity = 0.2);