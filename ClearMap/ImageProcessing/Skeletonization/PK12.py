# -*- coding: utf-8 -*-
"""
3d Skeletonization Pk12

This module implements the 3d parallel 12-Subiteration thinning algorithm by Palagy & Kuba via
parallel convolution of the data with a base template and lookup table matching

Reference:
  Palagyi & Kuba, A Parallel 3D 12-Subiteration Thinning Algorithm, Graphical Models and Image Processing 61, 199-221 (1999)
"""


import os
import numpy as np

from multiprocessing import Pool, cpu_count

from ClearMap.ImageProcessing.Topology.Topology3D import cubeBool, cubeBase, rotations12
 

def match(cube):
  """Match one of the masks in the algorithm 
  
  Arguments:
    cube (3x3x3 bool array): the local binary image
  
  Returns:
    bool: True if one of the masks matches
  
  Note:
    algorithm as in Palagyi & Kuba
  """
  #T1
  T1 = (cube[1,1,0] & cube[1,1,1] & 
        (cube[0,0,0] or cube[1,0,0] or cube[2,0,0] or
         cube[0,1,0] or cube[2,1,0] or
         cube[0,2,0] or cube[1,2,0] or cube[2,2,0] or
         cube[0,0,1] or cube[1,0,1] or cube[2,0,1] or
         cube[0,1,1] or cube[2,1,1] or 
         cube[0,2,1] or cube[1,2,1] or cube[2,2,1]) &
        (not cube[0,0,2]) & (not cube[1,0,2]) & (not cube[2,0,2]) &
        (not cube[0,1,2]) & (not cube[1,1,2]) & (not cube[2,1,2]) &
        (not cube[0,2,2]) & (not cube[1,2,2]) & (not cube[2,2,2]));
  if T1:
    return True;
  
  #T2
  T2 = (cube[1,1,1] & cube[1,2,1] & 
        (cube[0,1,0] or cube[1,1,0] or cube[2,1,0] or
         cube[0,2,0] or cube[1,2,0] or cube[2,2,0] or
         cube[0,1,1] or cube[2,1,1] or
         cube[0,2,1] or cube[2,2,1] or
         cube[0,1,2] or cube[1,1,2] or cube[2,1,2] or
         cube[0,2,2] or cube[1,2,2] or cube[2,2,2]) &
        (not cube[0,0,0]) & (not cube[1,0,0]) & (not cube[2,0,0]) &
        (not cube[0,0,1]) & (not cube[1,0,1]) & (not cube[2,0,1]) &
        (not cube[0,0,2]) & (not cube[1,0,2]) & (not cube[2,0,2]));
  if T2: 
    return True;
    
  #T3
  T3 = (cube[1,1,1] & cube[1,2,0] & 
        (cube[0,1,0] or cube[2,1,0] or
         cube[0,2,0] or cube[2,2,0] or
         cube[0,1,1] or cube[2,1,1] or
         cube[0,2,1] or cube[2,2,1]) &
        (not cube[0,0,0]) & (not cube[1,0,0]) & (not cube[2,0,0]) &
        (not cube[0,0,1]) & (not cube[1,0,1]) & (not cube[2,0,1]) &
        (not cube[0,0,2]) & (not cube[1,0,2]) & (not cube[2,0,2]) &
        (not cube[0,1,2]) & (not cube[1,1,2]) & ( not cube[2,1,2]) &
        (not cube[0,2,2]) & (not cube[1,2,2]) & (not cube[2,2,2]));
  if T3:
    return True;
  
  #T4
  T4 = (cube[1,1,0] & cube[1,1,1] & cube[1,2,1] & 
        ((not cube[0,0,1]) or (not cube[0,1,2])) &
        ((not cube[2,0,1]) or (not cube[2,1,2])) &
        (not cube[1,0,1]) & 
        (not cube[0,0,2]) & (not cube[1,0,2]) & (not cube[2,0,2]) &
        (not cube[1,1,2]));
  if T4:
    return True;    
  
  #T5
  T5 = (cube[1,1,0] & cube[1,1,1] & cube[1,2,1] & cube[2,0,2] &
        ((not cube[0,0,1]) or (not cube[0,1,2])) &
        (((not cube[2,0,1]) & cube[2,1,2]) or (cube[2,0,1] & (not cube[2,1,2]))) &
        (not cube[1,0,1]) & 
        (not cube[0,0,2]) & (not cube[1,0,2]) &
        (not cube[1,1,2]));
  if T5:
    return True;
    
  #T6
  T6 = (cube[1,1,0] & cube[1,1,1] & cube[1,2,1] & cube[0,0,2] &
        ((not cube[2,0,1]) or (not cube[2,1,2])) &
        (((not cube[0,0,1]) & cube[0,1,2]) or (cube[0,0,1] & (not cube[0,1,2]))) &
        (not cube[1,0,1]) & 
        (not cube[1,0,2]) & (not cube[2,0,2]) &
        (not cube[1,1,2]));
  if T6:
    return True;
    
  #T7
  T7 = (cube[1,1,0] & cube[1,1,1] & cube[2,1,1] &  cube[1,2,1] &
        ((not cube[0,0,1]) or (not cube[0,1,2])) &
        (not cube[1,0,1]) & 
        (not cube[0,0,2]) & (not cube[1,0,2]) &
        (not cube[1,1,2]));
  if T7:
    return True;
  
  #T8
  T8 = (cube[1,1,0] & cube[0,1,1] & cube[1,1,1] & cube[1,2,1] &
        ((not cube[2,0,1]) or (not cube[2,1,2])) &
        (not cube[1,0,1]) & 
        (not cube[1,0,2]) & (not cube[2,0,2]) &
        (not cube[1,1,2]));
  if T8:
    return True; 
    
  #T9
  T9 = (cube[1,1,0] & cube[1,1,1] & cube[2,1,1] & cube[0,0,2] & cube[1,2,1] &
        (((not cube[0,0,1]) & cube[0,1,2]) or (cube[0,0,1] & (not cube[0,1,2]))) &
        (not cube[1,0,1]) & 
        (not cube[1,0,2]) &
        (not cube[1,1,2]));
  if T9:
    return True;   
    
  #T10
  T10= (cube[1,1,0] & cube[0,1,1] & cube[1,1,1] & cube[2,0,2] & cube[1,2,1] &
        (((not cube[2,0,1]) & cube[2,1,2]) or (cube[2,0,1] & (not cube[2,1,2]))) &
        (not cube[1,0,1]) & 
        (not cube[1,0,2]) &
        (not cube[1,1,2]));
  if T10:
    return True;  
    
  #T11
  T11= (cube[2,1,0] & cube[1,1,1] & cube[1,2,0] &
        (not cube[0,0,0]) & (not cube[1,0,0]) & 
        (not cube[0,0,1]) & (not cube[1,0,1]) &
        (not cube[0,0,2]) & (not cube[1,0,2]) & (not cube[2,0,2]) &
        (not cube[0,1,2]) & (not cube[1,1,2]) & (not cube[2,1,2]) &
        (not cube[0,2,2]) & (not cube[1,2,2]) & (not cube[2,2,2]));
  if T11: 
    return True;
    
  #T12
  T12= (cube[0,1,0] & cube[1,2,0] & cube[1,1,1] &
        (not cube[1,0,0]) & (not cube[2,0,0]) & 
        (not cube[1,0,1]) & (not cube[2,0,1]) &
        (not cube[0,0,2]) & (not cube[1,0,2]) & (not cube[2,0,2]) &
        (not cube[0,1,2]) & (not cube[1,1,2]) & (not cube[2,1,2]) &
        (not cube[0,2,2]) & (not cube[1,2,2]) & (not cube[2,2,2]));
  if T12: 
    return True;
    
  #T13
  T13= (cube[1,2,0] & cube[1,1,1] & cube[2,2,1] &
        (not cube[0,0,0]) & (not cube[1,0,0]) & (not cube[2,0,0]) & 
        (not cube[0,0,1]) & (not cube[1,0,1]) & (not cube[2,0,1]) & 
        (not cube[0,0,2]) & (not cube[1,0,2]) & (not cube[2,0,2]) &
        (not cube[0,1,2]) & (not cube[1,1,2]) &
        (not cube[0,2,2]) & (not cube[1,2,2]));
  if T13: 
    return True; 
    
  #T14
  T14= (cube[1,2,0] & cube[1,1,1] & cube[0,2,1] &
        (not cube[0,0,0]) & (not cube[1,0,0]) & (not cube[2,0,0]) & 
        (not cube[0,0,1]) & (not cube[1,0,1]) & (not cube[2,0,1]) & 
        (not cube[0,0,2]) & (not cube[1,0,2]) & (not cube[2,0,2]) &
        (not cube[1,1,2]) & (not cube[2,1,2]) &
        (not cube[1,2,2]) & (not cube[2,2,2]));
  if T14: 
    return True; 
    
  return False;
 

def matchIndex(index, verbose = True):
  if verbose and index % 2**14 == 0:
    print('PK12 LUT: %d / %d' % (index, 2**26));
  cube = cubeBool(index = index, center = True);
  return match(cube);


def matchNonRemovable(index, verbose = True):
  if verbose and index % 2**14 == 0:
    print('PK12 LUT non-removables: %d / %d' % (index, 2**26));
  cube = top.cubeBool(index, center = False);
  n = cube.sum();
  if n < 2:
    return True;
  if n > 3:
    return False;
  x,y,z = np.where(cube); 
  if n == 2:
    if np.any(np.abs([x[1]-x[0], y[1]-y[0], z[1]-z[0]]) == 2):
      return True;
    else:
      return False;
  else:
     if np.any(np.abs([x[1]-x[0], y[1]-y[0], z[1]-z[0]]) == 2) and np.any(np.abs([x[2]-x[0], y[2]-y[0], z[2]-z[0]]) == 2) and np.any(np.abs([x[1]-x[2], y[1]-y[2], z[1]-z[2]]) == 2):
       return True;
     else:
       return False;


def generateLookupTable(function = matchIndex, verbose = True):
  """Generates lookup table for templates"""
   
  pool = Pool(cpu_count());
  lut = pool.map(function, range(2**26),chunksize=2**26/8/cpu_count());
  
  return np.array(lut, dtype = bool);


filename = "PK12.npy";
"""Filename for the look up table mapping a cube configuration to the deleatability of the center pixel"""

def initializeLookupTable(function = matchIndex, filename = filename):
  """Initialize the lookup table"""
  
  fn = os.path.join(os.path.dirname(os.path.abspath(__file__)), filename);
  if os.path.exists(fn):
    return np.load(fn);
  else:
    lut = generateLookupTable(function = function);
    np.save(fn, lut);
    return lut;


base = cubeBase(center = False);
"""Base kernel to multiply with cube to obtain index of cube"""

delete = initializeLookupTable();
"""Lookup table mapping cube index to its deleteability"""

keep = np.logical_not(delete);
"""Lookup table mapping cube index to its non-deleteability"""


filenameNonRemovable = "PK12nr.npy";
"""Filename for the lookup table mapping a cube configuration to the non-removeability of the center pixel"""

nonremovable = initializeLookupTable(filename = filenameNonRemovable, function = matchNonRemovable);
"""Lookup table mapping cube index to its non-removeability"""

consider = np.logical_not(nonremovable);
"""Lookup table mapping cube index to whether it needs to be considered further"""

rotations = rotations12(base);
"""Rotations of the base cube for the sub-iterations"""

#rotations_flat = [r.flatten()[np.hstack([range(13), range(14,27)])] for r in rotations];
#"""Rotations of the base cube as flat array with center removed for convolution in the sub-iterations"""


### Skeletonization

import time
import ClearMap.ImageProcessing.Topology.Topology3D as top
import ClearMap.DataProcessing.LargeData as ld
import ClearMap.DataProcessing.ConvolvePointList as cpl


def skeletonize3D(data, points = None, steps = None, removals = False, radii = True, verbose = True):
    """Skeletonize a binary 3d array using PK12 algorithm
    
    Arguments:
      data (array): binary image
      steps (int or None): number of maximal iteration steps (if None maximal reduction)
      removals (bool): if True returns the steps pixels in the input data where removed 
      radii (bool): if True the estimate of the local radius is returned
      memmap (bool): if True memmory map the large intermediate data structures
      verbose (bool): if True print progress info
      
    Returns:
      array: the skeleton
      array: the point coordinates of the skeleton nx3
    """
    
    if verbose:    
      print('#############################################################'); 
      print('Skeletonization PK12 [convolution]');
    tstart = time.time();
    
    # detect points
    #points = np.array(np.nonzero(data)).T;
    if points is None:
      points = ld.where(data);
    
    if verbose:
      print('Foreground points: %d, time %0.2f s' % (points.shape[0], time.time() - tstart));

    if removals is True or radii is True:
      #birth = np.zeros(data.shape, dtype = 'uint16');
      death = np.zeros(data.shape, dtype = 'uint16');
      with_info = True;
    else:
      with_info = False;
    
    # iterate
    step = 1;
    removed = 0;
    while True:
      if verbose:
        print('#############################################################');
        print('Iteration %d' % step);
      titer = time.time();
    
      border = cpl.convolve_3d_points(data, top.n6, points) < 6;
      borderpoints = points[border];
      borderids    = np.nonzero(border)[0];
      keep         = np.ones(len(border), dtype = bool);
      if verbose:  
        print('Border points: %d, time %0.2f s' % (len(borderpoints), time.time() - titer));
      
      #if info is not None:
      #  b = birth[borderpoints[:,0], borderpoints[:,1], borderpoints[:,2]];
      #  bids = b == 0;
      #  birth[borderpoints[bids,0], borderpoints[bids,1], borderpoints[bids,2]] = step;
        
      # sub iterations
      remiter = 0;
      for i in range(12):
        if verbose:
          print('-------------------------------------------------------------');
          print('Sub-Iteration %d' % i);
        tsubiter = time.time();
        tstep = tsubiter;
        
        remborder = delete[cpl.convolve_3d_points(data, rotations[i], borderpoints)];
        rempoints = borderpoints[remborder];
        if verbose:
          print('Matched points: %d, %0.2f s' % (len(rempoints), time.time() - tstep));
          tstep = time.time();
        

        data[rempoints[:,0], rempoints[:,1], rempoints[:,2]] = 0;
        keep[borderids[remborder]] = False;
        rem = len(rempoints);
        remiter += rem;
        removed += rem;
        if verbose:
          print('Deleted points: %d' % (rem));
          print('Sub-Iteration %d time: %0.2f s' % (i, time.time() - tstep));
          
        #death times
        if with_info is True:
          #remo = np.logical_not(keep);
          death[rempoints[:,0], rempoints[:,1], rempoints[:,2]] = 12 * step + i;

      #update foreground
      points = points[keep];
      if verbose:
        print('Foreground points: %d' % points.shape[0]);  
      
      if verbose:
        print('-------------------------------------------------------------');
        print('Iteration time: %0.2f s' % (time.time() - titer));
      
      step += 1;
      if step >= steps:
        break
      if remiter == 0:
        break
  
    if verbose:
      print('#############################################################');
      print('Skeletonization time %0.2f s' % (time.time()-tstart));
      print('Total removed:   %d' % (removed));
      print('Total remaining: %d' % (len(points)));
    
    
    if radii is True:
      #calculate average diameter as average death of neighbourhood
      radii = cpl.convolve_3d(death, np.array(top.n18, dtype = 'uint16'), points);
      if removals is True:
        return data, points, death, radii
      else:
        return data, points, radii
    else:
      if removals is True:
        return data, points, death
      else:
        return data, points;



def skeletonize3DIndex(data, points = None, steps = None, removals = False, radii = True, returnPoints = True, verbose = True):
    """Skeletonize a binary 3d array using PK12 algorithm using element index coordinates
    
    Arguments:
      data (array): binary image
      steps (int or None): number of maximal iteration steps (if None maximal reduction)
      removals (bool): if True returns the steps pixels in the input data where removed 
      radii (bool): if True the estimate of the local radius is returned
      verbose (bool): if True print progress info
      
    Returns:
      array: the skeleton
      array: the point coordinates of the skeleton nx3
    """
    
    if verbose:    
      print('#############################################################'); 
      print('Skeletonization PK12 [convolution, index]');
    tstart = time.time();
    
    dataflat = data.reshape(-1, order = 'A');
    
    # detect points
    if points is None:
      points = ld.where(dataflat);  
    npoints = points.shape[0];
    
    if verbose:
      print('Foreground points: %d, time %0.5f s' % (points.shape[0], time.time() - tstart));

    if removals is True or radii is True:
      #birth = np.zeros(data.shape, dtype = 'uint16');
      order = 'C';
      if data.flags.f_contiguous:
        order = 'F';
      death = np.zeros(data.shape, dtype = 'uint16', order = order);
      deathflat = death.reshape(-1, order = 'A')
      with_info = True;
    else:
      with_info = False;
    
    # iterate
    if steps is None:
      steps = -1;

    step = 1;
    nnonrem = 0;
    while True:
      if verbose:
        print('#############################################################');
        print('Iteration %d' % step);
      titer = time.time();
    
      border = cpl.convolve_3d_indices_if_smaller_than(data, top.n6, points, 6);
      borderpoints = points[border];
      #borderids    = np.nonzero(border)[0];
      borderids    = ld.where(border);
      keep         = np.ones(len(border), dtype = bool);
      if verbose:  
        print('Border points: %d, time %0.5f s' % (len(borderpoints), time.time() - titer));
      
      #if info is not None:
      #  b = birth[borderpoints[:,0], borderpoints[:,1], borderpoints[:,2]];
      #  bids = b == 0;
      #  birth[borderpoints[bids,0], borderpoints[bids,1], borderpoints[bids,2]] = step;
        
      # sub iterations
      remiter = 0;
      for i in range(12):
        if verbose:
          print('-------------------------------------------------------------');
          print('Sub-Iteration %d' % i);
        tsubiter = time.time();
        tstep = tsubiter;
        
        remborder = delete[cpl.convolve_3d_indices(data, rotations[i], borderpoints)];
        rempoints = borderpoints[remborder];
        if verbose:
          print('Matched points  : %d, %0.5f s' % (len(rempoints), time.time() - tstep));
          tstep = time.time();
        

        dataflat[rempoints] = 0;
        keep[borderids[remborder]] = False;
        rem = len(rempoints);
        remiter += rem;

        #death times
        if with_info is True:
          #remo = np.logical_not(keep);
          deathflat[rempoints] = 12 * step + i;
          
        if verbose:
          print('Sub-Iteration %d time: %0.5f s' % (i, time.time() - tsubiter));

      if verbose:
        print('-------------------------------------------------------------');

      #update foregroud
      points = points[keep];
      
      if step % 3 == 0:   
        npts = len(points);
        points = points[consider[cpl.convolve_3d_indices(data, base, points)]]; 
        nnonrem += npts - len(points)
        if verbose:
          print('Non-removable points: %d' % (npts - len(points)));
      
      if verbose:
        print('Foreground points   : %d' % points.shape[0]);  
      
      if verbose:
        print('-------------------------------------------------------------');
        print('Iteration time: %0.5f s' % (time.time() - titer));
      
      step += 1;
      if steps >= 0 and step >= steps:
        break
      if remiter == 0:
        break
    
    if verbose:
      print('#############################################################');
      print('Skeletonization done %0.5f s' % (time.time()-tstart));
    
   
    if verbose:
      print('Total time:      %0.5f s' % (time.time()-tstart));
      print('Total removed:   %d' % (npoints - (len(points) + nnonrem)));
      print('Total remaining: %d' % (len(points) + nnonrem));
    
    
    if radii is True or returnPoints is True:
      points = ld.where(dataflat)
    
    if radii is True:
      #calculate average diameter as death average death of neighbourhood     
      radii = cpl.convolve_3d_indices(death, top.n18, points, out_dtype = 'uint16');
    else:
      radii = None;
    
    if verbose:
      print('Final time:      %0.5f s' % (time.time()-tstart));
    
    ret = [data];
    if returnPoints:
      ret.append(points);
    if removals is True:
      ret.append(death);
    if radii is not None:
      ret.append(radii);
      
    if len(ret) > 1:
      return tuple(ret);
    else:
      return ret[0];



if __name__ == "__main__":
  import numpy as np;
  import ClearMap.ImageProcessing.Skeletonization.PK12 as PK12;
  import ClearMap.ImageProcessing.Skeletonization.Topology3D as top
  reload(PK12);
  #lut = PK12.generateLookupTable();
  #np.save(PK12.filename, lut);
  #lut.sum()
  
  #lutnr = PK12.generateLookupTable(function=PK12.matchNonRemovable, verbose = True);
  #np.save(PK12.filenameNonRemovable, lutnr);
  #lutnr.sum()
  
  #c = top.cubeBool(0, center = False)
  
  data = np.load('./ClearMap/Test/Skeletonization/test_bin.npy');
  data = top.deleteBorder(data);
  #data = np.zeros((5,5,5), dtype = bool);
  #data[1:4,1:4,1:4] = 1
  
  reload(PK12);
  skel, pts, radii = PK12.skeletonize3DIndex(data.copy(), verbose = True, radii = True);
  
  skel2, pts, radii = PK12.skeletonize3DIndexSlow(data.copy(), verbose = True, radii = True);  
  
  #plot the result
  import ClearMap.GUI.DataViewer as dv;
  dv.DataViewer([data, skel])  
  
  from mayavi import mlab
  mlab.figure(bgcolor=(1,1,1));
  mlab.contour3d(np.array(data, dtype = int), color = (0.5, 0.5, 0.5), opacity = 0.2);
  mlab.contour3d(np.array(skel, dtype = int), color = (1,0,0), opacity = 0.9);
  

  #looks very good !