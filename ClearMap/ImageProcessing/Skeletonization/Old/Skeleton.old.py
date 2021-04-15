# -*- coding: utf-8 -*-
"""
3d Skeletonization

Reference:
  Palagyi & Kuba, A Parallel 3D 12-Subiteration Thinning Algorithm, Graphical Models and Image Processing 61, 199-221 1999,
"""

import time
import numpy as np
import PK12

from ConvolvePointList import convolve3D
from Topology import n6


def deleteBorder(data, value = 0):
  data[[0,-1],:,:] = 0;
  data[:,[0,-1],:] = 0;
  data[:,:,[0,-1]] = 0;
  return data;

def neighbourXYZToIndex(x,y,z):
  return x + 3 * y + 9 * z;

def neighbourXYZFromIndex(index):
  index9 = index % 9;
  x = index9 % 3;
  y = index9 / 3;
  z = index / 9;
  return x,y,z

def neighbour6Indices(dtype = 'uint8'):
  return np.array([4, 10, 12, 14, 16, 22], dtype = dtype);

def neighbourList(img, dtype = 'int64', verbose = False): # memory intensive version
  """Return a list of x,y,z and list indices of the 26 neighbours"""
  if verbose:
    print("Generating neighbourhood...");    
  
  x,y,z = np.where(img);
  npts = len(x);
  
  if verbose:
    print('Creating labels...');
  
  label = np.full(img.shape, -1, dtype = dtype); # for 5000*5000*1500 need 2**64 indices
  label[x,y,z] = np.arange(npts);
  
  # calculate indices (if many voxels this is only 27 loops!)
  nhood = np.full((x.shape[0],27), -1, dtype = dtype);
  rg = range(3);
  direct = 0;
  for xx in rg:
    for yy in rg:
      for zz in rg:
        if verbose:
          print 'Filling direction %d / 27' % direct;
          direct+=1;
        w = xx + yy * 3 + zz * 9;
        idx = x+xx-1; idy = y+yy-1; idz = z+zz-1;
        nhood[:,w] = label[idx,idy,idz];
        
  return (x,y,z,nhood);
  

def neighboursOpposingDirectionList():
  dirs = np.zeros((27, 2), dtype = int);
  opdir = np.array([2,1,0]);
  rg = range(3);
  i = 0;
  for xx in rg:
    for yy in rg:
      for zz in rg:
        #w = _xyz_to_neighbourhood[xx,yy,zz];
        w = xx + yy * 3 + zz * 9;
        w2 = opdir[xx] + opdir[yy] * 3 + opdir[zz] * 9;
        dirs[i] = [w,w2];
        i+=1;
  return dirs;


def neighbourhoodListDelete(nhl, ids, changed = True):
  """Delete points in a neighbourhood list"""
  odirs = neighboursOpposingDirectionList();
  for i,oi in odirs:
    nhs = nhl[ids,i];
    nhi = nhs[nhs >= 0];
    nhl[nhi,oi] = -1;
 
  if changed: #collect non-deleted nodes with at least one neighbour deleted
    change = np.unique(nhl[ids].flatten());
    if len(change) > 0 and change[0] == -1:
      change = change[1:];
  
  # delete remaining neighbours of deleted pixels  
  nhl[ids] = -1;
  
  if changed:
    return nhl, change;
  return nhl;



def skeletonize3DList(data, steps = None, verbose = True, details = False):
    """Skeletonize a binary 3d array using a neighbour list based approach
    
    Arguments:
      data (array): binary image
      steps (int or None): number of maximal iteration steps
      verbose (bool): if True print progrss info
      
    Returns:
      array,array,array,array: positions and neighbourhood of skeleton pixels (x,y,z,n)
    """
    if verbose:    
      print('#############################################################'); 
      print('Skeletonization PK12 [neighbour list]');
    tstart = time.time();
    
    # create list of nozero points and their neighbourhood
    x,y,z,nhood = neighbourList(data, dtype = 'int64', verbose = verbose);
    
    # border points
    if verbose: 
      print('Initialize border pixel...');
    n6idx = neighbour6Indices(); 
    center = 1 + 3 + 9;
    border = np.where((nhood[:, n6idx] >= 0).sum(axis = 1) < 6)[0];
    
    # iterate
    if steps is None:
      steps = -1;

    step = 1;
    removed = 0;
    if details:
      remtime = np.full(len(x), 0, dtype = 'uint16');
      proctime = np.full(len(x), 0, dtype = 'uint16');
    
    while True:
      if verbose:
        print('#############################################################');
        print('Iteration %d' % step);
        print('Border points: %d' % len(border));
      titer = time.time();
      
      # sub iterations
      remiter = 0;
      consider = [];
      for i in range(12):
        if verbose:
          print('-------------------------------------------------------------');
          print('Sub-Iteration %d' % i);
        tsubiter = time.time();
        
        index = np.sum(PK12.rotations_flat[i] * (nhood[border] >= 0), axis = 1);
        match = PK12.delete[index];
        remids = border[match];
        rem = len(remids);
        remiter += rem;
        removed += rem;
        
        if verbose:
          print('Matched points: %d' % (rem));
          
        # delete points in neighbourhood
        nhood, changed = neighbourhoodListDelete(nhood, remids, changed = True);
        consider.append(changed);
        if details:
          remtime[remids] = step;
        
        if verbose:
          print('Deleted points: %d' % (rem));        
          print('Changed neighbours: %d' % len(changed));
        
        # detect new border ids to check
        #changedborder = np.where((nhood[changed, n6idx] >= 0).sum(axis = 1) < 6)[0];
        #newborder.append(changed[changedborder]);
        
        if verbose:
          print('Sub-Iteration %d time: %0.2f s' % (i, time.time() - tsubiter));
       
       
      if verbose:
        print('-------------------------------------------------------------');
        print('Iteration time: %0.2f s' % (time.time() - titer));
      
      
      step += 1;
      if steps >= 0 and step >= steps:
        break
      if remiter == 0:
        break
      
      # only consider cells whose neighbourhood has changed
      consider = np.unique(np.concatenate(consider));
      consider = consider[nhood[consider,center] >= 0];
      if details:
        proctime[consider] = step;
      if len(consider) == 0:
        break;
      border = consider[(nhood[consider][:, n6idx] >= 0).sum(axis = 1) < 6];
      if len(border) == 0:
        break;
        
  
    # cleanup
    valid = nhood[:,center] >=0;
    x = x[valid]; y = y[valid]; z = z[valid]; 
    if details:    
      proctime = proctime[valid];
    nhood = nhood[valid];
  
    if verbose:
      print('#############################################################');
      print('Skeletonization time %0.2f s' % (time.time()-tstart));
      print('Total removed:   %d' % (removed));
      print('Total remaining: %d' % (len(x)));
    
    if details:
      return x,y,z,nhood, proctime;
    else:
      return x,y,z,nhood
    


def skeletonize3DConvolve(data, steps = None, verbose = True, info = None):
    """Skeletonize a binary 3d array using a convolution approach
    
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
      print('Skeletonization PK12 [convolution]');
    tstart = time.time();
    
    # detect points
    points = np.array(np.nonzero(data)).T;
    if verbose:
      print('Foreground points: %d' % points.shape[0]);    
    
    if info is not None:
      #birth = np.zeros(data.shape, dtype = 'uint16');
      death = np.zeros(data.shape, dtype = 'uint16');
    
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
    
      border = convolve3D(data, n6, points) < 6;
      borderpoints = points[border];
      borderids    = np.nonzero(border)[0];
      keep         = np.ones(len(border), dtype = bool);
      if verbose:  
        print('Border points: %d' % len(borderpoints));
      
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
        
        remborder = PK12.delete[convolve3D(data, PK12.rotations[i], borderpoints)];
        rempoints = borderpoints[remborder];
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
      if info is not None:
        #remo = np.logical_not(keep);
        death[points[:,0], points[:,1], points[:,2]] = step;
      
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
  
    if info is not None:
      #return data, points, birth, death
      return data, points, death
    else:
      return data, points;



def test():
  import numpy as np;
  import Skeletonize as s;
  data = np.load('test_bin.npy');
  data = s.deleteBorder(data);
  
  #data = np.zeros((5,5,5), dtype = bool);
  #data[1:4,1:4,1:4] = 1
  
  reload(s);
  x,y,z,nh = s.skeletonize3DList(data, verbose = True);
  skel = np.zeros(data.shape, dtype = bool);
  skel[x,y,z] = True;

  reload(s);
  skel2, pts, death = s.skeletonize3DConvolve(data.copy(), verbose = True, info = True);

  #plot the result
  import os;
  cmpath = '/home/ckirst/Science/Projects/WholeBrainClearing/Vasculature/Analysis/ClearVessel'
  skpath = '/home/ckirst/Science/Projects/WholeBrainClearing/Vasculature/Analysis/ClearVessel/ClearMap/ImageProcessing/Skeletonization'
  os.chdir(cmpath);
  import ClearMap.GUI.DataViewer as dv;
  dv.DataViewer([data, skel1])
  os.chdir(skpath)
  
  
  from mayavi import mlab
  mlab.figure(bgcolor=(1,1,1));
  mlab.contour3d(np.array(data, dtype = int), color = (0.5, 0.5, 0.5), opacity = 0.2);
  mlab.contour3d(np.array(skel2, dtype = int), color = (1,0,0), opacity = 0.2);
  
  import os;
  cmpath = '/home/ckirst/Science/Projects/WholeBrainClearing/Vasculature/Analysis/ClearVessel'
  skpath = '/home/ckirst/Science/Projects/WholeBrainClearing/Vasculature/Analysis/ClearVessel/ClearMap/ImageProcessing/Skeletonization'
  os.chdir(cmpath);
  import ClearMap.GUI.DataViewer as dv;
  dv.DataViewer([data, skel, death])
  os.chdir(skpath)
  
  