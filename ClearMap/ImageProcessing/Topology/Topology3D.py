# -*- coding: utf-8 -*-
"""
Topology module

Defines basic discrete topology utils
"""
__author__    = 'Christoph Kirst <ckirst@rockefeller.edu>'
__license__   = 'MIT License <http://www.opensource.org/licenses/mit-license.php>'
__copyright__ = 'Copyright (c) 2017 by Christoph Kirst, The Rockefeller University, New York City'


import numpy as np

### Neighbourhoods
n6 =  np.array([[[0,0,0],[0,1,0],[0,0,0]],
                [[0,1,0],[1,0,1],[0,1,0]], 
                [[0,0,0],[0,1,0],[0,0,0]]], dtype = bool);
"""6-Neighborhood excluding center"""

n18 = np.array([[[0,1,0],[1,1,1],[0,1,0]],
                [[1,1,1],[1,0,1],[1,1,1]], 
                [[0,1,0],[1,1,1],[0,1,0]]], dtype = bool);
"""18-Neighborhood excluding center"""

n26 = np.array([[[1,1,1],[1,1,1],[1,1,1]],
                [[1,1,1],[1,0,1],[1,1,1]],
                [[1,1,1],[1,1,1],[1,1,1]]], dtype = bool);
"""26-Neighborhood excluding center"""


### Label and Base 2 representations / indexing

def cubeLabel(center = True):
  """Returns an array with labels on the cube"""
  cube = np.zeros((3,3,3), dtype = int);
  if center:
    i = 0;
  else:
    i = 1;
  for z in range(3):
    for y in range(3):
      for x in range(3):
        if center is not None and x == 1 and y == 1 and z == 1:
          cube[x,y,z] = center;
        else:
          cube[x,y,z] = i;
          i+=1;
  return cube;


def cubeBase(center = True):
  """Returns an array with base 2 numbers on the cube for convolution and lut matching"""
  cube = np.zeros((3,3,3), dtype = int);
  k = 0;
  for z in range(3):
    for y in range(3):
      for x in range(3):
        if center is not None and x == 1 and y ==1 and z == 1:
          cube[x,y,z] = center;
        else:
          cube[x,y,z] = 2**k;
          k+=1;
  return cube;


def cubeBool(index, center = True):
  """Returns a boolean cube for the corresponding index"""
  cube = np.zeros((3,3,3), dtype = bool);
  d = 0;
  for z in range(3):
    for y in range(3):
      for x in range(3):
        if center is not None and x == 1 and y == 1 and z == 1:
          cube[x,y,z] = center;
        else:
          cube[x,y,z] = (index >> d) & 0x01;
          d += 1;
  return cube;
    

def cubeIndex(cube, center = False):
  """Returns index for a boolean cube"""
  return (cubeBase(center = center) * np.array(cube)).sum()

 
def nCubes():
  """Number of different cubes"""
  return np.sum(cubeBase(center = False)) + 1;

### neighbourhood lists

def neighbourXYZToIndex(x,y,z):
  return x + 3 * y + 9 * z;


def neighbourXYZFromIndex(index):
  index9 = index % 9;
  x = index9 % 3;
  y = index9 / 3;
  z = index / 9;
  return x,y,z


def neighbour6Indices(dtype = 'uint8'):
  """Indices as 6 Neighbourhood"""
  return np.array([4, 10, 12, 14, 16, 22], dtype = dtype); 


def neighbourhoodList(img, dtype = 'int64', verbose = False):
  """Return a list of x,y,z and list indices of the 26 neighbours"""
  if verbose:
    print("Generating neighbourhood...");    
  
  x,y,z = np.where(img);
  npts = len(x);
  
  if verbose:
    print('Creating labels...');
  
  label = np.full(img.shape, -1, dtype = dtype);
  label[x,y,z] = np.arange(npts);
  
  # calculate indices (if many voxels this is only 27 loops!)
  nhood = np.full((x.shape[0],27), -1, dtype = dtype);
  rg = range(3);
  direct = 0;
  for xx in rg:
    for yy in rg:
      for zz in rg:
        if verbose:
          print('Filling direction %d / 27' % direct);
          direct+=1;
        w = xx + yy * 3 + zz * 9;
        idx = x+xx-1; idy = y+yy-1; idz = z+zz-1;
        nhood[:,w] = label[idx,idy,idz];
        
  return (x,y,z,nhood);

 
def neighbourhoodListDelete(nhl, ids, changed = True):
  """Delete points in a neighbourhood list"""
  odirs = neighbourhoodOpposingDirections();
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
  else:
    return nhl;
  

def neighbourhoodOpposingDirections():
  """Returns a list of neighbour indices that are opposite of each other"""
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


def extractNeighbourhood(img,x,y,z):
  """Return the neighbourhoods of the indicated voxels
  
  Arguments:
    img (array): the 3d image
    x,y,z (n array): coordinates of the voxels to extract neighbourhoods from
  
  Returns:
    array (nx27 array): neighbourhoods
    
  Note:
    Assumes borders of the image are zero so that 0<x,y,z<w,h,d !
  """
  nhood = np.zeros((x.shape[0],27), dtype = bool);
  
  # calculate indices (if many voxels this is only 27 loops!)
  for xx in range(3):
    for yy in range(3):
      for zz in range(3):
        #w = _xyz_to_neighbourhood[xx,yy,zz];
        w = xx + 3 * yy + 9 * zz;
        idx = x+xx-1; idy = y+yy-1; idz = z+zz-1;
        nhood[:,w]=img[idx, idy, idz];
  
  nhood.shape = (nhood.shape[0], 3, 3, 3);
  nhood[:, 1, 1, 1] = 0;
  return nhood;







### Rotations

def rotate(cube, axis = 2, steps = 0):
  """Rotate a cube around an axis in 90 degrees steps"""
  cube = cube.copy();  
  
  steps = steps % 4;
  if steps == 0:
    return cube;
  
  elif axis == 0:
    if steps == 1:
      return cube[:, ::-1, :].swapaxes(1, 2)
    elif steps == 2:  # rotate 180 degrees around x
      return cube[:, ::-1, ::-1]
    elif steps == 3:  # rotate 270 degrees around x
      return cube.swapaxes(1, 2)[:, ::-1, :]
      
  elif axis == 1:
    if steps == 1:
      return cube[:, :, ::-1].swapaxes(2, 0)
    elif steps == 2:  # rotate 180 degrees around x
      return cube[::-1, :, ::-1]
    elif steps == 3:  # rotate 270 degrees around x
      return cube.swapaxes(2, 0)[:, :, ::-1]
      
  if axis == 2: # z axis rotation
    if steps == 1:
      return cube[::-1, :, :].swapaxes(0, 1)
    elif steps == 2:  # rotate 180 degrees around z
      return cube[::-1, ::-1, :]
    elif steps == 3:  # rotate 270 degrees around z
      return cube.swapaxes(0, 1)[::-1, :, :]


def rotations6(cube):
  """Generate rotations in 6 main directions"""
  
  rotU = cube.copy();
  rotD = rotate(cube, axis = 0, steps = 2); 
  rotN = rotate(cube, axis = 0, steps = 1);
  rotS = rotate(cube, axis = 0, steps = 3);
  rotE = rotate(cube, axis = 1, steps = 3);
  rotW = rotate(cube, axis = 1, steps = 1);
  
  return [rotU, rotN, rotW, rotD, rotS, rotE];


def rotations12(cube):
  """Generate rotations in 12 diagonal directions"""
  
  rotUS = cube.copy();
  rotUW = rotate(cube, axis = 2, steps = 1);  
  rotUN = rotate(cube, axis = 2, steps = 2); 
  rotUE = rotate(cube, axis = 2, steps = 3);  

  rotDS = rotate(cube,  axis = 1, steps = 2);
  rotDW = rotate(rotDS, axis = 2, steps = 1); 
  rotDN = rotate(rotDS, axis = 2, steps = 2); 
  rotDE = rotate(rotDS, axis = 2, steps = 3);

  rotSW = rotate(cube, axis = 1, steps = 1);   
  rotSE = rotate(cube, axis = 1, steps = 3); 

  rotNW = rotate(rotUN, axis = 1, steps = 1);
  rotNE = rotate(rotUN, axis = 1, steps = 3);
  
  return [rotUS, rotNE, rotDW,  rotSE, rotUW, rotDN,  rotSW, rotUN, rotDE,  rotNW, rotUE, rotDS];


def orientations():
  """Generate 13 orientations"""
  
  cube0 = np.zeros((3,3,3), dtype = bool);
  
  o1 = cube0.copy(); o1[0,1,1] = True;
  o2 = cube0.copy(); o2[0,2,1] = True;
  o3 = cube0.copy(); o3[1,2,1] = True;
  o4 = cube0.copy(); o4[2,2,1] = True;
  
  o5 = cube0.copy(); o5[0,0,0] = True;
  o6 = cube0.copy(); o6[1,0,0] = True;
  o7 = cube0.copy(); o7[2,0,0] = True;
  o8 = cube0.copy(); o8[0,1,0] = True;
  o9 = cube0.copy(); o9[1,1,0] = True;
  o10= cube0.copy(); o10[2,1,0] = True;
  o11= cube0.copy(); o11[0,2,0] = True;
  o12= cube0.copy(); o12[1,2,0] = True;
  o13= cube0.copy(); o13[2,2,0] = True;
  
  return [o1, o2, o3,  o4, o5, o6, o7, o8, o9,  o10, o11, o12, o13];


###  Uitility

def deleteBorder(data, value = 0):
  data[[0,-1],:,:] = value;
  data[:,[0,-1],:] = value;
  data[:,:,[0,-1]] = value;
  return data;

 
 
### Printing 

def printCube(cube):
  """Print the cube for debugging"""
  #for z in range(3):
  for y in range(2,-1,-1):
    print('{} {} {}'.format(cube[:,y,0], cube[:,y,1], cube[:,y,2]));
    #print ""
  print("---")



### Testing

if __name__ == "__main__":
  from importlib import reload
  import numpy as np
  import Topology as top
  
  reload(top)
  label = top.cubeLabel();
  top.printCube(label)
  
  # Test rotations
  c = np.zeros((3,3,3), dtype= bool);
  c[1,0,0] = True;
  top.printCube(c)
  
  cs = [top.rotate(c, axis = 2, steps = r) for r in range(4)];
  [top.printCube(cc) for cc in cs]
  
  reload(top)
  l = top.cubeLabel();
  rts = top.rotations6(l);
  
  [top.printCube(r) for r in rts]
  
  reload(top);
  b = top.cubeBool(6);
  i = top.cubeIndex(b);
  print(i,6)
  
  
  us = np.zeros((3,3,3), dtype = int);
  us[1,1,2] = 1;
  us[1,0,1] = 1;
  us[1,2,0] = 2;
  #top.printCube(us)
            
  r12 = top.rotations12(us);
  [top.printCube(c) for c in r12]
  