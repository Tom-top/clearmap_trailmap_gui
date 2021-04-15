# -*- coding: utf-8 -*-
"""
2D Topology module

Defines basic discrete topology utils
"""
import numpy as np

# Topology numbers T8 and T4bar in 2D
#Configuration:
# 3 2 1			
# 4 X 0
# 5 6 7
# 
# configuration

t4barLUT = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 1, 1, 1, 1, 1, 2, 2, 2, 1, 1, 1,
       1, 1, 2, 2, 2, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 2, 2,
       2, 2, 1, 2, 2, 2, 1, 1, 1, 1, 1, 2, 2, 2, 1, 1, 1, 1, 1, 1, 2, 1, 2,
       1, 2, 1, 2, 2, 3, 2, 2, 1, 2, 1, 1, 1, 2, 1, 1, 0, 1, 0, 1, 1, 2, 1,
       1, 0, 1, 0, 1, 1, 2, 1, 2, 1, 2, 1, 2, 2, 3, 2, 2, 1, 2, 1, 1, 1, 2,
       1, 1, 0, 1, 0, 1, 1, 2, 1, 1, 0, 1, 0, 1, 1, 2, 1, 2, 1, 2, 1, 2, 2,
       3, 2, 2, 1, 2, 1, 2, 2, 3, 2, 2, 1, 2, 1, 2, 2, 3, 2, 2, 1, 2, 1, 2,
       2, 3, 2, 3, 2, 3, 2, 3, 3, 4, 3, 3, 2, 3, 2, 2, 2, 3, 2, 2, 1, 2, 1,
       2, 2, 3, 2, 2, 1, 2, 1, 1, 1, 2, 1, 2, 1, 2, 1, 2, 2, 3, 2, 2, 1, 2,
       1, 1, 1, 2, 1, 1, 0, 1, 0, 1, 1, 2, 1, 1, 0, 1, 0, 1, 1, 2, 1, 2, 1,
       2, 1, 2, 2, 3, 2, 2, 1, 2, 1, 1, 1, 2, 1, 1, 0, 1, 0, 1, 1, 2, 1, 1,
       0, 1, 0]);
  
t8LUT = np.array([0, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 1, 1, 1, 1, 1, 2, 2, 2, 1, 1, 1,
       1, 1, 2, 2, 2, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 2, 2,
       2, 2, 1, 2, 2, 2, 1, 1, 1, 1, 1, 2, 2, 2, 1, 1, 1, 1, 1, 1, 2, 1, 2,
       1, 2, 1, 2, 2, 3, 2, 2, 1, 2, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 2, 1,
       1, 1, 1, 1, 1, 1, 2, 1, 2, 1, 2, 1, 2, 2, 3, 2, 2, 1, 2, 1, 1, 1, 2,
       1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 2, 1, 2, 1, 2, 1, 2, 2,
       3, 2, 2, 1, 2, 1, 2, 2, 3, 2, 2, 1, 2, 1, 2, 2, 3, 2, 2, 1, 2, 1, 2,
       2, 3, 2, 3, 2, 3, 2, 3, 3, 4, 3, 3, 2, 3, 2, 2, 2, 3, 2, 2, 1, 2, 1,
       2, 2, 3, 2, 2, 1, 2, 1, 1, 1, 2, 1, 2, 1, 2, 1, 2, 2, 3, 2, 2, 1, 2,
       1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 2, 1, 2, 1,
       2, 1, 2, 2, 3, 2, 2, 1, 2, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1,
       1, 1, 1])

#def planeBase(center = True):
#  """Returns an array with base 2 numbers on the plane for convolution and lut matching"""
#  plane = np.zeros((3,3), dtype = int);
#  k = 0;
#  for y in range(3):
#    for x in range(3):
#      if center is not None and x == 1 and y ==1:
#        plane[x,y] = center;
#      else:
#        plane[x,y] = 2**k;
#        k+=1;
#  return plane;

def planeBase(center = False):
  if center is not None:
    base = np.array([[3,2,1],[4,0,0],[5,6,7]]);
    base = np.power(2, base);
    base[1,1] = center;
  else:
    base = np.array([[4,3,2],[5,0,1],[6,7,8]]);
    base = np.power(2, base);
  return base;


#def planeBool(index, center = False):
#  """Returns a boolean plane for the corresponding index"""
#  plane = np.zeros((3,3), dtype = bool);
#  d = 0;
#  for y in range(3):
#    for x in range(3):
#      if center is not None and x == 1 and y == 1 and z == 1:
#        plane[x,y] = center;
#      else:
#        plane[x,y] = (index >> d) & 0x01;
#        d += 1;
#  return plane;
    
def planeIndex(plane, center = False):
  """Returns index for a boolean cube"""
  return (planeBase(center = center) * np.array(plane)).sum()



if __name__ == "__main__":
  import Topology2D as t2;
  
  p = [[1,0,0],[0,0,1],[1,0,0]];
  pi = t2.planeIndex(p)
  print('t4b: %d, t8: %d' % (t2.t4barLUT[pi], t2.t8LUT[pi]))
  
  p = [[1,1,0],[0,0,1],[1,0,0]];
  pi = t2.planeIndex(p)
  print('t4b: %d, t8: %d' % (t2.t4barLUT[pi], t2.t8LUT[pi]))
  


###Neighbourhoods 


def extractNeighbourhood(img,x,y):
  """Return the neighbourhoods of the indicated voxels
  
  Arguments:
    img (array): the 2d image
    x,y (n array): coordinates of the voxels to extract neighbourhoods from
  
  Returns:
    array (nx9 array): neighbourhoods
    
  Note:
    Assumes borders of the image are zero so that 0<x,y<w,h !
  """
  nhood = np.zeros((x.shape[0],9), dtype = bool);
  
  # calculate indices (if many voxels this is only 9 loops!)
  for xx in range(3):
    for yy in range(3):
        #w = _xyz_to_neighbourhood[xx,yy,zz];
        w = 3 * xx + yy;
        idx = x+xx-1; idy = y+yy-1;
        nhood[:,w]=img[idx, idy];
  
  nhood.shape = (nhood.shape[0], 3, 3);
  nhood[:, 1, 1] = 0;
  return nhood;