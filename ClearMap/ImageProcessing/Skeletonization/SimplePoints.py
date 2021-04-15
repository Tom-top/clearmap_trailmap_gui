# -*- coding: utf-8 -*-
"""
Calculation of Simple Points

The calculation is done by summing the Euler characteristics in each octant.

Reference: ???

Alternative: via G_26, see Bertrand Coupire Powerful Parallel Symmetric 2d Thinning Schemes Based on Criical Kernels, 2014
"""

import os
import numpy as np

def generateEulerLookupTable():
  """Euler charactersitic lookup table"""
  LUT = np.zeros(256, dtype = int);
  
  LUT[1]  =  1;
  LUT[3]  = -1;
  LUT[5]  = -1;
  LUT[7]  =  1;
  LUT[9]  = -3;
  LUT[11] = -1;
  LUT[13] = -1;
  LUT[15] =  1;
  LUT[17] = -1;
  LUT[19] =  1;
  LUT[21] =  1;
  LUT[23] = -1;
  LUT[25] =  3;
  LUT[27] =  1;
  LUT[29] =  1;
  LUT[31] = -1;
  LUT[33] = -3;
  LUT[35] = -1;
  LUT[37] =  3;
  LUT[39] =  1;
  LUT[41] =  1;
  LUT[43] = -1;
  LUT[45] =  3;
  LUT[47] =  1;
  LUT[49] = -1;
  LUT[51] =  1;
  
  LUT[53] =  1;
  LUT[55] = -1;
  LUT[57] =  3;
  LUT[59] =  1;
  LUT[61] =  1;
  LUT[63] = -1;
  LUT[65] = -3;
  LUT[67] =  3;
  LUT[69] = -1;
  LUT[71] =  1;
  LUT[73] =  1;
  LUT[75] =  3;
  LUT[77] = -1;
  LUT[79] =  1;
  LUT[81] = -1;
  LUT[83] =  1;
  LUT[85] =  1;
  LUT[87] = -1;
  LUT[89] =  3;
  LUT[91] =  1;
  LUT[93] =  1;
  LUT[95] = -1;
  LUT[97] =  1;
  LUT[99] =  3;
  LUT[101] =  3;
  LUT[103] =  1;
  
  LUT[105] =  5;
  LUT[107] =  3;
  LUT[109] =  3;
  LUT[111] =  1;
  LUT[113] = -1;
  LUT[115] =  1;
  LUT[117] =  1;
  LUT[119] = -1;
  LUT[121] =  3;
  LUT[123] =  1;
  LUT[125] =  1;
  LUT[127] = -1;
  LUT[129] = -7;
  LUT[131] = -1;
  LUT[133] = -1;
  LUT[135] =  1;
  LUT[137] = -3;
  LUT[139] = -1;
  LUT[141] = -1;
  LUT[143] =  1;
  LUT[145] = -1;
  LUT[147] =  1;
  LUT[149] =  1;
  LUT[151] = -1;
  LUT[153] =  3;
  LUT[155] =  1;
  
  LUT[157] =  1;
  LUT[159] = -1;
  LUT[161] = -3;
  LUT[163] = -1;
  LUT[165] =  3;
  LUT[167] =  1;
  LUT[169] =  1;
  LUT[171] = -1;
  LUT[173] =  3;
  LUT[175] =  1;
  LUT[177] = -1;
  LUT[179] =  1;
  LUT[181] =  1;
  LUT[183] = -1;
  LUT[185] =  3;
  LUT[187] =  1;
  LUT[189] =  1;
  LUT[191] = -1;
  LUT[193] = -3;
  LUT[195] =  3;
  LUT[197] = -1;
  LUT[199] =  1;
  LUT[201] =  1;
  LUT[203] =  3;
  LUT[205] = -1;
  LUT[207] =  1;
  
  LUT[209] = -1;
  LUT[211] =  1;
  LUT[213] =  1;
  LUT[215] = -1;
  LUT[217] =  3;
  LUT[219] =  1;
  LUT[221] =  1;
  LUT[223] = -1;
  LUT[225] =  1;
  LUT[227] =  3;
  LUT[229] =  3;
  LUT[231] =  1;
  LUT[233] =  5;
  LUT[235] =  3;
  LUT[237] =  3;
  LUT[239] =  1;
  LUT[241] = -1;
  LUT[243] =  1;
  LUT[245] =  1;
  LUT[247] = -1;
  LUT[249] =  3;
  LUT[251] =  1;
  LUT[253] =  1;
  LUT[255] = -1;
  
  return LUT;

euler = generateEulerLookupTable();


def euler_invariant(nbhoods):
  """Calculate Euler characteristic/invariance for a voxel in a 26-neighbourhood
    
  Arguments:
    nbhoods (nx27 array): 27-neighbourhoods to check
  
  Returns:
    array: True at entry n if n-th pixel is Euler invariant
  
  Note:
    The calculation is done by summing the Euler characteristics in each octant.
  """
  k = nbhoods.shape[0]; 
  #nbhoods = np.array(nbhoods, dtype = bool);
  
  eulerChar = np.zeros(k);
  # Octant SWU
  n = np.ones(k, dtype = int);
  n[nbhoods[:,24]] += 128;
  n[nbhoods[:,25]] += 64;
  n[nbhoods[:,15]] += 32;
  n[nbhoods[:,16]] += 16;
  n[nbhoods[:,21]] += 8;
  n[nbhoods[:,22]] += 4;
  n[nbhoods[:,12]] += 2;
  eulerChar += euler[n];
  
  # Octant SEU
  n = np.ones(k, dtype = int);
  n[nbhoods[:,26]] += 128;
  n[nbhoods[:,23]] += 64;
  n[nbhoods[:,17]] += 32;
  n[nbhoods[:,14]] += 16;
  n[nbhoods[:,25]] += 8;
  n[nbhoods[:,22]] += 4;
  n[nbhoods[:,16]] += 2;
  eulerChar += euler[n];
  
  # Octant NWU
  n = np.ones(k,dtype = int);
  n[nbhoods[:,18]] += 128;
  n[nbhoods[:,21]] += 64;
  n[nbhoods[:, 9]] += 32;
  n[nbhoods[:,12]] += 16;
  n[nbhoods[:,19]] += 8;
  n[nbhoods[:,22]] += 4;
  n[nbhoods[:,10]] += 2;
  eulerChar += euler[n];
  
  # Octant NEU
  n = np.ones(k, dtype = int);
  n[nbhoods[:,20]] += 128;
  n[nbhoods[:,23]] += 64;
  n[nbhoods[:,19]] += 32;
  n[nbhoods[:,22]] += 16;
  n[nbhoods[:,11]] += 8;
  n[nbhoods[:,14]] += 4;
  n[nbhoods[:,10]] += 2;
  eulerChar += euler[n];

  # Octant SWB
  n = np.ones(k, dtype = int);
  n[nbhoods[:, 6]] += 128;
  n[nbhoods[:,15]] += 64;
  n[nbhoods[:, 7]] += 32;
  n[nbhoods[:,16]] += 16;
  n[nbhoods[:, 3]] += 8;
  n[nbhoods[:,12]] += 4;
  n[nbhoods[:, 4]] += 2;
  eulerChar += euler[n];
  
  # Octant SEB
  n = np.ones(k, dtype = int);
  n[nbhoods[:, 8]] += 128;
  n[nbhoods[:, 7]] += 64;
  n[nbhoods[:,17]] += 32;
  n[nbhoods[:,16]] += 16;
  n[nbhoods[:, 5]] += 8;
  n[nbhoods[:, 4]] += 4;
  n[nbhoods[:,14]] += 2;
  eulerChar += euler[n];

  # Octant NWB
  n = np.ones(k, dtype = int);
  n[nbhoods[:, 0]] += 128;
  n[nbhoods[:, 9]] += 64;
  n[nbhoods[:, 3]] += 32;
  n[nbhoods[:,12]] += 16;
  n[nbhoods[:, 1]] += 8;
  n[nbhoods[:,10]] += 4;
  n[nbhoods[:, 4]] += 2;
  eulerChar += euler[n];

  # Octant NEB
  n = np.ones(k, dtype = int);
  n[nbhoods[:, 2]] += 128;
  n[nbhoods[:, 1]] += 64;
  n[nbhoods[:,11]] += 32;
  n[nbhoods[:,10]] += 16;
  n[nbhoods[:, 5]] += 8;
  n[nbhoods[:, 4]] += 4;
  n[nbhoods[:,14]] += 2;
  eulerChar += euler[n];

  #return eulerChar
  return eulerChar == 0;


#def pinfo(o, l, c):
#  print 'octant=%d, l=%r, c=%r' % (o,l,c.T);
  
#def pinfo(o,l,c):
#  pass;

def _octree_label_1(label, cube):
    #pinfo(1, label, cube)
    idx = cube[0,:] == 1;
    cube[0,idx] = label[idx];
        
    idx = cube[1,:] == 1;
    if np.any(idx):
        cube[1,idx] = label[idx];
        cube[:,idx] = _octree_label_2(label[idx],cube[:,idx]);
    
    idx = cube[3,:] == 1;
    if np.any(idx):
        cube[3,idx] = label[idx];
        cube[:,idx] = _octree_label_3(label[idx],cube[:,idx]);    
    
    idx = cube[4,:] == 1;
    if np.any(idx):
        cube[4,idx] = label[idx];
        cube[:,idx] = _octree_label_2(label[idx],cube[:,idx]);   
        cube[:,idx] = _octree_label_3(label[idx],cube[:,idx]);  
        cube[:,idx] = _octree_label_4(label[idx],cube[:,idx]);  
                
    idx = cube[9,:] == 1;
    if np.any(idx):
        cube[9,idx] = label[idx];
        cube[:,idx] = _octree_label_5(label[idx],cube[:,idx]);   
    
    idx = cube[10,:] == 1;
    if np.any(idx):
        cube[10,idx] = label[idx];
        cube[:,idx] = _octree_label_2(label[idx],cube[:,idx]);   
        cube[:,idx] = _octree_label_5(label[idx],cube[:,idx]);  
        cube[:,idx] = _octree_label_6(label[idx],cube[:,idx]);  

    idx = cube[12,:] == 1;
    if np.any(idx):
        cube[12,idx] = label[idx];
        cube[:,idx] = _octree_label_3(label[idx],cube[:,idx]);   
        cube[:,idx] = _octree_label_5(label[idx],cube[:,idx]);  
        cube[:,idx] = _octree_label_7(label[idx],cube[:,idx]); 
    
    return cube;


def _octree_label_2(label, cube):
    #pinfo(2, label, cube)
    idx = cube[1,:] == 1;
    if np.any(idx):
        cube[1,idx] = label[idx];
        cube[:,idx] = _octree_label_1(label[idx],cube[:,idx]);
    
    idx = cube[4,:] == 1;
    if np.any(idx):
        cube[4,idx] = label[idx];
        cube[:,idx] = _octree_label_1(label[idx],cube[:,idx]);   
        cube[:,idx] = _octree_label_3(label[idx],cube[:,idx]);  
        cube[:,idx] = _octree_label_4(label[idx],cube[:,idx]);  

    idx = cube[10,:] == 1;
    if np.any(idx):
        cube[10,idx] = label[idx];
        cube[:,idx] = _octree_label_1(label[idx],cube[:,idx]);   
        cube[:,idx] = _octree_label_5(label[idx],cube[:,idx]);  
        cube[:,idx] = _octree_label_6(label[idx],cube[:,idx]);  

    idx = cube[2,:] == 1;
    cube[2,idx] = label[idx];

    idx = cube[5,:] == 1;
    if np.any(idx):
        cube[5,idx] = label[idx];
        cube[:,idx] = _octree_label_4(label[idx],cube[:,idx]);   

    idx = cube[11,:] == 1;
    if np.any(idx):
        cube[11,idx] = label[idx];
        cube[:,idx] = _octree_label_6(label[idx],cube[:,idx]);   
    
    idx = cube[13,:] == 1;
    if np.any(idx):
        cube[13,idx] = label[idx];
        cube[:,idx] = _octree_label_4(label[idx],cube[:,idx]);   
        cube[:,idx] = _octree_label_6(label[idx],cube[:,idx]);  
        cube[:,idx] = _octree_label_8(label[idx],cube[:,idx]);  
    
    return cube;


def _octree_label_3(label, cube):
    #pinfo(3, label, cube)
    idx = cube[3,:] == 1;
    if np.any(idx):
        cube[3,idx] = label[idx];
        cube[:,idx] = _octree_label_1(label[idx],cube[:,idx]);   

    idx = cube[4,:] == 1;
    if np.any(idx):
        cube[4,idx] = label[idx];
        cube[:,idx] = _octree_label_1(label[idx],cube[:,idx]);   
        cube[:,idx] = _octree_label_2(label[idx],cube[:,idx]);  
        cube[:,idx] = _octree_label_4(label[idx],cube[:,idx]);  

    idx = cube[12,:] == 1;
    if np.any(idx):
        cube[12,idx] = label[idx];
        cube[:,idx] = _octree_label_1(label[idx],cube[:,idx]);   
        cube[:,idx] = _octree_label_5(label[idx],cube[:,idx]);  
        cube[:,idx] = _octree_label_7(label[idx],cube[:,idx]);  

    idx = cube[6,:] == 1;
    cube[6,idx] = label[idx];

    idx = cube[7,:] == 1;
    if np.any(idx):
        cube[7,idx] = label[idx];
        cube[:,idx] = _octree_label_4(label[idx],cube[:,idx]);   
    
    idx = cube[14,:] == 1;
    if np.any(idx):
        cube[14,idx] = label[idx];
        cube[:,idx] = _octree_label_7(label[idx],cube[:,idx]);  
    
    idx = cube[15,:] == 1;
    if np.any(idx):
        cube[15,idx] = label[idx];
        cube[:,idx] = _octree_label_4(label[idx],cube[:,idx]);   
        cube[:,idx] = _octree_label_7(label[idx],cube[:,idx]);  
        cube[:,idx] = _octree_label_8(label[idx],cube[:,idx]);
    
    return cube;


def _octree_label_4(label, cube):
    #pinfo(4, label, cube)
    idx = cube[4,:] == 1;
    if np.any(idx):
        cube[4,idx] = label[idx];
        cube[:,idx] = _octree_label_1(label[idx],cube[:,idx]);   
        cube[:,idx] = _octree_label_2(label[idx],cube[:,idx]);  
        cube[:,idx] = _octree_label_3(label[idx],cube[:,idx]);  
    
    idx = cube[5,:] == 1;
    if np.any(idx):
        cube[5,idx] = label[idx];
        cube[:,idx] = _octree_label_2(label[idx],cube[:,idx]);
    
    idx = cube[13,:] == 1;
    if np.any(idx):
        cube[13,idx] = label[idx];
        cube[:,idx] = _octree_label_2(label[idx],cube[:,idx]);   
        cube[:,idx] = _octree_label_6(label[idx],cube[:,idx]);  
        cube[:,idx] = _octree_label_8(label[idx],cube[:,idx]);  

    idx = cube[7,:] == 1;
    if np.any(idx):
        cube[7,idx] = label[idx];
        cube[:,idx] = _octree_label_3(label[idx],cube[:,idx]);   
    
    idx = cube[15,:] == 1;
    if np.any(idx):
        cube[15,idx] = label[idx];
        cube[:,idx] = _octree_label_3(label[idx],cube[:,idx]);   
        cube[:,idx] = _octree_label_7(label[idx],cube[:,idx]);  
        cube[:,idx] = _octree_label_8(label[idx],cube[:,idx]);  

    idx = cube[8,:] == 1;
    cube[8,idx] = label[idx];

    idx = cube[16,:] == 1;
    if np.any(idx):
        cube[16,idx] = label[idx];
        cube[:,idx] = _octree_label_8(label[idx],cube[:,idx]);   
    
    return cube;


def _octree_label_5(label, cube):
    #pinfo(5, label, cube)
    idx = cube[9,:] == 1;
    if np.any(idx):
        cube[9,idx] = label[idx];
        cube[:,idx] = _octree_label_1(label[idx],cube[:,idx]);   

    idx = cube[10,:] == 1;
    if np.any(idx):
        cube[10,idx] = label[idx];
        cube[:,idx] = _octree_label_1(label[idx],cube[:,idx]);   
        cube[:,idx] = _octree_label_2(label[idx],cube[:,idx]);  
        cube[:,idx] = _octree_label_6(label[idx],cube[:,idx]);  

    idx = cube[12,:] == 1;
    if np.any(idx):
        cube[12,idx] = label[idx];
        cube[:,idx] = _octree_label_1(label[idx],cube[:,idx]);   
        cube[:,idx] = _octree_label_3(label[idx],cube[:,idx]);  
        cube[:,idx] = _octree_label_7(label[idx],cube[:,idx]);  
    
    idx = cube[17,:] == 1;
    cube[17,idx] = label[idx];

    idx = cube[18,:] == 1;
    if np.any(idx):
        cube[18,idx] = label[idx];
        cube[:,idx] = _octree_label_6(label[idx],cube[:,idx]);   

    idx = cube[20,:] == 1;
    if np.any(idx):
        cube[20,idx] = label[idx];
        cube[:,idx] = _octree_label_7(label[idx],cube[:,idx]);  

    idx = cube[21,:] == 1;
    if np.any(idx):
        cube[21,idx] = label[idx];
        cube[:,idx] = _octree_label_6(label[idx],cube[:,idx]);   
        cube[:,idx] = _octree_label_7(label[idx],cube[:,idx]);  
        cube[:,idx] = _octree_label_8(label[idx],cube[:,idx]);  
    
    return cube;


def _octree_label_6(label, cube):
    #pinfo(6, label, cube)
    idx = cube[10,:] == 1;
    if np.any(idx):
        cube[10,idx] = label[idx];
        cube[:,idx] = _octree_label_1(label[idx],cube[:,idx]);   
        cube[:,idx] = _octree_label_2(label[idx],cube[:,idx]);  
        cube[:,idx] = _octree_label_5(label[idx],cube[:,idx]);  
    
    idx = cube[11,:] == 1;
    if np.any(idx):
        cube[11,idx] = label[idx];
        cube[:,idx] = _octree_label_2(label[idx],cube[:,idx]);   

    idx = cube[13,:] == 1;
    if np.any(idx):
        cube[13,idx] = label[idx];
        cube[:,idx] = _octree_label_2(label[idx],cube[:,idx]);   
        cube[:,idx] = _octree_label_4(label[idx],cube[:,idx]);  
        cube[:,idx] = _octree_label_8(label[idx],cube[:,idx]);  

    idx = cube[18,:] == 1;
    if np.any(idx):
        cube[18,idx] = label[idx];
        cube[:,idx] = _octree_label_5(label[idx],cube[:,idx]);   

    idx = cube[21,:] == 1;
    if np.any(idx):
        cube[21,idx] = label[idx];
        cube[:,idx] = _octree_label_5(label[idx],cube[:,idx]);   
        cube[:,idx] = _octree_label_7(label[idx],cube[:,idx]);  
        cube[:,idx] = _octree_label_8(label[idx],cube[:,idx]);  

    idx = cube[19,:] == 1;
    cube[19,idx] = label[idx];

    idx = cube[22,:] == 1;
    if np.any(idx):
        cube[22,idx] = label[idx];
        cube[:,idx] = _octree_label_8(label[idx],cube[:,idx]);  
    
    return cube;


def _octree_label_7(label, cube):
    #pinfo(7, label, cube)
    idx = cube[12,:] == 1;
    if np.any(idx):
        cube[12,idx] = label[idx];
        cube[:,idx] = _octree_label_1(label[idx],cube[:,idx]);   
        cube[:,idx] = _octree_label_3(label[idx],cube[:,idx]);  
        cube[:,idx] = _octree_label_5(label[idx],cube[:,idx]);  
    
    idx = cube[14,:] == 1;
    if np.any(idx):
        cube[14,idx] = label[idx];
        cube[:,idx] = _octree_label_3(label[idx],cube[:,idx]);   

    idx = cube[15,:] == 1;
    if np.any(idx):
        cube[15,idx] = label[idx];
        cube[:,idx] = _octree_label_3(label[idx],cube[:,idx]);   
        cube[:,idx] = _octree_label_4(label[idx],cube[:,idx]);  
        cube[:,idx] = _octree_label_8(label[idx],cube[:,idx]);  

    idx = cube[20,:] == 1;
    if np.any(idx):
        cube[20,idx] = label[idx];
        cube[:,idx] = _octree_label_5(label[idx],cube[:,idx]);   

    idx = cube[21,:] == 1;
    if np.any(idx):
        cube[21,idx] = label[idx];
        cube[:,idx] = _octree_label_5(label[idx],cube[:,idx]);   
        cube[:,idx] = _octree_label_6(label[idx],cube[:,idx]);  
        cube[:,idx] = _octree_label_8(label[idx],cube[:,idx]);  

    idx = cube[23,:] == 1;
    cube[23,idx] = label[idx];

    idx = cube[24,:] == 1;
    if np.any(idx):
        cube[24,idx] = label[idx];
        cube[:,idx] = _octree_label_8(label[idx],cube[:,idx]);
      
    return cube;


def _octree_label_8(label, cube):
    #pinfo(8, label, cube)
    idx = cube[13,:] == 1;
    if np.any(idx):
        cube[13,idx] = label[idx];
        cube[:,idx] = _octree_label_2(label[idx],cube[:,idx]);   
        cube[:,idx] = _octree_label_4(label[idx],cube[:,idx]);  
        cube[:,idx] = _octree_label_6(label[idx],cube[:,idx]);  
    
    idx = cube[15,:] == 1;
    if np.any(idx):
        cube[15,idx] = label[idx];
        cube[:,idx] = _octree_label_3(label[idx],cube[:,idx]);   
        cube[:,idx] = _octree_label_4(label[idx],cube[:,idx]);  
        cube[:,idx] = _octree_label_7(label[idx],cube[:,idx]);  

    idx = cube[16,:] == 1;
    if np.any(idx):
        cube[16,idx] = label[idx];
        cube[:,idx] = _octree_label_4(label[idx],cube[:,idx]);

    idx = cube[21,:] == 1;
    if np.any(idx):
        cube[21,idx] = label[idx];
        cube[:,idx] = _octree_label_5(label[idx],cube[:,idx]);   
        cube[:,idx] = _octree_label_6(label[idx],cube[:,idx]);  
        cube[:,idx] = _octree_label_7(label[idx],cube[:,idx]);  

    idx = cube[22,:] == 1;
    if np.any(idx):
        cube[22,idx] = label[idx];
        cube[:,idx] = _octree_label_6(label[idx],cube[:,idx]);   

    idx = cube[24,:] == 1;
    if np.any(idx):
        cube[24,idx] = label[idx];
        cube[:,idx] = _octree_label_7(label[idx],cube[:,idx]);   

    idx = cube[25,:] == 1;
    cube[25,idx] = label[idx];

    return cube;


# octree index assignment
#                         1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26
#                         0  1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 
_octree_index = np.array([0, 0, 1, 0, 0, 1, 2, 2, 3, 0, 0, 1, 0, 1, 2, 2, 3, 4, 4, 5, 4, 4, 5, 6, 6, 7], dtype = int) + 1;

_octree_index_to_fun = {1 : _octree_label_1, 
                        2 : _octree_label_2, 
                        3 : _octree_label_3,
                        4 : _octree_label_4, 
                        5 : _octree_label_5,
                        6 : _octree_label_6, 
                        7 : _octree_label_7,
                        8 : _octree_label_8
                        };

_octree_fun = [_octree_index_to_fun[i] for i in _octree_index];


def simple_point(nbhoods, verbose = False):
  """Checks if point is a simple point
  
  Arguments:
    nbhoods (nx27 array): neighbourhoods of the points to check
  
  Returns:
    n array: entry n is True if point in n-th neighbourhood is simple
  """
  
  # copy neighbors for labeling
  n_p = nbhoods.shape[0];
  p_is_simple = np.ones(n_p, dtype = bool);

  #neigbourhood without point
  cube = np.zeros((26, n_p));
  cube[0:13,:] = nbhoods[:,0:13].T;
  cube[13:26,:]= nbhoods[:,14:27].T;
  
  label = 2 * np.ones(n_p, dtype = int);

  for i in range(26): #loop over neighbours
    if verbose:
      print 'simple point iteration %d' % i;
    idx = np.logical_and(cube[i,:] == 1, p_is_simple);
    #print 'i=%d, idx=%r' % (i, idx)
    if np.any(idx):
      # start recursion with any octant that contains the point i
      cube[:,idx] = _octree_fun[i](label[idx], cube[:,idx]);
      label[idx] += 1;
      p_is_simple[label-2 >= 2] = False;
      # label-2; in [Lee94] is the number of connected compontents

  return p_is_simple;



def matchIndex(index, verbose = True):
  if verbose and index % 2**14 == 0:
    print('Simple Point LUT: %d / %d' % (index, 2**26));
  cube = cubeBool(index = index, center = True);
  
  return issimple(cube);


def generateLookupTable(base = None, verbose = True):
  """Generates lookup table for templates"""
   
  pool = Pool(cpu_count());
  lut = pool.map(matchIndex, range(2**26),chunksize=2**26/8/cpu_count());
  
  return np.array(lut, dtype = bool);


filename = "SimplePoints.npy";


def initializeLookupTable():
  """Initialize the lookup table"""
  
  fn = os.path.join(os.path.dirname(os.path.abspath(__file__)), filename);
  if os.path.exists(fn):
    return np.load(fn);
  else:
    lut = generateLookupTable();
    np.save(fn, lut);
    return lut;





_simple_point_filename = os.path.join(os.path.split(os.path.abspath(__file__))[0], 'skeletonize_3d_simple_points.npy');
_simple_point_cache = None;

def precompute_simple_points(filename = _simple_point_filename, verbose = False):
  """Precomputes weather a point is simple in a 3^3 cube"""
  global _simple_point_cache;
  
  if _simple_point_cache is not None:
    return _simple_point_cache;
  if filename is not None and os.path.isfile(filename):
    return load_simple_points(filename);
    
  if verbose:
    print "precomputing simple points and saving to %r" % filename;
  
  nh = np.zeros((2**26, 27), dtype = bool)
  ids = np.arange(2**26);
  #ids = ids[ids & (2**13) > 0];
  k = 1;
  for x in range(3):
    for y in range(3):
      for z in range(3):
        if x == 1 and y == 1 and z == 1:
          nh[:, x+3*y+9*z] = True;
        else:
          nh[:, x+3*y+9*z] = ids & k
          k *=2;
  
  _simple_point_cache = simple_point(nh, verbose = verbose);
  if filename is not None:
    np.save(filename, _simple_point_cache);

  return _simple_point_cache;


def get_neighbourhood(img,x,y,z):
  """Return the neighbourhoods of the indicated voxels
  
  Arguments:
    img (array): the 3d image
    x,y,z (n array): coordinates of the voxels to extract neighbourhoods from
  
  Returns:
    array (nx27 array): neighbourhoods
    
  Note:
    Assumes x,y,z coordinates lie not on border, i.e. 0<x,y,z<w,h,d !
  """
  #nhood = np.zeros((x.shape[0],27), dtype = bool);
  nhood = np.zeros((x.shape[0],27), dtype = img.dtype);
  
  # calculate indices (if many voxels this is only 27 loops!)
  rg = range(3);
  for xx in rg:
    for yy in rg:
      for zz in rg:
        #w = _xyz_to_neighbourhood[xx,yy,zz];
        w = xx + yy * 3 + zz * 9;
        idx = x+xx-1; idy = y+yy-1; idz = z+zz-1;
        nhood[:,w]=img[idx, idy, idz];
  
  return nhood;






def simple_point_fast(nh):
  """Calculates simple points using the precomupted simple points
  
  Arguments:
    nbhoods (nx27 array): neighbourhoods of the points to check
  
  Returns:
    n array: entry n is True if point in n-th neighbourhood is simple
    
  Note:
    Assumes the simple points are precomputed and loaded into _simple_points_cache. 
  """
  idx = np.zeros(nh.shape[0], dtype = int);
  for i in range(13):
    idx[nh[:,i]] += 2**i;
  for i in range(14,27):
    idx[nh[:,i]] += 2**(i-1);
  return _simple_point_cache[idx];


#_xyz_to_neighbourhood = np.zeros((3,3,3), dtype = int);
#for xx in range(3):
#  for yy in range(3):
#    for zz in range(3):
#      _xyz_to_neighbourhood[xx,yy,zz] = xx + yy * 3 + zz * 9;


#import ClearMap.ImageProcessing.Topology.Topology3D as top3d;
