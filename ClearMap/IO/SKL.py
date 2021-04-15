# -*- coding: utf-8 -*-
"""
Interface to write binary skeleton files 

Example:
    >>> import os
    >>> import ClearMap.Settings as settings
    >>> import ClearMap.IO.SKL as skl
    >>> filename = os.path.join(settings.ClearMapPath, 'Test/Data/NPY/skeleton.skl');
    >>> skeleton = skl.read(filename);
    >>> print skeleton.shape
    (10,20,30)

Note:
  Effectively saves a 3d binary array as list of non-zero entries 
"""


import numpy as np

#import ClearMap.IO as io;
#import ClearMap.IO.Region as rg;

import ClearMap.ParallelProcessing.SharedMemory as shm

import ClearMap.DataProcessing.LargeData as ld

def pointsToArray(points, shape = None, dtype = bool, order = 'F', shared = False):
  if shape is None and hasattr(points, '__len__'):
    points, shape = points;
  if shared:
    data = shm.zeros(shape, dtype = dtype, order = order);
  else:
    data = np.zeros(shape, dtype = dtype, order = order);
  
  data_f = np.reshape(data, -1, order = order);
  data_f[points] = 1;
  
  return data;


def arrayToPoints(data):
  points = ld.where(np.reshape(data, -1, order = 'A'));
  return (points, data.shape);


def dataType(filename, dtype = bool, **args):
  return np.dtype(dtype);


def read(filename, dtype = bool, order = 'F', shared = False, array = True, **args):
  if not isinstance(filename, str):
    raise NotImplementedError('reading only from files implemented');
  
  points = ld.load(filename);
  shape = points[:3];
  points = points[3:];
  
  if array:
    return pointsToArray(points, shape, dtype = dtype, order = order, shared = shared);
  else:
    return (points, shape)
      

def write(filename, data, shape = None, array = False):
  if not isinstance(filename, str):
    raise NotImplementedError('writing only to files implemented');
  
  if array and data.ndim != 3:
    raise RuntimeError('Expected a binary array');
  if not array and (data.ndim !=1 or shape is None or len(shape) != 3):
    raise RuntimeError('Expected a list of nonzero flat indices and shape');

  if array:
    points,shape = arrayToPoints(data);
  else:
    points = data;

  ld.save(filename, np.hstack([shape, points]));
  
  return filename
    

def dataSize(filename, **args):
  if not isinstance(filename, str):
    raise NotImplementedError('reading data size only implemented for files');
  
  data = np.lib.format.open_memmap(filename, mode = 'r');
  shape = np.copy(data[:3]);
  data = [];
  
  return shape;


def test():    
    """Test SKL module"""
    import os 
    import numpy as np
    import ClearMap.IO.SKL as skl
    reload(skl)
    
    fn = os.path.split(skl.__file__);
    fn = os.path.join(fn[0], 'skeleton.skl');
    
    skeleton = np.random.rand(10,20,30) > 0.7;
    points, shape = skl.arrayToPoints(skeleton);
    
    skl.write(fn, points, shape);  
    skeleton_2 = skl.read(fn, order = 'C');

    if np.all(skeleton == skeleton_2):
      print('SKL module working!');
    
    os.remove(fn);

if __name__ == "__main__":
    test();
    
