# -*- coding: utf-8 -*-
"""
Convolution on a subset of points 

Paralllel convolution of a kernel at specified points of the data only.
Useful to speed up processing in large arrays and only a smaller number of convolution points

Notes:
  This module is heavily used for the skeletonization and sparse pont measurements

"""
__author__    = 'Christoph Kirst <ckirst@rockefeller.edu>'
__license__   = 'MIT License <http://www.opensource.org/licenses/mit-license.php>'
__copyright__ = 'Copyright 2017 by Christoph Kirst, The Rockefeller University, New York City'


import numpy as np;
from multiprocessing import cpu_count;

import pyximport; 
pyximport.install(setup_args={"include_dirs": [np.get_include()]}, reload_support=True)
 
import ClearMap.DataProcessing.ConvolvePointListCode as code


#def kernel_array_to_list(kernel):
#  """Converts a 3d kernel to a list of offsets and non-zero kernel values"""
#  kernel = np.array(kernel);
#  if kernel.ndim != 3:
#    raise ValueError('Kernel dimension %d != 3' % kernel.ndim);
#  
#  kpos = np.array(np.where(kernel));
#  x,y,z = kpos;
#  kval = kernel[x,y,z];
#  
#  return kpos, kval
  

def convolve_3d(data, kernel, points = None, indices = None, x = None, y = None, z = None, out = None, out_dtype = None, strides = None, check_border = True, processes = cpu_count()):
  """Convolves data with a specified kernel at specific points only
    
  Arguments:
        data (array): 3d binary array to convolve
        kernel (array): binary orinteger kernel to convolve
        points (array): list of points to convolve
        indices (array): indices to convolve
        x,y,z (array): arrays of x,y,z coordinates of points to convolve on
        out (array) : optional arry to write data into
        out_dtype (dtype): optional data type of the output array, the kernel type is used as default
        strides (array): the straides of the data in case its given as a 1d list
        check_border (bool): if True checks each kernel element to be inside the data array shape
        processes (int): number of processes
        
  Returns:
        array: list of results of convolution at specified points
    
  Note:
        Either points x,y,z or an index array needs to be given. This function wraps more specialized functions
        
  See also:
        convolve_3d_points, convolve_3d_xyz, convolve_3d_indices
  """
    
  if points is not None and points.ndim == 1:
      indices = points;
      points = None;
    
  if points is not None:
    return convolve_3d_points(data, kernel, points, out = out, out_dtype = out_dtype, check_border = check_border, processes = processes);
  elif indices is not None:
    return convolve_3d_indices(data, kernel, indices, out = out, out_dtype = out_dtype, check_border = check_border, processes = processes);
  elif x is not None and y is not None and z is not None:
    return convolve_3d_xyz(data, kernel, x, y, z, out = out, out_dtype = out_dtype, check_border = check_border, processes = processes);
  else:
    raise RuntimeError('Positions expected to be given as points, index or x,y,z arrays');
  


def convolve_3d_points(data, kernel, points, out = None, out_dtype = None, check_border = True, processes = cpu_count()):
  """Convolves data with a specified kernel at specific points only
    
  Arguments:
        data (array): 3d binary array to convolve
        kernel (array): binary orinteger kernel to convolve
        points (array): list of points to convolve
        processes (int): number of processors
    
  Returns:
        array: list of results of convolution at specified points
    
  Note:
        cython does not support bools -> use view on uint8 as numpy does
  """
  
  if data.dtype == bool:
    d = data.view('uint8');
  else:
    d = data;
  
  npts = points.shape[0];  
    
  if out is None:
    if out_dtype is None:
      out_dtype = kernel.dtype;
    out = np.zeros(npts, dtype = out_dtype);
  
  if out.shape[0] != npts:
     raise RuntimeError('The output has not the expected size of %d but %d' % (npts, out.shape[0]));
  
  if out.dtype == bool:
    o = out.view('uint8');
  else:
    o = out;

  if kernel.dtype == bool:
    k = np.array(kernel, 'uint8');
  else:
    k = kernel;
  
  if processes is None:
    processes = cpu_count();
  
  if check_border:
    code.convolve_3d_points(d, k, points, o, processes);
  else:
    code.convolve_3d_points_no_check(d, k, points, o, processes);
  
  return out;


def convolve_3d_xyz(data, kernel, x, y, z, out = None, out_dtype = None, check_border = True, processes = cpu_count()):
  """Convolves data with a specified kernel at specific points only
    
  Arguments:
        data (array): 3d binary array to convolve
        kernel (array): binary orinteger kernel to convolve
        points (array): list of points to convolve
        processes (int): number of processors
    
  Returns:
        array: list of results of convolution at specified points
    
  Note:
        cython does not support bools -> use view on uint8 as numpy does
  """
  
  if data.dtype == bool:
    d = data.view('uint8');
  else:
    d = data;
    
  npts = len(x);
    
  if out is None:
    if out_dtype is None:
      out_dtype = kernel.dtype;
    out = np.zeros(npts, dtype = out_dtype);
  
  if out.shape[0] != npts or len(y) != npts or len(z) != npts:
     raise RuntimeError('The output has size %d and does not match the x,y,z coordinates of sizes: %d = %d = %d' % (out.shape[0], len(x), len(y), len(z)));
  
  if out.dtype == bool:
    o = out.view('uint8');
  else:
    o = out;

  if kernel.dtype == bool:
    k = np.array(kernel, 'uint8');
  else:
    k = kernel;
    
  if processes is None:
    processes = cpu_count();
  
  if check_border:
    code.convolve_3d_xyz(d, k, x, y, z, o, processes);
  else:
    code.convolve_3_xyz_no_check(d, k, x, y, z, o, processes);
  
  return out;


def convolve_3d_indices(data, kernel, indices, out = None, out_dtype = None, strides = None, check_border = True, processes = cpu_count()):
  """Convolves data with a specified kernel at specific points given by a flat array index
    
  Arguments:
        data (array): 3d binary array to convolve
        kernel (array): binary orinteger kernel to convolve
        indices (array): list of points to convolve given by the flat array coordinates
        x,y,z (array): arrays of x,y,z coordinates of points to convolve on
        processes (int): number of processors
    
  Returns:
        array: list of results of convolution at specified points
    
  Note:
        cython does not support bools -> use view on uint8 as numpy does
  """
  d = data.reshape(-1, order = 'A');
  if data.dtype == bool:
    d = d.view('uint8');
    
  npts = indices.shape[0];
    
  if out is None:
    if out_dtype is None:
      out_dtype = kernel.dtype;
    out = np.zeros(npts, dtype = out_dtype);
  
  if out.shape[0] != npts:
     raise RuntimeError('The output has not the expected size of %d but %d' % (npts, out.shape[0]));
  
  if out.dtype == bool:
    o = out.view('uint8');
  else:
    o = out;

  if kernel.dtype == bool:
    k = np.array(kernel, 'uint8');
  else:
    k = kernel; 

  if processes is None:
    processes = cpu_count();
  
  if strides is None:
    strides = np.array(data.strides) / np.array(data.itemsize);
  
  #print d.dtype, strides.dtype, kernel.dtype, o.dtype
  if check_border:
    code.convolve_3d_indices(d, strides, k, indices, o, processes);
  else:
    code.convolve_3d_indices_no_check(d, strides, k, indices, o, processes);
  
  return out;



def convolve_3d_indices_if_smaller_than(data, kernel, indices, max_value, out = None, strides = None, check_border = True, processes = cpu_count()):
  """Convolves data with a specified kernel at specific points given by a flat array indx under conditon the value is smaller than a number
    
  Arguments:
        data (array): 3d binary array to convolve
        kernel (array): binary orinteger kernel to convolve
        points (array): list of points to convolve given by the flat array coordinates
        maxValue (int): if result of convolution is smaller then this value return True otherwise False in the out array
        processes (int): number of processors
    
  Returns:
        array: list of results of convolution at specified points
    
  Note:
        cython does not support bools -> use view on uint8 as numpy does
  """
  d = data.reshape(-1, order = 'A');
  if data.dtype == bool:
    d = d.view('uint8');
    
  npts = indices.shape[0];
  
  if out is None:
    out = np.zeros(npts, dtype = bool);
  
  if out.shape[0] != npts:
     raise RuntimeError('The output has not the expected size of %d but %d' % (npts, out.shape[0]));
  
  if out.dtype == bool:
    o = out.view('uint8');
  else:
    o = out;

  if kernel.dtype == bool:
    k = np.array(kernel, 'uint8');
  else:
    k = kernel;
  
  if processes is None:
    processes = cpu_count();
  
  if strides is None:
    strides = np.array(data.strides) / np.array(data.itemsize);
  
  #print d.dtype, strides.dtype, kernel.dtype, o.dtype
  if check_border:
    print d.dtype, strides.dtype, k.dtype, indices.dtype, np.array(max_value).dtype, o.dtype
    code.convolve_3d_indices_if_smaller_than(d, strides, k, indices, max_value, o, processes);
  else:
    code.convolve_3d_indices_if_smaller_than_no_check(d, strides, k, indices, max_value, o, processes);
  
  return out;



def test():
  import numpy as np
  import ClearMap.DataProcessing.ConvolvePointList as cpl
  reload(cpl);
  reload(cpl.code);
  
  #data = np.random.rand(2000,1000,1000) > 0.5;
  data = np.random.rand(100,100,500) > 0.5;
  data[[0,-1],:,:] = 0;
  data[:,[0,-1],:] = 0;
  data[:,:,[0,-1]] = 0;
  pts = np.where(data);
  pts = np.array(pts, dtype = int).T;
  
  #np.save('data.npy', data);
  #np.save('points.npy', pts)
  
  from ClearMap.ImageProcessing.Topology.Topology3D import n6
  n6i = np.array(n6, dtype = int)
  
  import time;
  t0 = time.time();
  result = cpl.convolve_3d(data, n6i, pts, processes=24, check_border = True)
  t1 = time.time();
  print('%f secs' % ((t1-t0)));
  
  import scipy.ndimage
  t0 = time.time();
  good = scipy.ndimage.filters.convolve(np.array(data, dtype = float), np.array(n6, dtype = float), mode = 'constant', cval = 0)
  t1 = time.time();
  print('%f secs' % ((t1-t0)));
  
  x,y,z = pts.T
  if np.all(good[x,y,z].astype(result.dtype) == result):
    print('works')
  else:
    print('error!')
  
  
  ptsi = np.where(data.reshape(-1, order = 'A'))[0];
  t0 = time.time();
  result = cpl.convolve_3d_index(data, n6i, ptsi, processes=24)
  t1 = time.time();
  print('%f secs' % ((t1-t0)));
  
  x,y,z = pts.T
  if np.all(good[x,y,z].astype(result.dtype) == result):
    print('works')
  else:
    print('error!')
  
  ptsi = np.where(data.reshape(-1, order = 'A'))[0];
  t0 = time.time();
  resultc = cpl.convolve_3d_index_if_smaller_than(data, n6i, ptsi, 3, processes=24)
  t1 = time.time();
  print('%f secs' % ((t1-t0)));
  
  np.all(resultc == (result < 3))
  
  
  x,y,z = pts.T
  if np.all(good[x,y,z].astype(result.dtype) == result):
    print('works');
  else:
    print('error!');
    
  #import os
  #os.remove('data.npy');
  #os.remove('points.npy')
