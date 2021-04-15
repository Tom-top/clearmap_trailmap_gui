# -*- coding: utf-8 -*-
"""
Interface to write binary files to a numpy memmap

The interface is based on the numpy library.

Example:
    >>> import os
    >>> import ClearMap.Settings as settings
    >>> import ClearMap.IO.MMP as mmp
    >>> filename = os.path.join(settings.ClearMapPath, 'Test/Data/NPY/points.npy');
    >>> points = mmp.read(filename);
    >>> print points.shape
    (5, 3)

Note:
  For image processing we use [x,y,z] order of arrays. To speed up acessing 
  the z planes memmaps are created in fortran order by default.
"""
__author__    = 'Christoph Kirst <ckirst@rockefeller.edu>'
__license__   = 'MIT License <http://www.opensource.org/licenses/mit-license.php>'


import numpy as np

import ClearMap.IO as io;
import ClearMap.IO.Region as rg;


def read(filename, **args):
  if isinstance(filename, basestring):
    data = np.lib.format.open_memmap(filename, mode = 'r');
    #data = np.load(filename, mmap_mode = 'r');
    #return io.dataToRegion(data, **args);
    return data;


def readData(filename, region = None, memmap = True, z = all):
  if isinstance(filename, basestring):
    data = np.lib.format.open_memmap(filename, mode = 'r');
    #data = np.load(filename, mmap_mode = 'r');
    
    if region is None and z is not all:
      region = rg.Region(z = z, source = filename);
    
    if isinstance(region, rg.Region):    
      data = data[region.slice()];
      #data = io.dataToRegion(data, **args);
      
    if memmap:
      return data;
    else:
      return np.array(data);


def dataType(filename):
  if isinstance(filename, basestring):
    data = np.lib.format.open_memmap(filename, mode = 'r');
    #data = np.load(filename, mmap_mode = 'r');
  else:
    data = filename;
  return data.dtype;



def create(filename, shape, dtype = float, fortran = True):
    if isinstance(shape, (int, long)):
      shape = (shape,);
    shape = tuple(shape);
    memmap = np.lib.format.open_memmap(filename, 'w+', shape = shape, dtype = dtype, fortran_order = fortran);
    return memmap;


def write(filename, data, region = all, fortran = True):
    if io.isFile(filename):
      memmap = np.lib.format.open_memmap(filename, 'r+');
      #memmap = np.load(filename, mmap_mode = 'r+');
    else:
      #handle region stuff
      if isinstance(region, rg.Region):
        shape = region.sourceSize();
      else:
        shape = data.shape;
      
      memmap = np.lib.format.open_memmap(filename, mode = 'w+', shape = shape, dtype = data.dtype, fortran_order = fortran);

    if isinstance(region, rg.Region):
      memmap[region.sourceSlice()] = data;
    else:
      memmap[:] = data[:];
    
    return filename
    
    
def writeData(filename, data, fortran = True):
  shape = data.shape;
  memmap = np.lib.format.open_memmap(filename, mode = 'w+', shape = shape, dtype = data.dtype, fortran_order = fortran);
  memmap[:] = data[:];
  return filename;


def dataSize(filename, x = None, y = None, z = None):
  data = read(filename);
  return data.shape;


def test():    
    """Test MMP module"""
    
    import os, numpy
    import ClearMap.IO.MMP as mmp
    reload(mmp)
    
    fn = os.path.split(mmp.__file__);
    fn = os.path.join(fn[0], 'points.npy');
    
    points = numpy.random.rand(5,3);
    mmp.writeData(fn, points);  
    print "Wrote points to " + fn;
    print "Points:"
    print points
    
    points2 = mmp.readData(fn);
    print "Read points: "
    print points2
    
    print "Difference: " + str(numpy.abs(points-points2).max())
    
    os.remove(fn);

if __name__ == "__main__":
    
    test();
    
