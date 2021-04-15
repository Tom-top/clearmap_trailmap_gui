# -*- coding: utf-8 -*-
"""
Data Region

This module provides basic handling of sub region specifications in data

See :mod:`ClearMap.IO` for details.
"""
__author__    = 'Christoph Kirst <ckirst@rockefeller.edu>'
__license__   = 'MIT License <http://www.opensource.org/licenses/mit-license.php>'
__copyright__ = 'Copyright (c) 2017 by Christoph Kirst, The Rockefeller University, New York City'


import numpy as np

class Region(object):
  """Class to handle sub region specifications, allows nested sub regions"""
  
  def __init__(self, region = all, x = None, y = None, z = None, ndim = None, source = None, sourceShape = None, storeSource = False, dtype = None, memory = None):
    """Region
    
    Arguments:
      region : all, tuple, list or Region
        region specifications
      x,y,z : all, int, tuple or None
        region specifications for specific dimensions
      ndim : None or int
        number of dimensions of the region, inferred from region, x, y, z, source or sourceShape
      source : None or a data source
        data source to wich this region referes to
      sourceShape : tuple
        the shape of the source used in case source not given
      storeSource : bool
        if True, keep a reference to the source
      dtype : dtype or None
        data type of the source, inferred from source if source given
      memory : str or None
        the memory type of the underlying array, 'S' for shared memory, 'M' for a memory map and None for other
        
    Note:
      storeSourc is set to False by default as array pointers durig parallel processing are problematic and might trigger 
      copying / pickling of large arrays
    """
    if isinstance(region, Region):
      self._slice       = region._slice;
      self._source      = region._source;
      self._sourceShape = region._sourceShape;
      self._dtype       = region._dtype;
      self._memory      = region._memory;
      ndim = len(self._slice);
    
    if ndim is None:
      if source is not None:
        ndim = io.ndim(source);
      elif sourceShape is not None:
        ndim = len(sourceShape);
      else:
        if x is not None:
          ndim = 1;
        if y is not None:
          ndim = 2;
        if z is not None:
          ndim = 3;
        if isinstance(region, tuple) or isinstance(region, list):
          ndim = len(region);
        if ndim is None:
          ndim = 1;
    
    if region is all or region is None:
      self._slice = tuple([slice(None) for d in range(ndim)]);
    elif isinstance(region, slice):
      self._slice = (region,); 
    elif isinstance(region, tuple) or isinstance(region, list):
      self._slice = tuple([self._parseSlice(r) for r in region]);
    else:
      self._slice = tuple([slice(None) for d in range(ndim)]);
    
    if x is not None or y is not None or z is not None:
      region = list(self._slice);
      for d,r in enumerate([x,y,z]):
        if r is not None:
          region[d] = self._parseSlice(r);
      self._slice = tuple(region);
    
    self._source = None;
    if storeSource:
      self._source = source;
    
    self._sourceShape  = sourceShape;
    if source is not None:
      self._sourceShape = io.shape(source);
    
    self._dtype = dtype;
    self._memory = memory;
    if source is not None:
      self._dtype = io.dtype(source);
      self._memory = io.memory(source);
    
    #simple consistency test
    if self._sourceShape is not None:
      assert len(self._slice) == len(self._sourceShape), "Inconsistent dimensions %d and %d" % (len(self._slice), len(self._sourceShape));
      
      for s,n,d in zip(self._slice, self._sourceShape, range(len(self._slice))):
        if s.stop is not None:
          assert s.stop <= n, "Region specification larger than data %d > %d in axis %d" % (s.stop, n, d);
  
  
  def slice(self):
    """Return region specifications as slice tuple for use with numpy"""
    return self._slice;
    
  def source(self):
    """Returns ultimate source of the region"""
    if isinstance(self._source, Region):
      return self._source.source();
    else:
      return self._source;
  
  def ndim(self):
    """Returns dimension of the region"""
    return len(self._slice);
    
  def dtype(self):
    """Returns data type of the region"""
    return self._dtype;
  
  def memory(self):
    """Returns memory type of the region"""
    return self._memory; 
  
  def shape(self, axis = None, undefined = all):
    """Tries to infer size of the region, returns undefined for undefined dimensions"""
    if axis is not None:
      s = self._slice[axis];
      
      stop = s.stop;
      if stop is None:
        if self._sourceShape is None:
          return undefined;
        else:
          stop = self._sourceShape[axis];
        
      start = s.start;
      if start is None:
        start = 0;
      
      step = s.step;
      if step is None:
        step = 1;
  
      return ((stop-1) - start)/step + 1;  
    
    else:     
      return tuple([self.shape(axis = d, undefined = undefined) for d in range(self.ndim())]);
  
  def sourceShape(self, depth = all):
    """Returns the shape of the underlying source"""
    if depth is all:
      depth = np.inf;
     
    if isinstance(self._source, Region) and depth > 0:
        return self._source.sourceShape(depth = depth -1);
    else:
      return self._sourceShape;
  
  def sourceSlice(self, depth = all):
    """Calculates the region in source data as slices"""
    if not isinstance(self._source, Region) or depth == 0:
      return self._slice
    
    if depth is all:
      depth = np.inf;
    
    superslice = self._source.sourceSlice(depth = depth -1);
    ndim = self.ndim();
    thisslice = self._slice;
    
    newslice = [];
    for d in range(ndim):
      r = thisslice[d];
      sr = superslice[d];

      sstep = sr.step;
      if sstep is None:
        sstep = 1;
      rstep = r.step;
      if rstep is None:
        rstep = 1;
      step = rstep * sstep;
      if step == 1:
        step = None;
          
      if r.start is None:
        start = sr.start;
      else:
        if sr.start is None:
          start = r.start;
        else:
          start = sstep * r.start + sr.start;
      
      if r.stop is None:
        stop = sr.stop;
      else:
        if sr.start is None:
          stop = r.stop * sstep;
        else:
          stop = (sr.start + r.stop) * sstep;
      
      newslice.append(slice(start, stop, step));
    
    return tuple(newslice);
  
  def read(self, source = None, **args):
    if source is None:
      source = self.source();
    return io.read(source, region = self, **args);

  def readData(self, **args):
     return self.read(**args);
  
  def write(self, data, sink = None, **args):
    if sink is None:
      sink = self.source();
    return io.write(data, sink, region = self, **args);
    
  def writeData(self, **args):
     return self.write(**args);
  
  
  def _parseSlice(self, r):
    """Parses a x,y,z argument to a slice object"""
    if isinstance(r,tuple) or isinstance(r,list):
      return slice(*r);
    elif isinstance(r,int):
      return slice(r);
    else:
      return slice(None);  
  
  def _strSlice(self, arg):
    """Formats printin of a slice object"""
    if isinstance(arg, slice):
      #if arg.start is None and arg.stop is None and arg.step is None:
      #  return None;
      return (arg.start, arg.stop, arg.step);
    if isinstance(arg, list):
      return str([self._strSlice(a) for a in arg]).replace(' ', '')
    if isinstance(arg, tuple):
      return str(tuple([self._strSlice(a) for a in arg])).replace(' ', '');
    return str(arg);
    
  def __str__(self):
    return "<Region: %s>" % self._strSlice(self.slice());
    
  def __repr__(self):
    return self.__str__() 


import ClearMap.IO.IO as io



if __name__ == '__main__':
  import numpy as np;
  import ClearMap.IO.Region as reg
  
  data = np.random.rand(40,50);
  
  r = reg.Region(source = data, x = (7,10));
    
 

 
#
#def readData(region, **args):
#  """Read the data from a region"""
#  return region.readData(**args);
#
#  
#def writeData(region, data, **args):
#  """write data to a specified region"""
#  return region.writeData(data, **args);
  
  #def dataSize(self):
  #  sourceSize = self.sourceSize();
  #  sourceSlice = self.sourceSlice();
    
  
#def _parseShape(shape, dim):
#  if shape is all:
#    return tuple([shape for d in range(dim)]);
#  if isinstance(shape, tuple) or isinstance(shape, list):
#    return tuple([shape[d % len(shape)] for d in range(dim)]);
#
#
#def cropRegionToShape(region, shape = all):
#  shape = _parseShape(shape, region.dim());
#  for s,r in zip(shape, region.slice()):
#    if s is not all or s is not None:
#      
#      
#def sliceFromRegion(region, shape = all):
#  if isinstance(region, slice):
#    return cropRegionToShape(region, shape);
#  elif isinstance(region, tuple) or isinstance(region, list):
#    return tuple([sliceFromRegion(r, shape)])


#class Region():
#  """SourceRegion class specifies rectangular sub regions of a data source"""
#  
#  def __init__(self, source = None, region = all, x = None, y = None, z = None, parent = None):
#    """Region class constructor
#        
#    Arguments:
#      source (str, array, tuple or None): the data source or tuple of data size to determine the size of the refrence data
#      region (tuple of tuples or all): region of the form ((xmin, xmax), (ymin, ymax), ...) or ((xmin,xmax,dx),...)
#      x,y,z (tuple optional): specific ranges for x,y,z coordinates
#      parent (Region): a parent region
#    """
#    self._region = slice(None);
#    
#    if isinstance(source, tuple):
#      self._sourceSize = source;
#      self._parent = None;
#      self._region = self.fromArgs(region = region, x = x, y = y, z = z);
#    elif isinstance(region, Region):
#      self._sourceSize = region._sourceSize;
#      self._parent = region._parent;
#      self._region = self.fromArgs(region = region._region, x = x, y = y, z = z);
#    elif io.isSource(source):
#      if parent is not None:
#        raise RuntimeError('source and parent specified simultaneously');
#      self._sourceSize = io.size(source);
#      self._parent = None;
#      self._region = self.fromArgs(region = region, x = x, y = y, z = z);
#    elif isinstance(parent, Region):
#      self._sourceSize = parent._sourceSize();
#      self._parent = parent;
#      self._region = self.fromArgs(region = region, x = x, y = y, z = z);
#    else:
#      raise RuntimeError('underlying size of reference data for region cannot be determined!');
#  
#  
#  def dim(self):
#    """Dimension of the underlying data source"""
#    if self._parent is not None:
#      return self._parent.dim();
#    elif self._sourceSize is not None:
#      return len(self._sourceSize);    
#    elif self._region is not None:
#      return len(self._region);
#    else:
#      raise RuntimeError('dimension of underlying data source in a region cannot be determined!');
#  
#  def sourceSize(self):
#    """Returns full size of the undelying data source"""
#    if self._sourceSize is not None:
#      return (self._sourceSize);
#    elif self._parent is not None:
#      return self._parent.sourceSize();
#  
#  def fromArgs(self, region = all, x = None, y = None, z = None, **args):
#    """Parses arguments to a region list [(xmin,xmax), (ymin, ymax), ...]
#    
#    Arguments:
#      region (list of tuples or all): region list of the form [(xmin, xmax), (ymin, ymax), ...]
#      x,y,z (tuple optional): specific ranges for x,y,z coordinates overwrites region specifications
#    """
#    
#    dim = self.dim();
#    
#    if region is all or region is None:
#      region = [(all,all) for d in range(dim)];
#    if (isinstance(region, tuple) or isinstance(region, list)) and len(region) == dim:
#      region = list(region);
#      for d in range(dim):
#        if region[d] is all:
#          region[d] = (all, all)
#        elif not isinstance(region[d], tuple) and len(region[d], 2):
#          raise RuntimeError('region specification at dimension %d not all or tuple of length 2!' % d);
#    else:
#        raise RuntimeError('region is not list of dimensions %d!' % dim);
#    #print region
#
#    for d in range(dim):
#      for i in range(2):
#        if not (isinstance(region[d][i], int) or region[d][i] is all):
#          raise RuntimeError('region at entry %d,%d is not a int or all!' % (d,i));
#    
#    for (i, (xyz, label)) in enumerate(zip((x,y,z), ('x', 'y', 'z'))):
#      if xyz is not None:
#        if dim <= i:
#           raise RuntimeError('region specification %s=%s invalid as dimension is <= %d!' % (label, str(xyz),dim));
#        if xyz is all:
#          region[i] = (all, all);
#        elif isinstance(xyz, tuple) and len(xyz) == 2:
#          region[i] = xyz;
#        else:
#          raise RuntimeError('region specification %s=%s is not a int or all!' % (label, str(xyz)));
#    
#    return tuple(region);
#  
#  
#  def size(self):
#    """Returns size of the region"""
#    region = self.absoluteRegion();
#    return tuple([r[1] - r[0] for r in region]);
#  
#  
#  def parent(self):
#    """Returns parent of this region"""
#    return self._parent;
#  
#  
#  def region(self):
#    """Returns region with potential all place holders"""
#    return self._region;
#  
#  
#  def setSourceSize(self, sourceSize):
#    """Sets the source size for this region
#    
#    Arguments:
#      sourcesize (str or tuple): the source size or a source
#    """
#    if io.isSource(sourceSize):
#      sourceSize = io.size(sourceSize);
#    
#    self._sourceSize = sourceSize;
#  
#  
#  def setParent(self, parent):
#    """Sets the parent region for this region
#    
#    Arguments:
#      parent (None or Region): the parent for this region
#    """
#    self._parent = parent;
#    if parent is not None:
#      self._sourceSize = parent.sourceSize();
#  
#  
#  def setRegion(self, region = None, x = None, y = None, z = None):
#    """Set the region
#    
#    Arguments:
#      region (list of tuples or all): region list of the form [(xmin, xmax), (ymin, ymax), ...]
#      x,y,z (tuple optional): specific ranges for x,y,z coordinates overwrites region specifications
#    """
#    if region is None:
#      region = self._region;
#    
#    self._region = self.fromArgs(region = region, x = x, y = y, z = z);
#  
#  
#  def absoluteRegion(self):
#    """Calculates the absolute region in source data as ints replacing all 'all' entries"""
#
#    dim = self.dim();
#    region = self._region;
#    
#    if self._parent is None:
#      superregion = None;
#      size = self.sourceSize();
#    else:
#      superregion = self._parent.absoluteRegion();
#      size = tuple([r[1] - r[0] for r in superregion]);
#
#    rg = numpy.zeros((dim,2), dtype = int);
#    for d in range(dim):
#      for i in range(2):
#        r = region[d][i];
#        if r is all:
#          if i == 0:
#            r = 0;
#          else:
#            r = size[d];
#            
#        if isinstance(r, int):
#          if r < 0:
#            r = size[d] + r;
#          if r < 0 or r > size[d]:
#            raise RuntimeError('region %d,%d = %d out of range!' % (d,i,r));
#        else:
#          raise RuntimeError('region at entry %d,%d is not a int or all!' % (d,i));
#        
#        rg[d][i] = r;
#    
#    if superregion is not None:
#      for d in range(dim):
#        for i in range(2):
#          rg[d][i] += superregion[d][0];
#    
#    region = tuple([tuple(rr) for rr in rg]);
#    return region;
#  
#  
#  def data(self, source):
#    """Returns data in the specified region
#    
#    Arguments:
#      source (str or array or None): source of the data
#    
#    Returns:
#      array: the data in the specified region
#    """
#  
#    if source is None:
#      return None;
#    elif isinstance(source, numpy.ndarray):
#      r = self.absoluteRegion();
#      data = source;
#      for d in range(self.dim()):
#        data = data[r[d][0]:r[d][1]];
#      return data;
#    else:
#      return io.read(source, region = self);
#  
#  
#  def flatten(self):
#    """Flattens the region hierarchy"""
#    if self._parent is not None:
#      self._region = self.absoluteRegion();
#      self._sourceSize = self._parent.sourceSize();
#      self._parent = None;
#  
#  
#  def coordinateShift(self):
#      """Calculate shift of coordinates from original origin of the source to origin of region
#      
#      Returns:
#          list: shift of points from original origin of data source to origin of region reduced data
#      """      
#      
#      return [r[0] for r in self.absoluteRegion()];
#
#  
#  def coordiantes(self, coordiantes, shift = False, data = None):
#      """Restrict array of point coordinates to this region
#      
#      Arguments:
#          coordiantes (array or str): point source
#          shift (bool): shift points to relative coordinates in the reduced image
#          data (array): optional array of meta data attached to the points to get reduced simultaneously
#      
#      Returns:
#          array or tuple: coordinates reduced in region and optionally shifted to the region's coordinates
#                          if data is not None then return tuple (reduced coordinates, reduced data)
#      """      
#      
#      if coordiantes is None:
#          return None;
#      elif isinstance(coordiantes, basestring):
#          coordinates = io.read(coordiantes);
#      
#      if isinstance(coordinates, numpy.ndarray):
#          n = coordinates.shape[1];
#          ra = self.absoluteRegion();
#          r = self.region();
#
#          ids = numpy.ones(n, dtype = 'bool');          
#          for d in range(self.dim()):
#            if r[d][0] is not all:
#              ids = numpy.logical_and(ids, coordinates[:,d] >= ra[d][0]);
#            if r[d][1] is not all:
#              ids = numpy.logical_and(ids, coordinates[:,d] <  ra[d][1]);
#
#          coordinates = coordinates[ids, :];
#              
#          if shift:
#              sh = [rr[0] for rr in ra];
#              coordinates = coordinates - sh;
#          
#          if data is None:
#              return coordinates;
#          else:
#              return (coordinates, data[ids]);
#          
#      else:
#          raise RuntimeError('coordiantes not None, str or numpy array!');
#
#
#
#
#
#def sizeFromRegion(source = None, region = None, **args):
#    """Converts data size to actual size given region specifications
#    
#    Arguments:
#        source (str, array or tuple): data source or 
#        region (tuple or Region): region specification
#        x,z,y (tuple or all): range specifications, ``all`` is full range
#        
#    Returns:
#        tuple: data size as tuple of integers
#    
#    See Also:
#        :func:`regionFromRegion`, :func:`dataFromRegion`
#    """
#    
#    if isinstance(source, Region):
#      region = source;
#      source = None;
#
#    if source is None and not isinstance(region, Region):
#       raise RuntimeError('need to specify region or source to determine region size');
#    
#    r = Region(source = source, region = region, **args);
#    
#    return r.size();
#    
#  
#def regionFromRegion(source = None, region = None, **args):
#    """Converts region into absolute region specifications
#    
#    Arguments:
#        source (str, array or tuple): data source or 
#        region (tuple or Region): region specification
#        x,z,y (tuple or all): range specifications, ``all`` is full range
#        
#    Returns:
#        tuple: region as tuple of tuple of integers
#    
#    See Also:
#        :func:`sizeFromRegion`, :func:`dataFromRegion`
#    """
#    
#    if isinstance(source, Region):
#      region = source;
#      source = None;
#
#    if source is None and not isinstance(region, Region):
#       raise RuntimeError('need to specify region or source to determine region');
#    
#    r = Region(source = source, region = region, **args);
#    
#    return r.absoluteRegion();
#
#
#def region(source = None, region = None, **args):
#    """returns standarized region from region specifications
#    
#    Arguments:
#        source (str, array or tuple): data source or 
#        region (tuple or Region): region specification
#        x,z,y (tuple or all): range specifications, ``all`` is full range
#        
#    Returns:
#        tuple: region as tuple of tuple of integers
#    
#    See Also:
#        :func:`sizeFromRegion`, :func:`dataFromRegion`
#    """
#    
#    if isinstance(source, Region):
#      region = source;
#      source = None;
#
#    if source is None and not isinstance(region, Region):
#       raise RuntimeError('need to specify region or source to determine region');
#    
#    r = Region(source = source, region = region, **args);
#    
#    return r.region();
#
#  
#def dataFromRegion(source, region = None, **args):
#    """Reduces data to specified ranges
#    
#    Arguments:
#        source (str, array or tuple): data source or 
#        region (tuple or Region): region specification
#        x,z,y (tuple or all): range specifications, ``all`` is full range
#        
#    Returns:
#        array: reduced data
#    
#    See Also:
#        :func:`dataSizeFromDataRange`
#    """  
#
#    if source is None:
#      return None;
#    
#    r = Region(source = source, region = region, **args);
#    
#    return r.data(source);    
#
#
#def coordinateShiftFromRegion(source = None, region = None, **args):
#    """Calculate shift of points given a specific range restriction
#    
#    Arguments:
#        source (str, array or tuple): data source or 
#        region (tuple or Region): region specification
#        x,z,y (tuple or all): range specifications, ``all`` is full range
#    
#    Returns:
#        list: shift of points from original origin of data to origin of range reduced data
#    """      
#    
#    if isinstance(source, Region):
#      region = source;
#      source = None;
#    
#    if source is None and not isinstance(region, Region):
#       raise RuntimeError('need to specify region or source to determine coordinate shifts');
#    
#    r = Region(source = source, region = region, **args);
#    
#    return r.coordinateShift();    
#
#
#
#
#def coordiantesFromRegion(coordinates, source = None, region = None, shift = False, data = None, **args):
#    """Restrict coordinates to a specific region
#    
#    Arguments:
#        coordinates (array or str): array source of coordinates
#        source (str, array or tuple): data source or data size specifications
#        region (tuple or Region): region specification
#        x,z,y (tuple or all): range specifications, ``all`` is full range
#        shift (bool): shift points to relative coordinates in the region
#        data (array): addiational meta data for each points to reduce as well
#    
#    Returns:
#        array: coordinates reduced in range and optionally shifted to the range reduced origin
#        or (array, array): reduced coordinates and data
#    """      
#
#    if coordinates is None:
#      return None;
#    
#    if isinstance(source, Region):
#      region = source;
#      source = None;
#      
#    if source is None and region is None:
#      source = tuple(coordinates.max(axis=0));
#    
#    if source is None and not isinstance(region, Region):
#      raise RuntimeError('need to specify region or source to determine coordinates within region');
#    
#    r = Region(source = source, region = region, **args);
#    
#    return r.coordinates(coordinates, shift = shift, data = data);
#
#
#
#def test():
#  import ClearMap.IO.Region as r;
#  pass
#
#
#if __name__ == '__main__':
#  test();
#      
