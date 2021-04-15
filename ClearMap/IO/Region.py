# -*- coding: utf-8 -*-
"""
Data Region

This module provides basic handling of sub region specifications in data

See :mod:`ClearMap.IO` for details.
"""
#:copyright: Copyright 2015 by Christoph Kirst, The Rockefeller University, New York City
#:license: GNU, see LICENSE.txt for details.

#TODO: conversion between physical space and regions

import numpy as np
import ClearMap.IO.IO as io


def dataSize(source):
  if isinstance(source, Region):
    return source.size();
  else:
    return io.dataSize(source);


class Region(object):
  """Class to handle sub region specifications, allows nested sub regions"""
  
  def __init__(self, region = all, source = None, x = None, y = None, z = None, dim = None):
    """Constructor"""
    if dim is None:
      if source is not None:
        dim = len(dataSize(source));
      else:
        if x is not None:
          dim = 1;
        if y is not None:
          dim = 2;
        if z is not None:
          dim = 3;
        if isinstance(region, tuple) or isinstance(region, list):
          dim = len(region);
        if dim is None:
          dim = 1;
    
    if region is all or region is None:
      self._region = tuple([slice(None) for d in range(dim)]);
    elif isinstance(region, slice):
      self._region = (region,); 
    elif isinstance(region, tuple) or isinstance(region, list):
      sldat = [];
      for r in region:
        sldat.append(self._createSlice(r));
      self._region = tuple(sldat);
    else:
      self._region = tuple([slice(None) for d in range(dim)]);
    
    region = list(self._region);
    for d,r in enumerate([x,y,z]):
      if r is not None:
        region[d] = self._createSlice(r);
    
    self._region = tuple(region);
      
    self._source = source;
  
  
  def _createSlice(self, r):
    if isinstance(r,tuple) or isinstance(r,list):
      return slice(*r);
    elif isinstance(r,int):
      return slice(r);
    elif isinstance(r, slice):
      return r;
    else:
      return slice(None); 
    
  
  def slice(self):
    """Return region specifications as slice tuple for use with numpy"""
    return self._region;
    
  def dim(self):
    """Returns dimension of the region"""
    return len(self._region);
  
  def source(self):
    """Returns source of the region"""
    if isinstance(self._source, Region):
      return self._source.source();
    else:
      return self._source;
    
  def sourceSize(self):
    """Returns the size of the underlying source"""
    source = self.source();
    return io.dataSize(source);
    
  def sourceShape(self):
     return self.sourceSize();
    
  def size(self):
    ssize = self.sourceSize();
    sslice = self.sourceSlice();
    sdat = [];
    for r,q in zip(sslice, ssize):
      s = r.start;
      e = r.stop;
      d = r.step;
      if s is None:
        s = 0;
      if e is None:
        e = q;
      if d is None:
        d = 1;
      sdat.append(int(np.ceil((e-s) * 1.0/d)));
    return tuple(sdat);
    
  def shape(self):
      return self.size();
  
  def sourceSlice(self, simplify = False):
    """Calculates the region in source data as slices"""
    if not isinstance(self._source, Region):
      return self._region
    
    superregion = self._source.sourceSlice();
    dim = self.dim();
    region = self._region;
    sourceShape = self.sourceShape();
    
    sldat = [];
    for d in range(dim):
      r = region[d];
      sr = superregion[d];

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
      
      if simplify:
        if start == 0:
          start = None;
          
        if stop == sourceShape[d]:
          stop = None;
      
      sldat.append(slice(start, stop, step));
    
    return tuple(sldat);
  
  def readData(self, **args):
    source = self.source();
    sslice = self.sourceSlice();
    return io.readData(source, sslice, **args);

  def read(self, **args):
     return self.readData(**args);
  
  def writeData(self, data, **args):
    source = self.source();
    return io.writeData(source, data, **args);
    
  def write(self, **args):
     return self.writeData(**args);
  
  def dataType(self):
    """Returns data type of the source"""
    return io.dataType(self.source());
    
    
  def shift(self, shift = None, x = None, y = None, z = None):
    """Shifts the region"""
    if shift is not None:
      x,y,z = shift;
      
    sshape = self.sourceShape(); 
    
    reg = list(self._region);
    for d,s in enumerate((x,y,z)):
      if s is not None and s != 0:
        r = self._region[d];
        start, stop = r.start, r.stop;
        if start is None:
          start = max(min(s, sshape[d]),0);
        else:
          start = max(min(start + s, sshape[d]),0);
        if stop is None:
          stop  = max(min(sshape[d] + s, sshape[d]),0);
        else:
          stop  = max(min(stop + s, sshape[d]),0);
        
        if start == 0:
          start = None;
        if stop == sshape[d]:
          stop = None;
        
        reg[d] = slice(start, stop, r.step);
    
    self._region = tuple(reg);
  
  def boundingBox(self, dim = None):
    """Returns bounding box in source"""
    ss = self.sourceSlice();
    sh = self.sourceShape();
    
    bb = [];
    for d,s in enumerate(ss):
      if s.start is None:
        s0 = 0;
      else:
        s0 = s.start;
        
      if s.stop is None:
        s1 = sh[d];
      else:
        s1 = s.stop;
      
      bb.append([s0,s1]);
    
    if dim is None:
      return bb;
    else:
      return bb[dim];
        
  def lowerBound(self):
    return [b[0] for b in self.boundingBox()];
    
  def upperBound(self):
    return [b[1] for b in self.boundingBox()];
  
  def __str__(self):
    return "<Region: %s>" % str(self.slice());
    
  def __repr__(self):
    return self.__str__() 



def readData(region, **args):
  """Read the data from a region"""
  return region.readData(**args);

  
def writeData(region, data, **args):
  """write data to a specified region"""
  return region.writeData(data, **args);
  

if __name__ == '__main__':
  import numpy as np
  import ClearMap.IO.Region as reg;
  reload(reg)
  
  data = np.random.rand(10,20,30);
  
  r = reg.Region(source = data, z = (5,20))
  print(r)
  
  r.shift(z = 10)
  print(r)

  r.shift(z = -20)
  print(r)





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



