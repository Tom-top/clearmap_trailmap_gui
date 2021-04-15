# -*- coding: utf-8 -*-
"""
Basic IO classes
"""

#io takes source and wraps class

from IO.MMP import MMP
from IO.NPY import NPY
from IO.TIF import TIF

sources = []


def sourceType(sorce):
  if isinstance(source, str):
    #check if file or file expression
  elif isinstance(source, np.ndarray):
    
  elif

def sourceInfo(source):
  soure = io.source(source);
  return source.shape, source.dtype, 
  
  

class Info(object):
  """Class provinding information on data"""
  
  def __init__(self, source = None, shape = None, dtype = None, order = None, memory = None, location = None):
    self._shape = shape;
    self._dtype = dtype;
    self._order = order;
    self._memory = memory
    self._location = location;
    
    if source is not None:
      self._shape, self._dtype = io.sourceInfo(source);


class Data(DataInfo):
  """Base class for all data"""
  
  def __init__(self, data):
    self.data = data;
  
  def read(self, region = all, sink = None):
    if sink is None:
      return self.data[region.slice()];
    else:
      sink.
      
  def write(self, data, region = all):
    if 
    self.data = data[region.slice()];