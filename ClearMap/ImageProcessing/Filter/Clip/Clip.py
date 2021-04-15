"""
Module to compute clipped images

Usefull to sace memory in large data sets
"""
__author__    = 'Christoph Kirst <ckirst@rockefeller.edu>'
__license__   = 'MIT License <http://www.opensource.org/licenses/mit-license.php>'
__copyright__ = 'Copyright (c) 2017 by Christoph Kirst, The Rockefeller University, New York City'


import numpy as np

import pyximport;
pyximport.install(setup_args={"include_dirs":np.get_include()}, reload_support=True)

import ClearMap.ImageProcessing.Filter.Clip.ClipCode as code


def clip(source, sink = None, clipmin = None, clipmax = None, norm = None, sink_dtype = None, processes = 1):
    """Clip and normalize image

    Parameters
    ----------
    source : 3-D array
        Input source.
    sink : 3-D array
        If None, a new array is allocated.
    clipmin : number
        Minimal number to clip to.
    clipmax : number
        Maximal number to clip to.
    norm : number
        Normalization constant.

    Returns
    -------
    sink : 3-D array
        clipped output.
    """

    if source.ndim != 3:
      raise ValueError('source assumed to be 3d found %dd' % source.ndim);
    
    if source is sink:
        raise NotImplementedError("Cannot perform operation in place.")

    if clipmin is None:
      clipmin = 0;
    if clipmax is None:
      clipmax = 255;
    if norm is None:
      norm = 255;

    #TODO: homogenize this in all processing functions: the info can be in sink -> single routine to creatre a new array from sink / source specs: e.g. order, shared mem etc. -> io
    if sink is None:
      if sink_dtype is None:
          sink_dtype = source.dtype
      
      order = 'C'; 
      if np.isfortran(source):
        order = 'F';
      sink = np.empty(source.shape, dtype = sink_dtype, order = order)
    
    if sink.dtype == bool:
      s = sink.view('uint8')
    else:
      s = sink;
      
    code.clip(source, s, clipmin, clipmax, norm, processes);
    
    return sink;
  