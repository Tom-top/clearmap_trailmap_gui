"""
Module to compute generic filters based on Hessian Matrix

Usefull for filtering vasculature data
"""
__author__    = 'Christoph Kirst <ckirst@rockefeller.edu>'
__license__   = 'MIT License <http://www.opensource.org/licenses/mit-license.php>'
__copyright__ = 'Copyright (c) 2017 by Christoph Kirst, The Rockefeller University, New York City'


import numpy as np

import pyximport;
pyximport.install(setup_args={"include_dirs":np.get_include()}, reload_support=True)


from . import HessianCode as code


__all__ = ['tubeness'];


def _handle_input(source, sink, sink_dtype=None, pixel_size=1):

    assert len(source.shape) == 3
    
    source = np.ascontiguousarray(source)
    
    if source is sink:
        raise NotImplementedError("Cannot perform operation in place.")

    if sink is None:
        if sink_dtype is None:
            sink_dtype = source.dtype
        sink = np.empty(source.shape+(pixel_size,), dtype = sink_dtype)
    else:
        if len(sink.shape) == 3:
            sink = sink.reshape(sink.shape+(pixel_size,))
    
    if sink.dtype == bool:
      sink = sink.view('uint8')
      is_bool = True;
    else:
      is_bool = False;
    
    return source, sink, is_bool


def _apply_scalar_per_pixel(func, source, sink, sink_dtype = None, par = 0):

    source, sink, is_bool = _handle_input(source, sink, sink_dtype)

    func(source, sink = sink, par = par)

    if is_bool:
      sink = sink.view('bool');

    return sink.reshape(sink.shape[:3])


#def _apply_vector_per_pixel(func, source, sink, sink_dtype=None, pixel_size=1, par = 0):
#
#    source, sink, is_bool = _handle_input(source, sink, sink_dtype, pixel_size=pixel_size)
#
#    func(source, sink = sink, par = par)
#
#    if is_bool:
#       sink = sink.view('bool');
#
#    return sink


def tubeness(source, sink = None, threshold = None):
    """Tubness mesure of source

    Parameters
    ----------
    srouce : 3-D array
        Input source.
    sink : 3-D array
        If None, a new array is allocated.
    threshold : float or None
        If not None the tubeness is thresholded at this level

    Returns
    -------
    sink : 3-D array
        Tubness output.
    """

    if threshold is None:
      return _apply_scalar_per_pixel(code.tubeness, source, sink)
    else: 
      return _apply_scalar_per_pixel(code.tubeness_threshold, source, sink, sink_dtype = bool, par = threshold)