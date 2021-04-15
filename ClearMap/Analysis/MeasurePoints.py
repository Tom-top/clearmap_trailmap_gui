#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
MeasurePoints module

Measure properites at specified points
"""

import numpy as np

import ClearMap.IO as io

#import ClearMap.DataProcessing.LargeData as ld
import ClearMap.DataProcessing.ConvolvePointList as cpl


def measure(data, kernel = None, points = None, indices = None, x = None, y = None, z = None, out = None, out_dtype = None, strides = None, check_border = True, processes = None):
  
  if isinstance(data, str):
    d = io.readData(data);
  else:
    d = data;
    
  if isinstance(points, str):
    p = io.readPoints(points);
  else:
    p = points;
  
  if kernel is None:
    kernel = np.ones((1,1,1), dtype = int);
  
  return cpl.convolve_3d(data = data, kernel = kernel, points = points, indices = indices, x = x, y = y, z = z, out = out, out_dtype = out_dtype, strides = strides, check_border = check_border, processes = processes);
