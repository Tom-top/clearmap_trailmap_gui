#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Specialized ColorMaps

Includes color maps for 3d orientation fields

Code based on dipy/viz/colormap.py
"""

import numpy as np

import vispy as vp
import vispy.color as vpc

import matplotlib as mpl
import matplotlib.pyplot as plt


### colors

def color(arg, alpha = None):
  """Converts color arguent to rgb values"""
  try:
    c = vpc.Color(arg);
    if alpha is False:
      c = c.rgb;
    else:
      c = c.rgba;
      if alpha is not None and alpha is not True:
        c[-1] = alpha;
    return c;
  except:
    pass;
  
  try:
    if alpha is False:
      c = mpl.colors.to_rgb(arg);
    else:
      c = mpl.colors.to_rgba(arg);
      if alpha is not None and alpha is not True:
        c[-1] = alpha;    
    return c;
  except:
    raise RuntimeError('color %r not recognized' % arg);



### color maps

def colormap(cm = None, alpha = True):
  """Convert argument to a color map function"""
  
  if cm is None:
    return colormap('gray');
  
  cmps = vpc.get_colormaps()
  if cm in cmps.keys():
    cm = vpc.colormap.get_colormap(cm);
    
  if isinstance(cm, vp.color.BaseColormap):
    if alpha is False:
      def cmap(x):
        return cm.map(np.array(x))[:,:3];
    elif alpha is True or alpha is None:
      def cmap(x):
        return cm.map(np.array(x));
    else:
      def cmap(x):
        c = cm.map(np.array(x));  
        c[:,-1] = alpha;
        return c;
    return cmap


  if cm in plt.cm.cmap_d.keys():
    cm = plt.get_cmap(cm);
    
  if isinstance(cm, mpl.colors.Colormap):
    if alpha is False:
      def cmap(x):
        return cm(np.array(x))[:,:3];
    elif alpha is True or alpha is None:
      def cmap(x):
        return cm(np.array(x));
    else:
      def cmap(x):
        c = cm(np.array(x));  
        c[:,-1] = alpha;
        return c;
    return cmap
    

  if isinstance(cm, (list, tuple, np.ndarray)):
    cols = [color(a) for a in cm];
    cm = vpc.ColorMap(cols, interpolation = 'linear');
    return colormap(cm, alpha = alpha);
  
  else:
    raise RuntimeError('cannot find colormap: %r' % cm);
  
  
  
### special color maps


def _cc(na, nd):
  return (na * np.cos(nd * np.pi / 180.0))

def _ss(na, nd):
  return na * np.sin(nd * np.pi / 180.0)

def boys2rgb(v):

    x = v[0]
    y = v[1]
    z = v[2]

    x2 = x ** 2
    y2 = y ** 2
    z2 = z ** 2

    x3 = x * x2
    y3 = y * y2
    z3 = z * z2

    z4 = z * z2

    xy = x * y
    xz = x * z
    yz = y * z

    hh1 = .5 * (3 * z2 - 1) / 1.58
    hh2 = 3 * xz / 2.745
    hh3 = 3 * yz / 2.745
    hh4 = 1.5 * (x2 - y2) / 2.745
    hh5 = 6 * xy / 5.5
    hh6 = (1 / 1.176) * .125 * (35 * z4 - 30 * z2 + 3)
    hh7 = 2.5 * x * (7 * z3 - 3 * z) / 3.737
    hh8 = 2.5 * y * (7 * z3 - 3 * z) / 3.737
    hh9 = ((x2 - y2) * 7.5 * (7 * z2 - 1)) / 15.85
    hh10 = ((2 * xy) * (7.5 * (7 * z2 - 1))) / 15.85
    hh11 = 105 * (4 * x3 * z - 3 * xz * (1 - z2)) / 59.32
    hh12 = 105 * (-4 * y3 * z + 3 * yz * (1 - z2)) / 59.32

    s0 = -23.0
    s1 = 227.9
    s2 = 251.0
    s3 = 125.0

    ss23 = _ss(2.71, s0)
    cc23 = _cc(2.71, s0)
    ss45 = _ss(2.12, s1)
    cc45 = _cc(2.12, s1)
    ss67 = _ss(.972, s2)
    cc67 = _cc(.972, s2)
    ss89 = _ss(.868, s3)
    cc89 = _cc(.868, s3)

    X = 0.0

    X = X + hh2 * cc23
    X = X + hh3 * ss23

    X = X + hh5 * cc45
    X = X + hh4 * ss45

    X = X + hh7 * cc67
    X = X + hh8 * ss67

    X = X + hh10 * cc89
    X = X + hh9 * ss89

    Y = 0.0

    Y = Y + hh2 * -ss23
    Y = Y + hh3 * cc23

    Y = Y + hh5 * -ss45
    Y = Y + hh4 * cc45

    Y = Y + hh7 * -ss67
    Y = Y + hh8 * cc67

    Y = Y + hh10 * -ss89
    Y = Y + hh9 * cc89

    Z = 0.0

    Z = Z + hh1 * -2.8
    Z = Z + hh6 * -0.5
    Z = Z + hh11 * 0.3
    Z = Z + hh12 * -2.5

    # scale and normalize to fit
    # in the rgb space

    w_x = 4.1925
    trl_x = -2.0425
    w_y = 4.0217
    trl_y = -1.8541
    w_z = 4.0694
    trl_z = -2.1899

    C = np.zeros((3,) + x.shape);

    C[0] = 0.9 * np.abs(((X - trl_x) / w_x)) + 0.05
    C[1] = 0.9 * np.abs(((Y - trl_y) / w_y)) + 0.05
    C[2] = 0.9 * np.abs(((Z - trl_z) / w_z)) + 0.05


    ids = np.logical_and(x == 0, np.logical_and(y == 0, z == 0));
    C[:,ids] = 0;

    return C;


def orient2rgb(v):
    """ standard orientation 2 rgb colormap
    v : array, shape (N, 3) of vectors not necessarily normalized
    Returns
    -------
    c : array, shape (N, 3) matrix of rgb colors corresponding to the vectors
           given in V.
    Examples
    --------
    >>> from dipy.viz import colormap
    >>> v = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    >>> c = colormap.orient2rgb(v)
    """

    if v.ndim == 1:
        orient = v
        orient = np.abs(orient / np.linalg.norm(orient))

    if v.ndim == 2:
        orientn = np.sqrt(v[:, 0] ** 2 + v[:, 1] ** 2 + v[:, 2] ** 2)
        orientn.shape = orientn.shape + (1,)
        orient = np.abs(v / orientn)

    return orient





### rgb look up tables

#create direct grid as 256**3 x 4 array 
def rgbLUT():
    xl = np.mgrid[0:256, 0:256, 0:256]
    lut = np.vstack((xl[0].reshape(1, 256**3),
                     xl[1].reshape(1, 256**3),
                     xl[2].reshape(1, 256**3),
                     255 * np.ones((1, 256**3)))).T
    return lut.astype('int32')

# indexing function to above grid
def rgbToLUT(rgb):
    r,g,b = rgb;
    return 256**2 * r + 256 * g + b

  
