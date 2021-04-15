# -*- coding: utf-8 -*-
"""
3d Skeletonization

Main routines for 3d skeletonization.

Supported algorithsm:

  * PK12 - parallel 3d 12 sub-iteration thinning algorithm by Palagyi and Kuba
  * RC6  - parallel 3d 6 sub-iteration istmus-based thinning algorithms

References:
  Palagyi & Kuba, A Parallel 3D 12-Subiteration Thinning Algorithm, Graphical Models and Image Processing 61, 199-221 (1999)
  B. Raynal and M. Couprie, Istmus-Based 6-Directional Parallel Thinning Algorithms
"""

import PK12

from ClearMap.ImageProcessing.Topology.Topology3D import deleteBorder


def skeletonize3D(data, points = None, method = 'PK12i', steps = None, verbose = True, info = None, **kwargs):
    """Skeletonize 3d binary arrays
    
    Arguments:
      data (array): binary image
      method (string): 'PK12' or "RC6' algorithms, faster index versions 'PK12i' and 'RC6i'
      steps (int or None): number of maximal iteration steps
      info (None, or True): if True return additional info obtained during thinning (e.g. thickness)
      verbose (bool): if True print progrss info
      
    Returns:
      array: the skeleton
      array: the point coordinates of the skeleton nx3
    """
    
    if method == 'PK12':
      return PK12.skeletonize3D(data, points = points, steps = steps, verbose = verbose, **kwargs)
    elif method == 'PK12i':
      return PK12.skeletonize3DIndex(data, points = points, steps = steps, verbose = verbose, **kwargs)
    else:
      raise RuntimeError('No method %s' % method);



def test():
  import numpy as np;
  import ClearMap.ImageProcessing.Skeletonization.Skeletonize as s;
  data = np.load('test_bin.npy');
  data = s.deleteBorder(data);
  
  #data = np.zeros((5,5,5), dtype = bool);
  #data[1:4,1:4,1:4] = 1
  
  reload(s);
  skel, pts, death = s.skeletonize3D(data.copy(), verbose = True, info = True);
  
  radius = np.zeros(skel.shape, dtype = int);
  radius[pts[:,0], pts[:,1], pts[:,2]] = death;
  
  #plot the result
  import os;
  cmpath = '/home/ckirst/Science/Projects/WholeBrainClearing/Vasculature/Analysis/ClearVessel'
  skpath = '/home/ckirst/Science/Projects/WholeBrainClearing/Vasculature/Analysis/ClearVessel/ClearMap/ImageProcessing/Skeletonization'
  os.chdir(cmpath);
  import ClearMap.GUI.DataViewer as dv;
  dv.DataViewer([data, skel])
  os.chdir(skpath)
  
  
  from mayavi import mlab
  mlab.figure(bgcolor=(1,1,1));
  mlab.contour3d(np.array(data, dtype = int), color = (0.5, 0.5, 0.5), opacity = 0.2);
  mlab.contour3d(np.array(skel, dtype = int), color = (1,0,0), opacity = 0.2);
  
  import os;
  cmpath = '/home/ckirst/Science/Projects/WholeBrainClearing/Vasculature/Analysis/ClearVessel'
  skpath = '/home/ckirst/Science/Projects/WholeBrainClearing/Vasculature/Analysis/ClearVessel/ClearMap/ImageProcessing/Skeletonization'
  os.chdir(cmpath);
  import ClearMap.GUI.DataViewer as dv;
  dv.DataViewer([data, skel, death])
  os.chdir(skpath)
  
  from mayavi import mlab
  mlab.figure(bgcolor=(1,1,1));
  mlab.contour3d(np.array(data, dtype = int), color = (0.5, 0.5, 0.5), opacity = 0.2);
  mlab.contour3d(np.array(skel, dtype = int), color = (1,0,0), opacity = 0.9);
  mlab.contour3d(np.array(death, dtype = int), color = (0,1,0), opacity = 0.5);
  
  
  from mayavi import mlab
  mlab.figure(bgcolor=(1,1,1));
  mlab.contour3d(np.array(data, dtype = int), color = (0.5, 0.5, 0.5), opacity = 0.2);
  mlab.contour3d(np.array(radius, dtype = int), colormap = 'jet', opacity = 1);
  
  import os;
  cmpath = '/home/ckirst/Science/Projects/WholeBrainClearing/Vasculature/Analysis/ClearVessel'
  skpath = '/home/ckirst/Science/Projects/WholeBrainClearing/Vasculature/Analysis/ClearVessel/ClearMap/ImageProcessing/Skeletonization'
  os.chdir(cmpath);
  import ClearMap.GUI.DataViewer as dv;
  dv.DataViewer([data, skel, radius])
  os.chdir(skpath)
  