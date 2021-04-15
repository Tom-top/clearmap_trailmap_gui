#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 16:34:16 2017

"""


import numpy as np

data_bin = np.load('bin.npy')

def zero_border(data):
  data[[0,-1],:,:] = False;
  data[:,[0,-1],:] = False;
  data[:,:,[0,-1]] = False;
  return data;

data_bin = zero_border(data_bin);

#%%
import ClearMap.ImageProcessing.Skeletonization.skeletonize_3d_list as s3d;
reload(s3d)
data_test = data_bin.copy();
skel = s3d.skeletonize(data_test, pad_with_zeros=False, verbose = True, fast = True);


#%% compare with working method
from ClearMap.Utils.Timer import Timer;
import ClearMap.ImageProcessing.Skeletonization.skeletonize_3d as s3d0;
reload(s3d0);
data_test = data_bin.copy();
timer = Timer();
skel2= s3d0.skeletonize(data_test, pad_with_zeros=False, verbose = True, fast = True);
timer.printElapsedTime();

#%%

print 'tot: %d' % data_bin.sum();
print 'new: %d' % skel.sum();
print 'old: %d' % skel2.sum()

skel_diff = np.array(skel2, dtype = int) - np.array(skel, dtype = int);



#%% Plot
import numpy as np

import ClearMap.ImageProcessing.Skeletonization.skeleton_graph as s2g;
import ClearMap.ImageProcessing.Skeletonization.plot_graph_3d as spg  
reload(s2g); reload(spg);

skel = np.load('skel_test.npy')
skel = np.load('skel.npy')
data_bin = np.load('bin.npy')

g = s2g.skeleton_to_nx_graph(s2g.ensure_zero_border(skel));

from mayavi import mlab;
mlab.figure(bgcolor = (0,0,0))
mlab.contour3d(np.array(data_bin, dtype = int) * 100, color = (0.5, 0.5, 0.5), opacity = 0.1)

spg.plot_graph_3d(g, colormap='hsv')



#%%



reload(s3d);
od = s3d.opposing_direction_list()


reload(s3d);
b = s3d.neighbourhood_to_block(nh[0])

reload(s3d)
timer = Timer();
x,y,z,nh = s3d.get_neighbourhood_list_2(data_bin);
timer.printElapsedTime();

# test data

data_bin = np.array(np.random.rand(5,5,5) > 0.5, dtype = bool);
data_bin = zero_border(data_bin);
x,y,z,nhl = s3d.get_neighbourhood_list_2(data_bin);

nhl2 = nhl.copy();
nhl3 = s3d.delete_in_neighbours([0,1], nhl2)

timer = Timer();
skel = s3d.skeletonize(data_bin, pad_with_zeros=True, verbose = True, fast = True);
timer.printElapsedTime();


