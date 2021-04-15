# -*- coding: utf-8 -*-
"""
Created on Tue Oct 25 12:36:33 2016

@author: ckirst
"""


import os;
import numpy as np;

dir_clearmap = '/data/Science/Projects/BrainActivityMap/Analysis/ClearMapUnstable';
dir_skeleton = '/data/Science/Projects/BrainVasculature/Analysis/Skeletonization';


dir_base = '/data/Science/Projects/BrainVasculature/Experiment/'
file_data   = os.path.join(dir_base, 'stitched/Z\d{4}.tif');
file_npy    = os.path.join(dir_base, 'stitched.npy');
file_skel   = os.path.join(dir_skeleton, 'stitched_skeleton.npy');
file_bin    = os.path.join(dir_base, 'stitched_thresholded.npy');
file_closed = os.path.join(dir_base, 'stitched_closed.npy');
file_graph  = os.path.join(dir_skeleton, 'vasculature_graph.gt');

file_bin_test = '/data/Science/Projects/BrainVasculature/Experiment/stitched_thresholded_test.tif';


### Prepare Data


os.chdir(dir_clearmap)
import ClearMap.IO as io

data_shape = io.dataSize(file_data)

### Load Data and write as nrrd file

data = np.load(file_skel)

#write data to an nrrd file

file_image = os.path.join(dir_skeleton, 'skeleton.nrrd');
data = np.array(data, dtype = 'uint8');
io.writeData(file_image, data);


### Load Graph

os.chdir(dir_skeleton)
import numpy as np
import skeleton_graph as sg

graph = sg.gt.load_graph(file_graph)


# get graph positions
x = np.array(graph.vertex_properties['x'].get_array(), dtype = 'int32');
y = np.array(graph.vertex_properties['y'].get_array(), dtype = 'int32');
z = np.array(graph.vertex_properties['z'].get_array(), dtype = 'int32');


xlim = [1300, 1500];
ylim = [0,all];
zlim = [0,all];


def sub_idx(x,y,z,xlim, ylim,zlim):
   idx = np.ones(x.shape[0], dtype = bool);
   for d,l in zip([x,y,z],[xlim,ylim,zlim]):
     if l[0] is not all:
       idx = np.logical_and(idx, l[0] <= d);
     if l[1] is not all:
       idx = np.logical_and(idx, d < l[1]);
   return idx;
   
idx = sub_idx(x,y,z,xlim,ylim,zlim);

xs = x[idx];
ys = y[idx];
zs = z[idx];

scalars = np.ones(xs.shape[0], dtype = 'float32');   
  
ids = np.where(idx)[0];
   

from graph_tool import GraphView, generation

# select some vertices
#vfilt = graph.new_vertex_property('bool');
#vfilt.fa = idx;
#sub = GraphView(g, vfilt)


import graph_tool.spectral as gs

m = gs.adjacency(g)
m_sub = m[ids,:][:,ids]

i,j = m_sub.nonzero();

edgelist = np.vstack([i,j]).T;

from mayavi import mlab

pts = mlab.pipeline.scalar_scatter(xs, ys, zs, scalars)

pts.mlab_source.dataset.lines = edgelist;
pts.update()    
  
lines = mlab.pipeline.stripper(pts);
colormap='jet';
line_width = 2;
opacity=.9;
mlab.pipeline.surface(lines, colormap = colormap, line_width = line_width, opacity = opacity)
  
mlab.axes();


