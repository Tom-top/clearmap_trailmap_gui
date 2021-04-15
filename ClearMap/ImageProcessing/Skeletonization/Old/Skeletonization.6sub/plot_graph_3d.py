# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 19:11:24 2016

@author: ckirst
"""

import networkx as nx
import numpy as np
from mayavi import mlab

import os

def plot_graph_3d(graph, colormap='jet', radii = None,
                  line_width = 2, opacity=.9):

    # get graph positions
    g2 = nx.convert_node_labels_to_integers(graph, label_attribute = 'xyz');
    xyz = np.array([x['xyz'] for x in g2.node.values()]);

    # scalar colors
    if radii is not None:
      scalars = [radii[tuple(x)] for x in xyz];
    else:
      #scalars = np.arange(5, xyz.shape[0]+5);
      scalars = np.ones(xyz.shape[0]);
    
    #pts = mlab.points3d(xyz[:,0], xyz[:,1], xyz[:,2],
    #                    scalars,
    #                    scale_factor=node_size,
    #                    scale_mode='none',
    #                    colormap=graph_colormap,
    #                    resolution=20)

    pts = mlab.pipeline.scalar_scatter(xyz[:,0], xyz[:,1], xyz[:,2], scalars)


    pts.mlab_source.dataset.lines = np.array(g2.edges())
    pts.update()    
    
    #tube = mlab.pipeline.tube(pts, tube_radius=edge_size)
    #lab.pipeline.surface(tube, color=edge_color)
    
    lines = mlab.pipeline.stripper(pts);
    mlab.pipeline.surface(lines, colormap = colormap, line_width = line_width, opacity = opacity)
    
    mlab.colorbar(orientation = 'vertical', title='Radius [pix]');    
    
    return lines