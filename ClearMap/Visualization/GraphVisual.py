#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Graph Visual to display a graph in 3d using tube rendering for the edges

Notes: The graph is assumed to be a reduced graph with each edge given as a list of coordinates and radii
"""

from __future__ import division

import numpy as np

from vispy.visuals.mesh import MeshVisual

from vispy.util.transforms import rotate
from vispy.color import ColorArray


class GraphVisual(MeshVisual):
    """Displays a graph in 3d using tube rendering for the edges
    """
    
    def __init__(self, graph, n_tube_points=8, shading = 'smooth', colors=None, mode='triangles'):
      pass
      
      
      
import ClearMap.ParallelProcessing.SharedMemoryManager as smm
    
def edgeGeometryToMesh(points, radii, start, end):
  
  points_hdl   = smm.insert(points);
  radii_hdl = smm.insert(radii);
  start_hdl = smm.insert(start);
  end_hdl   = smm.insert(end);
  
  func = partial(_parallel_mesh, points_hdl = points_hdl, radii_hdl = radii_hdl, start_hdl  = start_hdl, end_hdl = end_hdl);
  argdata = np.arange(len(start));
  
  # process in parallel
  pool = Pool(processes = processes);    
  results = pool.map(func, argdata);
  pool.close();
  pool.join();
  
  smm.free(points_hdl);
  smm.free(radii_hdl);
  smm.free(start_hdl);
  smm.free(end_hdl);  


def _parallel_mesh(i, points_hdl, radii_hdl, start_hdl, end_hdl):

  points = smm.get(points_hdl);
  radii  = smm.get(radii_hdl);
  start  = smm.get(start_hdl)[i];
  end    = smm.get(end_hdl)[i];
  
  points = points[start:end];
  radii  = radii[start:end];
  
  return _mesh(points, radii)


def _frenet_frames(points):
    """Calculates and returns the tangents, normals and binormals for a chain of points"""
    npoints = len(points);
    
    epsilon = 0.0001

    # compute tangent vectors for each segment
    tangents = np.zeros((npoints, 3));
    tangents[:-1,:] = np.diff(points, axis = 0);
    tangents[-1] = points[-1] - points[-2]
    tangents = (tangents.T / np.linalg.norm(tangents, axis = 1)).T;

    # get initial normal and binormal
    t = np.abs(tangents[0])

    smallest = np.argmin(t)
    normal = np.zeros(3)
    normal[smallest] = 1.

    vec = np.cross(tangents[0], normal)

    normals = np.zeros((npoints, 3));
    normals[0] = np.cross(tangents[0], vec)

    # compute changea along trajectory 
    theta = np.arccos(np.clip(np.sum(tangents[:1] * tangents[1:], axis = 1),-1,1));
    vec = np.cross(tangents[:1], tangents[1:]);
    nrm = np.linalg.norm(vec, axis = 1);

    # compute normal and binormal vectors along the path
    for i in range(npoints-1):
      normals[i+1] = normals[i]
      
      if nrm[i] > epsilon:
        v = vec[i] / nrm[i];
        normals[i+1] = rotate(-np.degrees(theta[i]), v)[:3, :3].dot(normals[i+1])
    
    binormals = np.cross(tangents, normals)

    return tangents, normals, binormals
  

def _mesh(points, radii, ntube = 8):
  npoints = len(points);
  tangents, normals, binormals = _frenet_frames_2(points)

  #circular tube
  v = np.arange(ntube, dtype = float) / ntube * 2 * np.pi;
  c = np.cos(v);
  s = np.sin(v);
  
  r = radii[:, np.newaxis];
  n = normals * r;
  b = binormals * r;
  
  #grid shape (npoints, ntube, 3)
  grid = points[:, np.newaxis, :] + c[np.newaxis, :, np.newaxis] * n[:, np.newaxis, :] + s[np.newaxis, :, np.newaxis] * b[:, np.newaxis, :];

  # construct the mesh
  nsegments = npoints - 1;
  jp = np.ones(nsegments*ntube, dtype = int);
  jp[ntube * np.arange(nsegments, dtype = int) - 1] -= ntube;
  
  i1 = np.arange(nsegments*ntube);
  i2 = i1 + ntube;
  i3 = i2 + jp;
  i4 = i1 + jp;
  
  indices = np.vstack([np.array([i1,i2,i4]).T, np.array([i2, i3, i4]).T]);
  
  return grid, indices

import vispy
import vispy.scene

def plotMesh(grid, indices, view = None):
   # build visuals
  mv = vispy.scene.visuals.create_visual_node(MeshVisual)

  if view is None:
      # build canvas
      canvas = vispy.scene.SceneCanvas(keys='interactive', title='plot3d', show=True)

      # Add a ViewBox to let the user zoom/rotate
      view = canvas.central_widget.add_view()
      view.camera = 'turntable'
      view.camera.fov = 100
      view.camera.distance = 0
      view.camera.elevation = 0
      view.camera.azimuth = 0

  return mv(grid, indices, parent=view.scene)





def _frenet_frames_2(points, closed = False):
    '''Calculates and returns the tangents, normals and binormals for
    the tube.'''
    tangents = np.zeros((len(points), 3))
    normals = np.zeros((len(points), 3))

    epsilon = 0.0001

    # Compute tangent vectors for each segment
    tangents = np.roll(points, -1, axis=0) - np.roll(points, 1, axis=0)
    if not closed:
        tangents[0] = points[1] - points[0]
        tangents[-1] = points[-1] - points[-2]
    mags = np.sqrt(np.sum(tangents * tangents, axis=1))

    tangents /= mags[:, np.newaxis]

    # Get initial normal and binormal
    t = np.abs(tangents[0])

    smallest = np.argmin(t)
    normal = np.zeros(3)
    normal[smallest] = 1.

    vec = np.cross(tangents[0], normal)

    normals[0] = np.cross(tangents[0], vec)

    # Compute normal and binormal vectors along the path
    for i in range(1, len(points)):
        normals[i] = normals[i - 1]

        vec = np.cross(tangents[i - 1], tangents[i])
        if np.linalg.norm(vec) > epsilon:
            vec /= np.linalg.norm(vec)
            theta = np.arccos(np.clip(tangents[i - 1].dot(tangents[i]), -1, 1))
            normals[i] = rotate(-np.degrees(theta),
                                vec)[:3, :3].dot(normals[i])

    if closed:
        theta = np.arccos(np.clip(normals[0].dot(normals[-1]), -1, 1))
        theta /= len(points) - 1

        if tangents[0].dot(np.cross(normals[0], normals[-1])) > 0:
            theta *= -1.

        for i in range(1, len(points)):
            normals[i] = rotate(-np.degrees(theta * i),
                                tangents[i])[:3, :3].dot(normals[i])

    binormals = np.cross(tangents, normals)

    return tangents, normals, binormals




if  __name__ == "__main__":
  import numpy as np  
  import ClearMap.Visualization.GraphVisual as gv
  reload(gv)
  
  ne = 1;
  n = 10;
  nn = n * ne;
  t = np.arange(nn);
  points = np.array([t, np.sin(t), np.cos(t)]).T;
  radii = np.sin(t) + 2;
  radii = hnp.ones(len(t));
  start = np.array([0]);
  end = np.array([nn]);
  
  reload(gv)
  grid, indices = gv._mesh(points, radii)
  grid = np.reshape(grid, (-1,3));
  
  from vispy.visuals.mesh import MeshVisual
  
  gv.plotMesh(grid, indices)