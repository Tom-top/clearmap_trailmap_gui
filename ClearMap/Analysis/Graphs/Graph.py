#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Graph analysis module

Wrapper interfacing graph_tool to ClearMap
"""

import copy 
import numpy as np

import graph_tool as gt

import ClearMap.IO as io
import ClearMap.Analysis.Label as lbl


###############################################################################
### Type Conversions
###############################################################################


def dtypeToGtType(dtype):
    """convert a data type to a graph_tool data type"""
    dtype = np.dtype(dtype);
    name = dtype.name;
    
    alias = { 'float64' : 'double', 
              'float32' : 'double',
              'int64'   : 'int',
              'int32'   : 'int'};
    if name in alias:
        name = alias[name];

    return gt._type_alias(name);


def dimToGtDim(data, dtype):
  if data.ndim == 2:
    dtype = "vector<%s>" % dtype;
  elif data.ndim > 2:
    raise ValueError('data for vertex properites can only be 1 or 2d');
  
  return dtype;

def dimToGtGraphDim(data, dtype):
  if len(data) > 0:
    dtype = "vector<%s>" % dtype;  
  return dtype;


def dataToGtType(data):
    """determines the graph_tool data type from the data"""
    if not isinstance(data, np.ndarray):
        data = np.array(data);
    dtype = data.dtype;
    dtype = dtypeToGtType(dtype);
    dtype = dimToGtDim(data, dtype);
    
    return dtype;


def dataToGtGraphType(data):
    """determines the graph_tool data type from the data"""
    if not isinstance(data, np.ndarray):
        data = np.array(data);
    dtype = data.dtype;
    dtype = dtypeToGtType(dtype);
    dtype = dimToGtGraphDim(data, dtype);
    
    return dtype;


###############################################################################
### Graph Class
###############################################################################

class Graph(object):
  """Graph class to handle graph construction and analysis
  
  Note: This is a wrapper class for a graph_tool class for use with ClearMap
  """
  
  def __init__(self, graph = None, name = None):
    self.name = name;
    if graph is None:
        graph = gt.Graph(directed = False);
    self.graph = graph;
    
  def copy(self):
    return Graph(name = copy.copy(self.name), graph = self.graph.copy())
  
  def nVertices(self):
    return self.graph.num_vertices();
    
  def nEdges(self):
    return self.graph.num_edges();
  
  def addVertices(self, n):
    self.graph.add_vertex(n = n);
      
  def addEdges(self, edges):
    self.graph.add_edge_list(edges);
  
  def edges(self):
    return  self.graph.get_edges()[:,:2];
  
  def vertices(self):
    return self.graph.vertices();
  
  def degrees(self):
    return self.graph.get_out_degrees(self.graph.get_vertices());
  
  def save(self, filename):
    self.graph.save(filename);
    
  def load(self, filename):
    self.graph = gt.load_graph(filename);
  
  def subGraph(self, vertexFilter = None, edgeFilter = None, view = False):
    gv = gt.GraphView(self.graph, vfilt = vertexFilter, efilt = edgeFilter);
    if view:
      return gv;
    else:    
      g = gt.Graph(gv, prune = True);
      return Graph(graph = g);
      
    
  
  def addVertexProperty(self, name, data, dtype = None):
    if dtype is None:
       dtype = dataToGtType(data);
    else:
      dtype = dtypeToGtType(dtype);
      dtype = dimToGtDim(data, dtype);

    p = self.graph.new_vertex_property(dtype);
    if data.ndim == 2:
      p.set_2d_array(data.T);
    else:
      p.a = data;
      
    self.graph.vertex_properties[name] = p;
  
  def setVertexProperty(self, name, data, add = True):
    if name not in self.graph.vertex_properties.keys():
      if add:
        self.addVertexProperty(name, data);
      else:
        raise ValueError('graph has not property %s' % name);
    else: 
      p = self.graph.vertex_properties[name];
      if data.ndim == 2:
        p.set_2d_array(data.T);
      else:
        p.a = data;
  
  def vertexProperty(self, name):
    p = self.graph.vertex_properties[name];
    if p.a is not None:
      return np.array(p.a);
    else:
      ndim = len(p[0]);
      return p.get_2d_array(range(ndim)).T;
    
  def vertexProperties(self):
    return self.graph.vertex_properties;  
  
  
  def addEdgeProperty(self, name, data, dtype = None):
    if dtype is None:
       dtype = dataToGtType(data);
    else:
      dtype = dtypeToGtType(dtype);
      if data.ndim == 2:
        dtype = "vector<%s>" % dtype;
      elif data.ndim > 2:
        raise ValueError('data for vertex properites can only be 1 or 2d');    
    
    p = self.graph.new_edge_property(dtype);
    if data.ndim == 2:
      p.set_2d_array(data.T);
    else:
      p.a = data;
    
    self.graph.edge_properties[name] = p;

  def setEdgeProperty(self, name, data, add = True):
    if name not in self.graph.edge_properties.keys():
      if add:
        self.addEdgeProperty(name, data);
      else:
        raise ValueError('graph has not property %s' % name);
    else: 
      p = self.graph.edge_properties[name];
      if data.ndim == 2:
        p.set_2d_array(data.T);
      else:
        p.a = data;
  
  def edgeProperty(self, name):
    p = self.graph.edge_properties[name];
    if p.a is not None:
      return np.array(p.a);
    else:
      ndim = len(p[0]);
      return p.get_2d_array(range(ndim)).T;
    
  def edgeProperties(self):
    return self.graph.edge_properties; 
  
  
  def addGraphProperty(self, name, data, dtype = None):
    if dtype is None:
       dtype = dataToGtGraphType(data);
    else:
      dtype = dataToGtGraphType(dtype);
      if len(data) > 0:
        dtype = "vector<%s>" % dtype;
      
    p = self.graph.new_graph_property(dtype);
    p.set_value(data);
    
    self.graph.graph_properties[name] = p; 
  
  def setGraphProperty(self, name, data, add = True):
    if name not in self.graph.graph_properties.keys():
      if add:
        self.addGraphProperty(name, data);
      else:
        raise ValueError('graph has not property %s' % name);
    else: 
      p = self.graph.graph_properties[name];
      p.set_value(data);
        
  def graphProperty(self, name):
    return self.graph.graph_properties[name];
    
  def graphProperties(self):
    return self.graph.graph_properties; 
  
  def removeSelfLoops(self):
    gt.stats.remove_self_loops(self.graph);
  
  
  ### Specialization to transport networks with vertices in space and edge radii
  def shape(self):
    """The shape of the underlying data space in which the graph is embedded"""
    return self.graphProperty('shape');
    
  def setShape(self, shape):
    return self.setGraphProperty('shape', shape);
  
  
  
  def vertexCoordinates(self, axis = None):
    if axis in [0,1,2]:
      axis = np.array(['x', 'y', 'z'])[axis];
    if axis in ['x', 'y', 'z']:
      return np.array(self.graph.vertex_properties[axis].a);
    
    return np.c_[self.graph.vertex_properties['x'].a, self.graph.vertex_properties['y'].a, self.graph.vertex_properties['z'].a];
  
  def setVertexCoordinates(self, coordinates, axis = None):
    if axis in [0,1,2]:
      axis = np.array(['x', 'y', 'z'])[axis];
    if axis in ['x', 'y', 'z']:
      self.setVertexProperty(axis,coordinates);
    else:
      x,y,z = coordinates.T;
      for a,n in zip([x,y,z], ['x', 'y', 'z']):
        self.setVertexProperty(n,a);
  
  def vertexRadius(self):
    return self.vertexProperty('r');
        
  def setVertexRadius(self, radii):
    self.setVertexProperty('r', radii);
  
  
  def edgeRadius(self):
    return self.edgeProperty('r');
  
  def setEdgeRadius(self, radii):
    self.setEdgeProperty('r', radii);
  
  def setEdgeRadiusFromVertexRadius(self):
    r = self.vertexRadius();
    i,j = self.edges().T;
    er = 0.5 * (r[i] + r[j]);
    self.setEdgeRadius(er);
  
  def edgeCoordinates(self, alpha = 0.5):
    vxyz = self.vertexCoordinates();
    i,j = self.edges();
    return alpha * vxyz[i] + (1-alpha) * vxyz[j];
  
  def edgeOrientations(self, normalize = False):
    vxyz = self.vertexCoordinates();
    i,j = self.edges();
    o = vxyz[i] - vxyz[j];
    if normalize:
      o = (o.T / np.linalg.norm(o, axis = 1)).T;
    return o;
  
  
  ### specialization to include detailed edge structures for reduced graphs
  def setEdgeGeometry(self, coordinates, radii, start, end):
    """Adds detailed geometrical information to each edge
    
    Arguments:
    ----------
      coordinates : nx3 array
        coordinates of all the individual edge points
      radii : array
        radii for tall the individual edge points
      start : array
        start index in data arrays above for a specific reduced branch edge
      end :
        end index in the data arrays above for a sepcific reduced branch edge
        
     Note
     ----
       The start end end arrays are treated as edge propeties, the detailed data as graph properties
    """
    
    self.setEdgeProperty('geometry_start', start);
    self.setEdgeProperty('geometry_end', end);
    
    self.setGraphProperty('edge_geometry_x', coordinates[:,0]);
    self.setGraphProperty('edge_geometry_y', coordinates[:,1]);
    self.setGraphProperty('edge_geometry_z', coordinates[:,2]);
    self.setGraphProperty('edge_geometry_r', radii);

  
  def edgeGeometry(self, edge = None):
    if edge is None:
      xyz = np.c_[self.graphProperty('edge_geometry_x'), self.graphProperty('edge_geometry_y'), self.graphProperty('edge_geometry_z')];
      return xyz, self.graphProperty('edge_geometry_r'), self.edgeProperty('geometry_start'), self.edgeProperty('geometry_end')
    else:
      s = self.edgeProperty('geometry_start')[edge];
      e = self.edgeProperty('geometry_end')[edge];
      xyz = np.c_[self.graphProperty('edge_geometry_x')[s:e], self.graphProperty('edge_geometry_y')[s:e], self.graphProperty('edge_geometry_z')[s:e]]
      r = self.graphProperty('edge_geometry_r')[s:e];
      return xyz,r
    
  def edgeGeometryCoordinates(self, edge = None):
    if edge is None:
      xyz = np.c_[self.graphProperty('edge_geometry_x'), self.graphProperty('edge_geometry_y'), self.graphProperty('edge_geometry_z')];
      return xyz
    else:
      s = self.edgeProperty('geometry_start')[edge];
      e = self.edgeProperty('geometry_end')[edge];
      xyz = np.c_[self.graphProperty('edge_geometry_x')[s:e], self.graphProperty('edge_geometry_y')[s:e], self.graphProperty('edge_geometry_z')[s:e]]
      return xyz

  
  ### specialization to labeled graphs aligned to an atlas
  def addLabel(self, annotation = lbl.AnnotationFile, key = 'id', value = 'order'):

    # label points
    aba = np.array(io.readData(annotation), dtype = int);
    
    # get vertex coordinates
    x,y,z = self.vertexCoordinates().T;
  
    ids = np.ones(len(x), dtype = bool);
    for a,s in zip([x,y,z], aba.shape):
      ids = np.logical_and(ids, a >= 0);
      ids = np.logical_and(ids, a < s);
  
    # label points
    g_ids = np.zeros(len(x), dtype = int);
    g_ids[ids] = aba[x[ids],y[ids],z[ids]];
  
    if value is not None:
      id_to_order = lbl.getMap(key = key, value = value)
      g_order = id_to_order[g_ids];
    else:
      value = key;
    
    self.addVertexProperty(value, g_order);
  
  
  def __str__(self):
    if self.name is not None:
      return 'Graph %r: %d vertices, %d edges' % (self.name, self.nVertices(), self.nEdges());
    else:
      return 'Graph: %d vertices, %d edges' % (self.nVertices(), self.nEdges());

    
  def __repr__(self):
    return self.__str__();


def load(filename):
  g = gt.load_graph(filename);
  return Graph(graph = g);

def save(filename, graph):
  g.save(filename);


if __name__ == '__main__':
  import ClearMap.Analysis.Graphs.Graph as gr
  reload(gr)
  
  g = gr.Graph('test');
  
  g.addVertices(10);
  
  el = [[1,3],[2,5],[6,7]];
  g.addEdges(el)
  
  print g
  coords = np.random.rand(10,3);
  g.setVertexCoordinates(coords)




  
  
  
  