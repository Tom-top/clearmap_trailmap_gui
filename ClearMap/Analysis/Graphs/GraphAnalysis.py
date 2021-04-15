# -*- coding: utf-8 -*-
"""
Graph Analysis Module

Routines to convert skeletons to graphs and analyse them
"""
__license__ = 'MIT License <http://www.opensource.org/licenses/mit-license.php>'
__author__ = 'Christoph Kirst <ckirst@rockefeller.edu>'
__docformat__ = 'rest'

import numpy as np
import networkx as nx
from mayavi import mlab

import ClearMap.ImageProcessing.Skeletonization.Topology3D as top3d



#TODO: join with SkeletonGraph routines -> homogenize to Graph class

### Skeleton To Graph Converson
def skeletonToAdjacencyMatrix(skeleton, verbose = False):
  """Converts a binary skeleton image to a adjacency matrix
  
  Arguments:
    skeleton (array): 2d/3d binary skeleton image
    
  Returns:
    dict: dict of adjacency information with entries node_id : [neighbours]
  """
  
  x,y,z,nh = top3d.neighbourhoodList(skeleton);
  ids = np.transpose([x,y,z]);
 
  adj = {}; 
  if len(ids) == 1:
    adj[tuple(ids[0])] = [];
  elif len(ids) > 1:
    for i,pos in enumerate(ids):
      if verbose and i % 1000 == 0:
        print('adjacency %d / %d' % (i, len(ids)));      
      posnh = np.where(nh[i]);
      adj[tuple(pos)] = [tuple(p + pos -1) for p in np.transpose(posnh)]
  
  return adj;
    
    
def skeletonToNxGraph(skeleton, convertIndices = True, verbose = False):
  """Converts a binary skeleton image to a networkx graph
  
  Arguments:
    skeleton (array): 2d/3d binary skeleton image
    values (array): 2d/3d array of same size as skeleton with values for the skeleton points (e.g. radius data)
    convertIndices (bool): if True convert the nodes to have integer labels adn xyz attributes
    verbose (bool): if True print info during proecssing
    
  Returns:
    graph: networkx graph object
  """
  
  x,y,z,nh = top3d.neighbourhoodList(skeleton);
  ids = np.transpose([x,y,z]);
  nh = np.reshape(nh,(nh.shape[0], 3,3,3), order = 'F')
  nh[:, 1, 1, 1] = -1;
  if verbose:
    print('Graph: information extracted from skeleton...')
 
  if len(ids) == 0:
     return nx.Graph();
  elif len(ids) == 1:
    adj = {};
    adj[tuple(ids[0])] = [];
    g = nx.from_dict_of_lists(adj);
  else:
    g = nx.Graph();
    for i,pos in enumerate(ids):
      if verbose and i %1000 == 0:
          print('Graph: %d/%d nodes constructed...' % (i, len(ids)));
      p = tuple(pos);
      g.add_node(p);
      posnh = np.where(nh[i]>=0);
      for pp in np.transpose(posnh):
          g.add_edge(p, tuple(pp+pos-1));
  
  if convertIndices is True:
    g = nx.convert_node_labels_to_integers(g, label_attribute = 'xyz');
  
  return g;



### Clean Graph from cliques arising at branch points in skeleton

def cleanGraph(graph, convertIndices = True, verbose = True):
  """Remove all 3 cliques to get pure branch structure"""
  g = graph.copy();
  deg = np.array([g.degree(i) for i in g]);
  ids = np.where(deg >= 3)[0];
  if verbose:
    print('Graph cleaning: found %d branch points' % len(ids));
  if len(ids) ==0 :
    return g;

  gsub = nx.subgraph(g, ids);
  cc = list(nx.connected_components(gsub));
  l = np.array([len(c) for c in cc]);
  cids = np.where(l > 1)[0];
  cadj = [cc[i] for i in cids];
  if verbose:
    print('Graph cleaning: found %d meta branch points' % len(cids));
  
  # remove meta cliques
  nnodes = g.number_of_nodes();
  for i,ci in enumerate(cadj):
    if verbose and i % 1000 == 0:
      print('Graph cleaning: reducing %d / %d' % (i, len(cadj)));
      
    # get nodes coordinates
    xyz = np.array([g.node[n]['xyz'] for n in ci]);
    
    # get neighbours 
    nbs = set(np.hstack([g.adj[n].keys() for n in ci]));
    nbs = nbs - ci;
    
    # remove all meta clique nodes
    g.remove_nodes_from(ci);
    
    # generate new center node
    newxyz = tuple(np.array(np.mean(xyz, axis = 0), dtype = int));
    g.add_node(nnodes, xyz = newxyz);   
    
    # reconnect
    g.add_edges_from([(c, nnodes) for c in nbs]);

    nnodes += 1;
  
  if convertIndices is True:
    g = nx.convert_node_labels_to_integers(g);  
  
  return g;



def relabelGraph(graph, startIndex = 0, saveOldIndices = None):
  """Relabels the nodes to consecutive integers"""
  g = nx.convert_node_labels_to_integers(graph, first_label = startIndex, label_attribute = saveOldIndices);
  return g;

### Reduce graph sturcture

def reduceGraph(graph, verbose = True, convertIndices = True, saveOldIndices = None):
  """Reduce graph by replacing all vertices with degree two"""
  
  #TODO: include all positions, path length and width estimate on each link
  g = graph.copy();
  flag = 0;
  step = 0;
  while flag == 0:
    if verbose:
      print('Graph reduction: iteration: %d' % step);
      step += 1;
    
    flag = 1
    for i in g.nodes():
      #if g.degree(i) < 2:
      #  g.remove_node(i)
      #  flag = 0
      if g.degree(i) == 2:
        u,v = g[i].keys();
        g.add_edge(u,v)
        g.remove_node(i)
        flag = 0;
  
  if convertIndices is True:
    g = nx.convert_node_labels_to_integers(g, label_attribute = saveOldIndices);  
  
  return g;


def subGraph(graph, region):
  """Extracts the graph specified by the region
  
  Arguments:
    graph: network nx graph
    region: list of tuples or slices or all for each direction, 
  
  Returns:
    nxgraph: reduced graph
  """
  
  xyz = nodeCoordinates(graph);
  
  valid = np.ones(len(xyz), dtype = bool);
  for d in range(3):
    r = region[d];
    if isinstance(r, slice):
      valid = np.logical_and(valid, r.start <= xyz[:,d]);
      valid = np.logical_and(valid, xyz[:,d] < r.end);
    elif isinstance(r, list) or isinstance(r, tuple):
      valid = np.logical_and(valid, r[0] <= xyz[:,d]);
      valid = np.logical_and(valid, xyz[:,d] < r[1]);
    #else: every thing is valid
  
  ids = np.where(valid)[0];
  return nx.subgraph(graph, ids);


#def extractGraph()




### Reduce graph sturcture
#
#def reduceGraph2(graph, verbose = True, convertIndices = True, maximalLength = None):
#  """Reduce graph by replacing all vertices with degree two"""
#  
#  #TODO: include all positions, path length and width estimate on each link
#  g = graph.copy();
#  nn = g.number_of_nodes();
#  
#  deg = np.array([g.degree(i) for i in g]);
#  np.arange(nn);
#  bids = np.where(deg != 2)[0];
#  cids = np.where(deg)
#  if verbose:
#    print('Graph reduction: found %d branch points' % len(ids));
#  if len(ids) == 0 :
#    return g;
#  
#  flag = 0;
#  step = 0;
#  while flag == 0:
#    if verbose:
#      print('Graph reduction: iteration: %d' % step);
#    
#    flag = 1
#    for i in g:
#      #if g.degree(i) < 2:
#      #  g.remove_node(i)
#      #  flag = 0
#      if g.degree(i) == 2:
#        u,v = g[i].keys();
#        g.add_edge(u,v)
#        g.remove_node(i)
#        flag = 0;
#  
#  if convertIndices is True:
#    g = nx.convert_node_labels_to_integers(g);  
#  
#  return g;


def nodeCoordinates(graph):
  """Returns the cooordinates of the nodes of a graph"""
  return np.array([graph.node[n]['xyz'] for n in graph]);


### Transform a graph onto a reference

import ClearMap.Alignment.Elastix as elx;

def transformGraph(graph, transformParameterFile = None, transformDirectory = None, indices = True, resultDirectory = None, tmpFile = None, verbose = True):
    """Transform coordinates of vertices in a graph math:`x` via elastix estimated transformation to :math:`T(x)`

    Note the transformation is from the fixed image coorindates to the moving image coordiantes.
    
    Arguments:
        graph (networkx graph): graph to be transformed, coordinates are given as 'xyz' node attribute
        transformParameterFile (str or None): parameter file for the primary transformation, if None, the file is determined from the transformDirectory.
        transformDirectory (str or None): result directory of elastix alignment, if None the transformParameterFile has to be given.
        indices (bool): if True use points as pixel coordinates otherwise spatial coordinates.
        resultDirectory (str or None): elastic result directory
        tmpFile (str or None): file name for the elastix point file.
        
    Returns:
        graph: graph with transformed point coordinates
    """
    
    points = nodeCoordinates(graph);
    
    if verbose:
      print('Transforming %d node coordinates' % len(points));
    newpoints = elx.transformPoints(points, sink = None, 
                                    transformParameterFile = transformParameterFile, transformDirectory = transformDirectory, 
                                    indices = indices, resultDirectory = resultDirectory, tmpFile = tmpFile);
                                    
    if verbose:
      print('Creating new graph...');
    values = dict(zip(graph.nodes(), newpoints));
    nx.set_node_attributes(graph, 'xyz', values);

    return graph;





### Some Measures from the graph

def nodeDegrees(graph):
  return np.array([graph.degree(i) for i in graph]);

def branchNodes(graph, degree = all):
  """Returns list of branch points in the graph"""
  deg = nodeDegrees(graph);
  if degree is all:
    return np.where(deg >= 3)[0];
  else:
    return np.where(deg = degree)[0];




### Plotting

def plotNxGraph3d(graph, radii = None,  cmap = 'jet', line_width = 2, opacity = .9, zfactor = 1.0, colorbar = True):
  """Plot a 3d graph of the skeleton
  
  Arguments:
    radii: radii of the edges used in color code, if None uniform color
    
  Returns:
    mayavi scence
  """
  # get graph positions
  xyz = nx.get_node_attributes(graph, 'xyz');
  xyz = np.array([xyz[n] for n in graph.nodes()], dtype = 'int32');
  #g2 = nx.convert_node_labels_to_integers(graph);

  # scalar colors
  if radii is not None:
    scalars = np.array([radii[tuple(x)] for x in xyz], dtype = 'float32');
  else:
    #scalars = np.arange(5, xyz.shape[0]+5);
    scalars = np.ones(xyz.shape[0], dtype = 'float32');
  
  #pts = mlab.points3d(xyz[:,0], xyz[:,1], xyz[:,2],
  #                    scalars,
  #                    scale_factor=node_size,
  #                    scale_mode='none',
  #                    colormap=graph_colormap,
  #                    resolution=20)

  pts = mlab.pipeline.scalar_scatter(xyz[:,0], xyz[:,1], zfactor * xyz[:,2], scalars)
  
  pts.mlab_source.dataset.lines = np.array(graph.edges(), dtype = 'int32')
  pts.update()    
  
  #tube = mlab.pipeline.tube(pts, tube_radius=edge_size)
  #lab.pipeline.surface(tube, color=edge_color)
  
  lines = mlab.pipeline.stripper(pts);
  mlab.pipeline.surface(lines, colormap = cmap, line_width = line_width, opacity = opacity)
  
  if radii is not None and colorbar is True:
      mlab.colorbar(orientation = 'vertical', title='Radius [pix]');    
      
  mlab.axes()
  
  return lines




### Tests

if __name__ == "__main__":
  import numpy as np
  import matplotlib.pyplot as plt
  from mayavi import mlab
  import ClearMap.Analysis.GraphAnalysis as sg
  
  import ClearMap.ImageProcessing.Skeletonization.Skeletonize as s;
  data = np.load('./ClearMap/ImageProcessing/Skeletonization/test_bin.npy');
  data = s.deleteBorder(data);
  skel, pts = s.skeletonize3D(data.copy(), verbose = True, info = False);
  
  mlab.figure(bgcolor=(1,1,1));
  mlab.contour3d(np.array(data, dtype = int), color = (0.5, 0.5, 0.5), opacity = 0.2);
  mlab.contour3d(np.array(skel, dtype = int), color = (1,0,0), opacity = 0.2);  
  
  reload(sg)
  print(skel.shape, skel.sum())
  g = sg.skeletonToNxGraph(skel);
  
  mlab.figure(bgcolor = (1,1,1))
  sg.plotNxGraph3d(g, cmap = 'hot');
  mlab.contour3d(np.array(data, dtype = int), color = (0.5, 0.5, 0.5), opacity = 0.2);
  
  #gr = sg.reduceGraph(g);
  gc = sg.cleanGraph(g);
  
  mlab.figure(bgcolor = (1,1,1))
  sg.plotNxGraph3d(gc, cmap = 'hot');
  mlab.contour3d(np.array(data, dtype = int), color = (0.5, 0.5, 0.5), opacity = 0.2);

  gs = sg.reduceGraph(gc);
  mlab.figure(bgcolor = (1,1,1))
  sg.plotNxGraph3d(gs, cmap = 'hot');
  mlab.contour3d(np.array(data, dtype = int), color = (0.5, 0.5, 0.5), opacity = 0.2);









#import graph_tool as gt;
#
#
#def skeleton_to_gt_graph(skeleton, with_coordinates = True, verbose = True):
#  """Converts a binary skeleton image to a networkx graph
#  
#  Arguments:
#    skeleton (array): 2d/3d binary skeleton image
#    
#  Returns:
#    dict: dict of adjacency information with entries node_id : [neighbours]
#  """
#  dim = skeleton.ndim;
#  shape =skeleton.shape;
#  
#  coords, nh = skeleton_to_list(skeleton, with_neighborhoods = True);
#  nnodes = coords.shape[0];
#  
#  coordids = np.ravel_multi_index(coords.T, shape);  
#  coordids2ids = { k:i for i,k in enumerate(coordids) };
#    
#  # create graph
#  if verbose:
#    print 'creating graph...';
#  g = gt.Graph(directed = False);
#  g.add_vertex(nnodes);
#  if with_coordinates:
#      if verbose:
#        print 'creating coordinate properties...'
#      vp = g.new_vertex_property('int', coords[:,0]);
#      g.vertex_properties['x'] = vp;
#      vp = g.new_vertex_property('int', coords[:,1]);
#      g.vertex_properties['y'] = vp;
#      if dim > 2:
#          vp = g.new_vertex_property('int', coords[:,2]);
#          g.vertex_properties['z'] = vp;
#  
#  for i, pos in enumerate(coords):
#    if verbose and i % 1000 == 0:
#      print '%d/%d nodes constructed...' % (i, len(coords));
#    #print 'nh'
#    #print nh[i]
#    posnh = np.transpose(np.where(nh[i]));
#    #print 'pos'
#    #print pos
#    #print posnh
#    if len(posnh) > 0:
#      posnh = np.array(posnh + pos - 1);
#      #print posnh
#      #print posnh.shape
#      ids = np.ravel_multi_index(posnh.T, shape);
#      #print ids
#      for j in [coordids2ids[k] for k in ids]:
#        if i < j:
#          g.add_edge(i,j);
#  
#  return g

#
#def plot_gt_graph_3d(graph, radii = None,  colormap='jet', line_width = 2, opacity=.9):
#  """Plot a 3d graph of the skeleton
#  
#  Arguments:
#    radii: radii of the edges used in color code, if None uniform color
#    
#  Returns:
#    mayavi scence
#  """
#  # get graph positions
#  x = np.array(graph.vertex_properties['x'].get_array(), dtype = 'int32');
#  y = np.array(graph.vertex_properties['y'].get_array(), dtype = 'int32');
#  z = np.array(graph.vertex_properties['z'].get_array(), dtype = 'int32');
#
#  # scalar colors
#  if radii is not None:
#    #scalars = [radii[tuple(x)] for x in xyz];
#    scalars = np.array(radii, dtype = 'float32');
#  else:
#    #scalars = np.arange(5, xyz.shape[0]+5);
#    scalars = np.ones(x.shape[0], dtype = 'float32');
#  
#  #pts = mlab.points3d(xyz[:,0], xyz[:,1], xyz[:,2],
#  #                    scalars,
#  #                    scale_factor=node_size,
#  #                    scale_mode='none',
#  #                    colormap=graph_colormap,
#  #                    resolution=20)
#
#  pts = mlab.pipeline.scalar_scatter(x, y, z, scalars)
#  
#  edgelist = np.vstack([np.array([e.source(), e.target()], dtype = 'int32') for e in graph.edges()]);
#  pts.mlab_source.dataset.lines = edgelist;
#  pts.update()    
#  
#  #tube = mlab.pipeline.tube(pts, tube_radius=edge_size)
#  #lab.pipeline.surface(tube, color=edge_color)
#  
#  lines = mlab.pipeline.stripper(pts);
#  mlab.pipeline.surface(lines, colormap = colormap, line_width = line_width, opacity = opacity)
#  
#  if radii is not None:
#      mlab.colorbar(orientation = 'vertical', title='Radius [pix]');    
#  
#  mlab.axes();
#  
#  return lines


### Carls code
#  g = graph.copy();
#  flag = 0
#  while flag == 0:
#    flag = 1
#    for i in g.copy():
#      #if g.degree(i) < 2:
#      #    g.remove_node(i)
#      #    flag = 0
#      if g.degree(i) == 2:
#        u,v = g[i].keys();
#        g.add_edge(u,v)
#        g.remove_node(i)
#        flag = 0
#            
#flag = 0
#while flag == 0:
#    flag = 1
#    flag2 = 0
#    for i in g.copy():
#        if g.has_node(i):
#          if g.degree(i) > 2:
#                for j in g[i].copy():
#                    if (((i[0]-j[0])**2 + (i[1]-j[1])**2 + (i[2]-j[2])**2) < 25) and g.degree(j) > 2:
#                        g.add_node(((i[0]+j[0])/2.0,(i[1]+j[1])/2.0,(i[2]+j[2])/2.0))
#                        #dictOfNodesAndRadius[((i[0]+j[0])/2.0,(i[1]+j[1])/2.0,(i[2]+j[2])/2.0)] = (dictOfNodesAndRadius[i] + dictOfNodesAndRadius[j])/2.0
#                        if g.has_edge(i,j):
#                            g.remove_edge(i,j)
#                        if g.has_node(j):
#                            for k in g[j].copy():
#                                g.add_edge(((i[0]+j[0])/2.0,(i[1]+j[1])/2.0,(i[2]+j[2])/2.0),k)
#                        if g.has_node(i):
#                            for k in g[i].copy():
#                                g.add_edge(((i[0]+j[0])/2.0,(i[1]+j[1])/2.0,(i[2]+j[2])/2.0),k)
#                        if g.has_node(j):
#                            g.remove_node(j)
#                            #del dictOfNodesAndRadius[j]
#                        if g.has_node(i):
#                            g.remove_node(i)
#                            #del dictOfNodesAndRadius[i]
#                        flag = 0
#                        flag2 = 1
#                        break
#        if flag2 == 1:
#            break
#            
#flag = 0
#while flag == 0:
#    flag = 1
#    for i in g.copy():
#        if g.degree(i) < 2:
#            g.remove_node(i)
#            flag = 0
#        elif g.degree(i) == 2:
#            loopcounter = 1
#            for j in g[i]:
#                if loopcounter == 1:
#                    u = j
#                    loopcounter = 2
#                else:
#                    v = j
#            g.add_edge(u,v)
#            g.remove_node(i)
#            flag = 0
#
#nx.write_edgelist(g, '/Users/carl/Desktop/vascgraph/graphtest', data=True)
##nx.write_edgelist(sub, '/Users/carl/Desktop/vascgraph/subgraphtest', data=True)


#triplett version
#def _findConnectedCliquesTriplet(c, tocheck, cliques):
#  """Helper to find connected Cliques"""
#  ncl = len(cliques);
#  cc = np.zeros(ncl, dtype = bool);
#  for n in c:
#    cc = np.logical_or([n in c2 for c2 in cliques], cc);
#  cc = np.logical_and(cc, tocheck);
#  cids = np.where(cc)[0];
#  tocheck[cids]= False;
#  
#  if len(cids) > 0:
#    cids2 = [];
#    for ci in cids:
#      cidsadd, tocheck = _findConnectedCliques(cliques[ci], tocheck, cliques);
#      cids2.extend(cidsadd);
#    cids = np.hstack([cids, cids2]);
#  return cids, tocheck;
#
#
#
#
#def cleanGraphTriplet(graph, verbose = True):
#  """Remove all 3 cliques to get pure branch structure"""
#  g = graph.copy();
#  clqs = list(nx.find_cliques(g));
#  l = np.array([len(c) for c in clqs]);
#  rem = np.where(l>=3)[0];
#  if len(rem) == 0:
#    return g;
#  clqs = [clqs[r] for r in rem];
#  ncl = len(clqs);
#  
#  if verbose:
#    print('Graph cleaning: found %d cliques' % ncl)
#
#  # identify meta cliques
#  ctest = np.ones(ncl, dtype = bool);
#  cadj = [];
#  for i, c in enumerate(clqs):
#    if verbose and i % 100 == 0:
#      print('Meta clique detection: %d / %d');
#    if ctest[i]:
#      ctest[i] = False;
#      cids, ctest = _findConnectedCliques(c, ctest, clqs);
#      cadj.append(np.array(np.hstack([[i],cids]), dtype = int));
#  
#  if verbose:
#    print('Graph cleaning: found %d meta cliques' % len(cadj));
#  
#  # remove meta cliques
#  for ci in cadj:
#    cs = [clqs[i] for i in ci];
#    cnds = [];
#    for c in cs:
#      cnds.extend([n for n in c if n not in cnds]);
#      
#    # find all unique neighbours
#    nbs = [];
#    for n in cnds:
#      nbs.extend([nn for nn in g[n] if nn not in cnds and nn not in nbs]);
#    
#    # remove all meta clique nodes
#    g.remove_nodes_from(cnds);
#    
#    # generate new center node
#    new = tuple(np.array(np.mean(np.array(cnds), axis = 0), dtype = int));
#    g.add_node(new);      
#    
#    # reconnect
#    g.add_edges_from([(c, new) for c in nbs]);
#  
#  return g;
