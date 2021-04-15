# -*- coding: utf-8 -*-
"""
Skeleton to Graph transforms

This module provides routines to convert 3d skeletons to graphs, optimize and plot them.
"""
__license__ = 'MIT License <http://www.opensource.org/licenses/mit-license.php>'
__author__ = 'Christoph Kirst <ckirst@rockefeller.edu>'
__docformat__ = 'rest'

import sys
import numpy as np


import ClearMap.Analysis.Graphs.Graph as grp;

import ClearMap.ImageProcessing.Topology.Topology3D as t3d

import ClearMap.DataProcessing.LargeData as ld

import ClearMap.Utils.Timer as Timer;


###############################################################################
### Transform Skeleton to Graph
###############################################################################

def skeletonToGraph(skeleton, points = None, coordinates = True, radii = None, method = 'gt', verbose = False):
  """Converts a binary skeleton image to a graph
    
  Arguments
  ---------
    skeleton : array
        2d/3d binary skeleton image
    points : array
        list of skeleton points (as 1d indices of flat array)to save time for large data sets if calculated previously
    radii : array
        list of radii associated with each vertex
    method : 'gt'
        specify the underlying graph backend, 'gt' for graph_tool is currently the only option.
    verbose : bool
        print progress
    
  Returns
  -------
    object: graph tool object
  """

  # create underlying gt graph
  g = skeletonToGtGraph(skeleton = skeleton, points = points, coordinates = coordinates, radii = radii, verbose = verbose);
  
  # create Graph and meta data
  gr = grp.Graph(graph = g);
  gr.setGraphProperty('shape', skeleton.shape);
  
  return gr


import graph_tool as gt


def skeletonToGtGraph(skeleton, points = None, coordinates = True, radii = None, verbose = False):
  """Converts a binary skeleton image to a graph-tool graph
  
  Arguments:
    skeleton (array): 2d/3d binary skeleton image
    points (array): list of skeleton points (as 1d indices of flat array)to save time for large data sets if calculated previously
    radii (array): list of radii associated with each vertex
    verbose (bool): print progress
    
  Returns:
    object: graph tool object
  """
  
  timer = Timer.Timer();
  timer_all = Timer.Timer();
  
  if points is None:
    points = ld.where(skeleton.reshape(-1, order = 'A'));
    
    if verbose:
      timer.printElapsedTime('Point list generation');
      timer.reset();
  
  #create labeled skeleton array for convolution
  nvertices = points.shape[0];
  #label = np.zeros_like(skeleton, dtype = 'int64', order = 'A');
  #labelflat = label.reshape(-1, order = 'A');
  #labelflat[points] = np.arange(1, nvertices+1);
  #if verbose:
  #    timer.printElapsedTime('Label generation');
  #    timer.reset();
  
  g = gt.Graph(directed = False);
  g.add_vertex(n=nvertices);
  
  if verbose:
    timer.printElapsedTime('Graph initialized with %d vertices' % nvertices);
    timer.reset();
  
  #add edges
  edges_all = np.zeros((0,2), dtype = int);
  for oi,o in enumerate(t3d.orientations()):
    # calculate off set
    offset = np.sum((np.hstack(np.where(o))-[1,1,1]) * skeleton.strides) 
    edges = ld.neighbours(points, offset);
    if len(edges) > 0:
      edges_all = np.vstack([edges_all, edges]);
      
    if verbose:
      timer.printElapsedTime('%d edges with orientation %d/13 found' % (edges.shape[0], oi+1));
      timer.reset();
   
  
  if edges_all.shape[0] > 0:
    g.add_edge_list(edges_all);
    
  if verbose:
    timer.printElapsedTime('Added %d edges to graph' % (edges_all.shape[0]));
    timer.reset();

  if coordinates:
    vx, vy, vz = g.new_vertex_property("int"), g.new_vertex_property("int"), g.new_vertex_property("int");
    order = 'C';
    if skeleton.flags['F_CONTIGUOUS']:
      order = 'F';
    vx.a, vy.a, vz.a = np.unravel_index(points, skeleton.shape, order = order);
    g.vertex_properties['x'] = vx;
    g.vertex_properties['y'] = vy;
    g.vertex_properties['z'] = vz;
  
  if radii is not None:
    vr = g.new_vertex_property("double");
    vr.a = radii;
    g.vertex_properties['r'] = vr;
  
  if verbose:
    timer_all.printElapsedTime('Skeleton to Graph');

  return g;



###############################################################################
### Reduce graph sturcture
###############################################################################

def reduceGraph(graph, remove_isolated_vertices = True, length = False, edge_geometry = False, verbose = False):
  """Reduce graph by replacing all vertices with degree two"""
   
  result = reduceGtGraph(graph.graph, remove_isolated_vertices = remove_isolated_vertices, length = length, edge_geometry = edge_geometry, verbose = verbose);

  if edge_geometry:
    gg, bids, bstart, bend = result;
    g = grp.Graph(graph = gg, name = graph.name);
    
    r = graph.vertexRadius();
    xyz = graph.vertexCoordinates();
    
    g.setEdgeGeometry(coordinates = xyz[bids], radii = r[bids], start = bstart, end = bend);
  
  else:
    g = grp.Graph(graph = result, name = graph.name);
    
  #g.removeSelfLoops();
  
  return g;


def _graphBranch(v0):
  """Follows a branch in a graph of vertices with degree 2"""

  es = [e for e in v0.out_edges()]; 
  assert(len(es) == 2)
  
  vs1, es1, ep1 = _follow_branch(es[0], v0);
  if ep1 != v0:
    vs2, es2, ep2 = _follow_branch(es[1], v0);
    
    vs1.reverse();
    vs1.append(v0);
    vs1.extend(vs2); 

    es1.reverse();
    es1.extend(es2);
    is_isolated_loop = False;
  else:
    ep2 = v0;
    is_isolated_loop = True;
  
  #vs = [int(vv) for vv in vs1];
  #if len(np.unique(vs)) != len(vs):
  #  print 'double %d' % int(v), len(vs), len(np.unique(vs))
  #  print vs
    
  return (vs1, es1, [ep1, ep2], is_isolated_loop);


def _follow_branch(e, v0):
  """follows branch from vertex v in direction of edge e"""
  #argument v0 to prevent infinte loops for isolated closed loops
  edges = [e];
  
  v = e.target();
  vertices = [v];
  
  #counter = 0;
  while v.out_degree() == 2 and v != v0: # and counter < 10000
    for enew in v.out_edges():
      if enew != e:
        break;
    e = enew;
    edges.append(e);
    v = e.target();
    vertices.append(v);
    #counter += 1;
  
  #if counter >= 1000:
  #  print 'large counter: %d' % int(v);
  
  return (vertices, edges, v)
        

def reduceGtGraph(graph, remove_isolated_vertices = True, length = False, edge_geometry = False, verbose = False):
  """Reduce graph by replacing all vertices with degree two"""
  
  timer = Timer.Timer();
  timer_all = Timer.Timer();
  
  #add ids to original graph
  g = graph.copy();
  p = graph.new_vertex_property('int');
  p.a = np.arange(g.num_vertices());
  g.vertex_properties['_id'] = p;
  
  # remove isolated vertices
  if remove_isolated_vertices:
    non_iso = g.get_out_degrees(g.get_vertices());
    non_iso = non_iso > 0;
    gv = gt.GraphView(g, vfilt = non_iso);
    g = gt.Graph(gv, prune = True);
    non_iso = [];
    
    if verbose:
      timer.printElapsedTime('Removed %d isolated nodes' % np.logical_not(non_iso).sum());
      timer.reset();
  
  #find vertices with deg 2    
  bpts = g.get_out_degrees(g.get_vertices());
  bpts = bpts == 2;
  ids  = np.where(bpts)[0];
  bpts = np.logical_not(bpts)
  ids_branch = np.where(bpts)[0];
  nnb = len(ids);
  nb = len(ids_branch);
  nn = len(bpts);
  bpts = [];
    
  if length is True:
    with_edge_length = True;
    if 'l' in g.edge_properties.keys():
      el = g.edge_properties['l']; 
    else:
      el = g.new_edge_property("double");
      el.a += 1;
      g.edge_properties['l'] = el;
  else:
    with_edge_length = False;
  
  if 'r' in graph.edge_properties.keys():
    with_edge_radii = True;
    er = g.edge_properties['r'];
  else:
    with_edge_radii = False;
  
  
  eadd = []; elen = []; erad = [];
  
  if edge_geometry:
    branch_start = [];
    branch_end   = [];
    branch_indices = [];
    index_pos = 0;
  
  if verbose:
    timer.printElapsedTime('Found %d branching and %d non-branching nodes' % (nb, nnb));
    timer.reset();
  
  
  # find direct edges between branch points 
  checked = np.zeros(g.num_edges(), dtype = bool);
  for ii,i in enumerate(ids_branch):
    v = g.vertex(i);
    es = v.out_edges();
    for e in es:
      if e.target().out_degree() != 2:
        eid = g.edge_index[e];
        if not checked[eid]:
          checked[eid] = True;
          eadd.append([e.source(), e.target()]);
          if with_edge_length:
            elen.append(el.a[eid]);
          if with_edge_radii:
            erad.append(er.a[eid]);
          
          if edge_geometry:
            branch_start.append(index_pos);
            index_pos += 2;
            branch_end.append(index_pos);
            branch_indices.extend([int(e.source()), int(e.target())]); 
          
    if verbose and (ii+1) % 250000 == 0:
      timer.printElapsedTime('Scanned %d/%d branching nodes found %d branches' % (ii+1, nb, len(eadd)));
  
  if verbose:
    timer.printElapsedTime('Scanned %d/%d branching nodes found %d branches' % (ii+1, nb, len(eadd)));
  
  
  # reduce longer edges
  checked = np.zeros(nn, dtype = bool);
  for ii,i in enumerate(ids):    
    if not checked[i]:
      # extract branch 
      checked[i] = True;
      v = g.vertex(i);
      vs, es, epts, is_isolated_loop = _graphBranch(v);
      if not is_isolated_loop:
        vids = [int(vv) for vv in vs];
        checked[vids[1:-1]] = True;
        
        eadd.append(epts);
        
        if with_edge_length or with_edge_radii:
          eids = [g.edge_index[e] for e in es];
          if with_edge_length:
            elen.append(np.sum(el.a[eids]));
          if with_edge_radii:
            erad.append(np.mean(er.a[eids]));
        
        if edge_geometry:
          branch_start.append(index_pos);
          index_pos += len(vids);
          branch_end.append(index_pos);
          branch_indices.append(vids);
      
    if verbose and (ii+1) % 250000 == 0:
      timer.printElapsedTime('Scanned %d/%d non-branching nodes found %d branches' % (ii+1, nnb, len(eadd)));     
  
  if verbose:
    timer.printElapsedTime('Scanned %d/%d non-branching nodes found %d branches' % (ii+1, nnb, len(eadd)));

   
  #update graph
  gv = gt.GraphView(g, efilt = np.zeros(g.num_edges()));
  g = gt.Graph(gv, prune = True);
    
  g.add_edge_list(eadd);
  if with_edge_length:
    el = g.edge_properties['l'];
    el.a[:] = elen;
  if with_edge_radii: 
    er = g.edge_properties['r'];
    er.a[:] = erad;
    
  if edge_geometry:
    branch_indices = np.hstack(branch_indices);
    branch_start = np.array(branch_start);
    branch_end   = np.array(branch_end);

    old_indices = g.vertex_properties['_id'].a
    branch_indices = old_indices[branch_indices];    
    del g.vertex_properties['_id'];   
    
  # remove non-branching points 
  gv = gt.GraphView(g, vfilt = np.logical_not(checked));
  g = gt.Graph(gv, prune = True); 
  
  if verbose:
    timer_all.printElapsedTime('Graph reduced from %d to %d nodes and %d to %d edges' % (graph.num_vertices(), g.num_vertices(), graph.num_edges(), g.num_edges()));

  if edge_geometry:
    return g, branch_indices, branch_start, branch_end;
  else:
    return g;


###############################################################################
### Clean Graph
###############################################################################

def cleanGraph(graph, verbose = True):
  """Remove all cliques to get pure branch structure"""
  
  gg = cleanGtGraph(graph.graph, verbose = verbose);
  g = grp.Graph(graph = gg, name = graph.name);
  g.removeSelfLoops();
  
  return g;
  

import graph_tool.topology as gtop

def _groupArray(array):
  idx = np.argsort(array);
  arr = array[idx];
  dif = np.diff(arr);
  res = np.split(idx, np.where(dif > 0)[0]+1);
  return res


def cleanGtGraph(graph, verbose = True):
  """Remove all cliques to get pure branch structure"""
  
  timer = Timer.Timer();
  timer_all = Timer.Timer();
  
  # find branch points
  deg = graph.get_out_degrees(graph.get_vertices());
  ids = deg >= 3;
  nbp = ids.sum();
  nvertices = len(deg);
  
  if verbose:
    timer.printElapsedTime('Graph cleaning: found %d branch points' % nbp);
    timer.reset();
  
  if nbp == 0:
    return graph.copy();
  
  # find cliques = connected components of branch points
  vfilt = np.zeros(nvertices, dtype = bool);
  vfilt[ids] = True;
  
  gsub = gt.GraphView(graph, vfilt);
  comp, hist = gtop.label_components(gsub);
  comp = np.array(comp.a); 
  comp[vfilt] += 1;
  comp = _groupArray(comp);
  comp = comp[1:];  
  clique_ids = np.where(hist > 1)[0];
  ncliques = len(clique_ids);
  
  if verbose:
    timer.printElapsedTime('Graph cleaning: detected %d cliques of branch points' % ncliques);
    timer.reset();
  
  # remove meta cliques
  vfilt = np.ones(nvertices + ncliques, dtype = bool)
  g = graph.copy();
  g.add_vertex(ncliques);
  
  if 'x' in g.vertex_properties.keys():
    with_xyz = True;
    vx = g.vertex_properties['x']  
    vy = g.vertex_properties['y'] 
    vz = g.vertex_properties['z'] 
  else:
    with_xyz = False;
    
  if 'r' in g.vertex_properties.keys():
    with_vertex_radii = True;
    vr = g.vertex_properties['r'];
  
  #if 'r' in g.edge_properties.keys():
  #  with_edge_radii = True;
  #  er = g.edge_properties['r'];  
  
  
  for i,ci in enumerate(clique_ids):
    cc = comp[ci];
    
    #remove clique vertices 
    vfilt[cc] = False;
    
    # get neighbours 
    nbs = np.hstack([graph.get_out_neighbours(c) for c in cc])
    nbs = np.setdiff1d(np.unique(nbs), cc);
    #print [[ni, i+nvertices] for ni in nbs]
    g.add_edge_list([[ni, i+nvertices] for ni in nbs]);
    
    if with_xyz:
      vx.a[i+nvertices] = np.mean(vx.a[cc]);
      vy.a[i+nvertices] = np.mean(vy.a[cc]);
      vz.a[i+nvertices] = np.mean(vz.a[cc]);
    
    if with_vertex_radii:
      vr.a[i+nvertices] = np.max(vr.a[cc]);
    
    if verbose and i % 10000 == 0:
      timer.printElapsedTime('Graph cleaning: reducing %d / %d' % (i, ncliques));
  
  #generate new graph
  gv = gt.GraphView(g, vfilt = vfilt);
  g = gt.Graph(gv, prune = True);
    
  if verbose:
    timer_all.printElapsedTime('Graph cleaning: removed %d cliques of branch points from %d to %d nodes and %d to %d edges' % (ncliques, graph.num_vertices(), g.num_vertices(), graph.num_edges(), g.num_edges()));
  
  return g;




### Tests

if __name__ == "__main__":
  #from importlib import reload
  import numpy as np
  import matplotlib.pyplot as plt
  
  import ClearMap.Analysis.Graphs.SkeletonGraph as sg
  reload(sg)
  
  data = np.load('./ClearMap/Test/Skeletonization/test_bin.npy');
  skel = np.load('./ClearMap/Test/Skeletonization/test_skel.npy');
  radii = np.load('./ClearMap/Test/Skeletonization/test_radii.npy');
  pts = np.load('./ClearMap/Test/Skeletonization/test_pts.npy');
  
  #mlab.figure(bgcolor=(1,1,1));
  #mlab.contour3d(np.array(data, dtype = int), color = (0.5, 0.5, 0.5), opacity = 0.2);
  #mlab.contour3d(np.array(skel, dtype = int), color = (1,0,0), opacity = 0.2);  
  
  radii[radii > 1000] = 1000;
  
  reload(sg)
  g = sg.skeletonToGtGraph(skel, radii = radii, verbose = True);
  
  
  #get branch points
  vs = [v for v in g.vertices() if v.out_degree() > 2]
  branches = np.zeros_like(skel);
  for v in vs:
    x,y,z = g.vertex_properties['x'][v],g.vertex_properties['y'][v],g.vertex_properties['z'][v]
    branches[x,y,z] = 1;
  
  
  dv.plot([data, skel, branches])
  
  sg.plotGtGraph3d(g)
  
  #mlab.figure(bgcolor = (1,1,1))
  #sg.plotGtGraph3d(g)
  #mlab.contour3d(np.array(data, dtype = int), color = (0.5, 0.5, 0.5), opacity = 0.2);
  
  
  gc = sg.cleanGtGraph(g)
  
  gr = sg.reduceGtGraph(gc, verbose = True, length = True)
  
  sg.plotGtGraph3d(gr)
  sg.plotGtGraph3d(gc)
  
  #mlab.figure(bgcolor = (1,1,1))
  #sg.plotGtGraph3d(gr, feature = 'l')




### Graph-tool Routines
#
#def skeletonToGtGraph(skeleton, points = None, coordinates = True, radii = None, verbose = False):
#  """Converts a binary skeleton image to a graph-tool graph
#  
#  Arguments:
#    skeleton (array): 2d/3d binary skeleton image
#    points (array): list of skeleton points (as 1d indices of flat array)to save time for large data sets if calculated previously
#    radii (array): list of radii associated with each vertex
#    verbose (bool): print progress
#    
#  Returns:
#    object: graph tool object
#  """
#  
#  timer = Timer.Timer();
#  timer_all = Timer.Timer();
#  
#  if points is None:
#    points = ld.where(skeleton.reshape(-1, order = 'A'));
#    
#    if verbose:
#        timer.printElapsedTime('Point list generation');
#        timer.reset();
#  
#  #create labeled skeleton array for convolution
#  nvertices = points.shape[0];
#  label = np.zeros_like(skeleton, dtype = 'int64', order = 'A');
#  labelflat = label.reshape(-1, order = 'A');
#  labelflat[points] = np.arange(1, nvertices+1);
#  
#  if verbose:
#      timer.printElapsedTime('Label generation');
#      timer.reset();
#  
#  
#  g = gt.Graph();
#  g = gt.Graph(directed = False);
#  g.add_vertex(n=nvertices);
#  
#  if verbose:
#    timer.printElapsedTime('Graph initialized with %d vertices' % nvertices);
#    timer.reset();
#  
#  #add edges
#  for oi,o in enumerate(top.orientations()):
#    nhd = cnv.convolve3DIndex(label, np.array(o, dtype = label.dtype), points);
#    ids = ld.where(nhd);
#    if len(ids) > 0:
#      g.add_edge_list(np.array([ids, nhd[ids]-1]).T);
#      
#    if verbose:
#      timer.printElapsedTime('%d edges with orientation %d/12 added' % (len(ids), oi));
#      timer.reset();
#
#  if coordinates:
#    vx, vy, vz = g.new_vertex_property("int"), g.new_vertex_property("int"), g.new_vertex_property("int");
#    order = 'C';
#    if skeleton.flags['F_CONTIGUOUS']:
#      order = 'F';
#    vx.a, vy.a, vz.a = np.unravel_index(points, skeleton.shape, order = order);
#    g.vertex_properties['x'] = vx;
#    g.vertex_properties['y'] = vy;
#    g.vertex_properties['z'] = vz;
#  
#  if radii is not None:
#    vr = g.new_vertex_property("double");
#    vr.a = radii;
#    g.vertex_properties['r'] = vr;
#    
#  if verbose:
#    timer_all.printElapsedTime('Skeleton to Graph');
#
#  return g;

#
#
#### Graph-tool Routines
#
#def skeletonToGtGraph(skeleton, points = None, coordinates = True, radii = None, verbose = False):
#  """Converts a binary skeleton image to a graph-tool graph
#  
#  Arguments:
#    skeleton (array): 2d/3d binary skeleton image
#    points (array): list of skeleton points (as 1d indices of flat array)to save time for large data sets if calculated previously
#    radii (array): list of radii associated with each vertex
#    verbose (bool): print progress
#    
#  Returns:
#    object: graph tool object
#  """
#  
#  timer = Timer.Timer();
#  timer_all = Timer.Timer();
#  
#  if points is None:
#    points = ld.where(skeleton.reshape(-1, order = 'A'));
#    
#    if verbose:
#      timer.printElapsedTime('Point list generation');
#      timer.reset();
#  
#  #create labeled skeleton array for convolution
#  nvertices = points.shape[0];
#  #label = np.zeros_like(skeleton, dtype = 'int64', order = 'A');
#  #labelflat = label.reshape(-1, order = 'A');
#  #labelflat[points] = np.arange(1, nvertices+1);
#  #if verbose:
#  #    timer.printElapsedTime('Label generation');
#  #    timer.reset();
#  
#  
#  g = gt.Graph();
#  g = gt.Graph(directed = False);
#  g.add_vertex(n=nvertices);
#  
#  if verbose:
#    timer.printElapsedTime('Graph initialized with %d vertices' % nvertices);
#    timer.reset();
#  
#  #add edges
#  for oi,o in enumerate(top.orientations()):
#    # calculate off set
#    offset = np.sum((np.hstack(np.where(o))-[1,1,1]) * skeleton.strides) 
#    edges  = ld.neighbours(points, offset);
#      
#    if edges.shape[0] > 0:
#      g.add_edge_list(edges);
#      
#    if verbose:
#      timer.printElapsedTime('%d edges with orientation %d/12 added' % (edges.shape[0], oi));
#      timer.reset();
#
#  if coordinates:
#    vx, vy, vz = g.new_vertex_property("int"), g.new_vertex_property("int"), g.new_vertex_property("int");
#    order = 'C';
#    if skeleton.flags['F_CONTIGUOUS']:
#      order = 'F';
#    vx.a, vy.a, vz.a = np.unravel_index(points, skeleton.shape, order = order);
#    g.vertex_properties['x'] = vx;
#    g.vertex_properties['y'] = vy;
#    g.vertex_properties['z'] = vz;
#  
#  if radii is not None:
#    vr = g.new_vertex_property("double");
#    vr.a = radii;
#    g.vertex_properties['r'] = vr;
#    
#  if verbose:
#    timer_all.printElapsedTime('Skeleton to Graph');
#
#  return g;
#
#
#
#def plotGtGraph3d(graph, zfactor = 1.0, feature = 'r', line_width = 2, opacity = .9, cmap = 'jet', colorbar = True):
#  """Plot a 3d graph of the skeleton
#  
#  Arguments:
#    graph (object): the graph to plot
#    zfactor (float): scaling in z direction
#    feature (string): the feature to color edges / vertices in
#    line_width (float): line width of edges
#    opacity (float): opacity of the edges
#    cmap (string): color map
#    colorbar (bool): add color bar to the plot
#  
#  Returns:
#    mayavi scence
#  """
#  
#  if feature in graph.edge_properties.keys():
#    return plotGtGraph3dEdges(graph, zfactor = zfactor, feature = feature, line_width = line_width, opacity = opacity, cmap = cmap, colorbar = colorbar);
#  else:
#    return plotGtGraph3dVertices(graph, zfactor = zfactor, feature = feature, line_width = line_width, opacity = opacity, cmap = cmap, colorbar = colorbar);
#
#
#def plotGtGraph3dVertices(graph, zfactor = 1.0, feature = 'r', line_width = 2, opacity = .9, cmap = 'jet', colorbar = True):
#  """Plot a 3d graph of the skeleton ysing the vertex radii for coloring
#  
#  Arguments:
#    graph (object): the graph to plot
#    zfactor (float): scaling in z direction
#    line_width (float): line width of edges
#    opacity (float): opacity of the edges
#    cmap (string): color map
#    colorbar (bool): add color bar to the plot
#  
#  Returns:
#    mayavi scence
#  """
#  # get graph positions
#  x = np.array(graph.vertex_properties['x'].a, dtype = 'int32');
#  y = np.array(graph.vertex_properties['y'].a, dtype = 'int32');
#  z = np.array(graph.vertex_properties['z'].a, dtype = 'int32');
#  
#  with_radii = feature in graph.vertex_properties.keys();
#  if with_radii:
#    scalars = np.array(graph.vertex_properties[feature].a, dtype = 'float32');
#  else:
#    scalars = np.ones(x.shape[0], dtype = 'float32');
#
#  pts = mlab.pipeline.scalar_scatter(x, y, zfactor * z, scalars)
#  pts.mlab_source.dataset.lines = np.array(graph.get_edges()[:,:2], dtype = 'int32')
#  pts.update()
#  
#  lines = mlab.pipeline.stripper(pts);
#  mlab.pipeline.surface(lines, colormap = cmap, line_width = line_width, opacity = opacity)
#  
#  if with_radii and colorbar is True:
#      mlab.colorbar(orientation = 'vertical', title='Radius [pix]');    
#      
#  mlab.axes()
#  
#  return lines
#
#
#def plotGtGraph3dEdges(graph, zfactor = 1.0, feature = 'r', line_width = 2, opacity = .9, cmap = 'jet', colorbar = True):
#  """Plot a 3d graph of the skeleton using a edge width property for coloring
#  
#  Arguments:
#    graph (object): the graph to plot
#    zfactor (float): scaling in z direction
#    line_width (float): line width of edges
#    opacity (float): opacity of the edges
#    cmap (string): color map
#    colorbar (bool): add color bar to the plot
#  
#  Returns:
#    mayavi scence
#  """
#  # get graph positions
#  x = np.array(graph.vertex_properties['x'].a, dtype = 'int32');
#  y = np.array(graph.vertex_properties['y'].a, dtype = 'int32');
#  z = np.array(graph.vertex_properties['z'].a, dtype = 'int32');
#  
#  with_radii = feature in graph.edge_properties.keys();
#  if with_radii:
#    scalars = np.array(graph.edge_properties[feature].a, dtype = 'float32');
#  else:
#    scalars = np.ones(graph.num_edges(), dtype = 'float32');
#  
#  edges = np.array(graph.get_edges(), dtype = 'int32');
#  nedges = edges.shape[0];
#  i, j = edges[:, 0], edges[:, 1];
#  x = np.hstack([x[i], x[j]]);
#  y = np.hstack([y[i], y[j]]);
#  z = np.hstack([z[i], z[j]]);
#  c = np.hstack([scalars, scalars]);
#  connections = np.vstack([np.arange(nedges), np.arange(nedges,2*nedges)]).T
#  
#  pts = mlab.pipeline.scalar_scatter(x, y, zfactor * z, c)
#  pts.mlab_source.dataset.lines = connections;
#  pts.update();
#  
#  lines = mlab.pipeline.stripper(pts);
#  surf = mlab.pipeline.surface(lines, colormap = cmap, line_width = line_width, opacity = opacity)
#  
#  if with_radii and colorbar is True:
#      mlab.colorbar(orientation = 'vertical', title='Radius [pix]');    
#      
#  mlab.axes()
#  
#  return surf
#
#
#
#def plotGtGraph3dTubes(graph, zfactor = 1.0, feature = 'r', line_width = 2, opacity = .9, cmap = 'jet', colorbar = True):
#  """Plot a 3d graph of the skeleton using feature for the radii
#  
#  Arguments:
#    graph (object): the graph to plot
#    zfactor (float): scaling in z direction
#    line_width (float): line width of edges
#    opacity (float): opacity of the edges
#    cmap (string): color map
#    colorbar (bool): add color bar to the plot
#  
#  Returns:
#    mayavi scence
#  """
#  # get graph positions
#  x = np.array(graph.vertex_properties['x'].a, dtype = 'int32');
#  y = np.array(graph.vertex_properties['y'].a, dtype = 'int32');
#  z = np.array(graph.vertex_properties['z'].a, dtype = 'int32');
#
#  with_radii = feature in graph.vertex_properties.keys();
#  if with_radii:
#    scalars = np.array(graph.vertex_properties[feature].a, dtype = 'float32');
#  else:
#    scalars = np.ones(x.shape[0], dtype = 'float32');
#
#  pts = mlab.pipeline.scalar_scatter(x, y, zfactor * z, scalars)
#  pts.mlab_source.dataset.lines = np.array(graph.get_edges()[:,:2], dtype = 'int32')
#  pts.update()
#
#  #pts = mlab.points3d(x, y, z, 1.5 * scalars.max() - scalars,scale_factor=0.015, resolution=10)
#  #pts.mlab_source.dataset.lines = np.array(connections)
#
#  # Use a tube fiter to plot tubes on the link, varying the radius with the scalar value
#  tube = mlab.pipeline.tube(pts, tube_radius=0.15)
#  tube.filter.radius_factor = 1.
#  tube.filter.vary_radius = 'vary_radius_by_scalar'
#  surf = mlab.pipeline.surface(tube, colormap = cmap, opacity = opacity);
#
#  # Visualize the local branch density
#  #mlab.pipeline.volume(mlab.pipeline.gaussian_splatter(pts))
#  
#  return surf;
#  
#  
#  
#  
#  # get graph positions
#  x = np.array(graph.vertex_properties['x'].a, dtype = 'int32');
#  y = np.array(graph.vertex_properties['y'].a, dtype = 'int32');
#  z = np.array(graph.vertex_properties['z'].a, dtype = 'int32');
#  
#  with_radii = feature in graph.vertex_properties.keys();
#  if with_radii:
#    scalars = np.array(graph.vertex_properties[feature].a, dtype = 'float32');
#  else:
#    scalars = np.ones(x.shape[0], dtype = 'float32');
#
#  pts = mlab.pipeline.scalar_scatter(x, y, zfactor * z, scalars)
#  pts.mlab_source.dataset.lines = np.array(graph.get_edges()[:,:2], dtype = 'int32')
#  pts.update()
#  
#  lines = mlab.pipeline.stripper(pts);
#  mlab.pipeline.surface(lines, colormap = cmap, line_width = line_width, opacity = opacity)
#  
#  if with_radii and colorbar is True:
#      mlab.colorbar(orientation = 'vertical', title='Radius [pix]');    
#      
#  mlab.axes()
#  
#  return lines
#
#
#
#
#
#### Reduce graph sturcture
#
#import sys
#
#def _graphBranch(v):
#  """Follows a branch in a graph of vertices with degree 2"""
#  vs = [v];
#  es = [e for e in v.out_edges()];
#  ls = [];
#
#  #if len(es) != 2: 
#  #  raise RuntimeWarning('more than two edges on vertex');
#  
#  vs, es, ls = _follow_branch_2(v, es[0], v, vs, es, ls);
#  vs, es, ls = _follow_branch_2(v, es[1], v, vs, es, ls);  
#  
#  return (vs, es, ls);
#
#def _follow_branch_2(v, e, v0, vertices, edges, endpoints):
#  vnew = e.target();
#  if vnew == v0:
#    return (vertices, edges, endpoints);
#  if vnew.out_degree() == 2:
#    vertices.append(vnew);
#    for enew in vnew.out_edges():
#      if enew != e:
#        edges.append(enew);
#        vertices, edges, endpoints = _follow_branch_2(vnew, enew, v0, vertices, edges, endpoints);
#        break
#  else:
#    endpoints.append(vnew);
#  return (vertices, edges, endpoints)
#        
#
#def addEdgeRadiusGtGraph(graph):
#  """Add edge radius property to the graph given radii of vertices"""
#  
#  if 'r' not in graph.vertex_properties.keys():
#    raise RuntimeError('no radii information on vertices to calculate edge radii');
#  
#  er = graph.new_edge_property("double");
#  graph.edge_properties['r'] = er;
#   
#  vr = graph.vertex_properties['r'];
#  
#  edges = graph.get_edges();
#  i, j = edges[:,0], edges[:,1];
#  er.a = 0.5 * (vr.a[i] + vr.a[j]);
#  
#  return graph;
#
#
#def removeIsolatedVerticesGtGraph(graph, verbose = True):
#  
#  timer = Timer.Timer();
#  degs = graph.get_out_degrees(graph.get_vertices());
#  degs = degs > 0;
#  gv = gt.GraphView(graph, vfilt = degs);
#  g = gt.Graph(gv, prune = True);
#  if verbose:
#    timer.printElapsedTime('Removed %d isolated nodes' % np.logical_not(degs).sum());
#    timer.reset();
#    
#  return g;
#  
#
#
#def reduceGtGraph(graph, removeIsolatedVertices = True, verbose = False, length = False):
#  """Reduce graph by replacing all vertices with degree two"""
#  
#  timer = Timer.Timer();
#  timer_all = Timer.Timer();
#  
#  if removeIsolatedVertices:
#    degs = graph.get_out_degrees(graph.get_vertices());
#    degs = degs > 0;
#    gv = gt.GraphView(graph, vfilt = degs);
#    g = gt.Graph(gv, prune = True);
#    if verbose:
#      timer.printElapsedTime('Removed %d isolated nodes' % np.logical_not(degs).sum());
#      timer.reset();
#    
#  else:
#    g = graph.copy();
#  
#  #find vertices with deg 2    
#  degs = g.get_out_degrees(g.get_vertices());
#  ids = np.where(degs == 2)[0];
#  nrem = len(ids);
#  checked = np.zeros(len(degs), dtype = bool);
#
#  if verbose:
#    timer.printElapsedTime('Found %d non-branching nodes' % nrem);
#    timer.reset();
#    
#  if length is True:
#    with_edge_length = True;
#    el = g.new_edge_property("double");
#    el.a += 1;
#    g.edge_properties['l'] = el;
#  else:
#    with_edge_length = False;
#   
#  if 'r' in graph.edge_properties.keys():
#    with_edge_radii = True;
#    er = g.edge_properties['r'];
#  else:
#    with_edge_radii = False;
#  
#  branchLimit = g.num_vertices();
#  reclimit = sys.getrecursionlimit();
#  sys.setrecursionlimit(branchLimit);
#
#  #vrem = []; erem = []; 
#  eadd = []; elength = []; erad = [];
#
#  nbranches = 0;
#  for ii,i in enumerate(ids):    
#    if not checked[i]:
#      # extract branch 
#      checked[i] = True;
#      nbranches += 1;
#      v = g.vertex(i);
#      vs, es, epts = _graphBranch(v);
#      vids = [int(vv) for vv in vs];
#      checked[vids] = True;
#      #vrem.extend(vs); erem.extend(es); 
#      eadd.append(epts);
#      
#      if with_edge_length:
#        elength.append(len(vids));
#        #elength.append(np.sum(el.a[vids]));
#        
#      if with_edge_radii:
#        eids = [g.edge_index[e] for e in es];
#        erad.append(np.mean(er.a[eids]));
#        #elength.append(np.sum(el.a[vids]));
#  
#    if verbose and ii % 250000 == 0:
#      timer.printElapsedTime('Scanned %d/%d non-branching nodes' % (ii, nrem));
#  
#  if verbose:
#    timer.printElapsedTime('Detected %d branches' % nbranches);
#  
#  #update graph
#  g.add_edge_list(eadd);
#  if with_edge_length:
#    el.a[-len(eadd):] = elength;
#  if with_edge_radii:
#    er.a[-len(eadd):] = erad;
#  
#  gv = gt.GraphView(g, vfilt = np.logical_not(checked));
#  g = gt.Graph(gv, prune = True);
#  
#  sys.setrecursionlimit(reclimit)
#  
#  if verbose:
#    timer_all.printElapsedTime('Graph reduced from %d to %d nodes and %d to %d edges' % (graph.num_vertices(), g.num_vertices(), graph.num_edges(), g.num_edges()));
#  
#  return g;
#
#
#
#
#import graph_tool.topology as gtop
#
#def _groupArray(array):
#  idx = np.argsort(array);
#  arr = array[idx];
#  dif = np.diff(arr);
#  res = np.split(idx, np.where(dif > 0)[0]+1);
#  return res
#
#
#def cleanGtGraph(graph, verbose = True):
#  """Remove all 3 cliques to get pure branch structure"""
#  
#  timer = Timer.Timer();
#  timer_all = Timer.Timer();
#  
#  # find branch points
#  deg = graph.get_out_degrees(graph.get_vertices());
#  ids = deg >= 3;
#  nbp = ids.sum();
#  nvertices = len(deg);
#  
#  if verbose:
#    timer.printElapsedTime('Graph cleaning: found %d branch points' % nbp);
#    timer.reset();
#  
#  if nbp == 0:
#    return graph.copy();
#  
#  # find cliques = connected components of branch points
#  vfilt = np.zeros(nvertices, dtype = bool);
#  vfilt[ids] = True;
#  
#  gsub = gt.GraphView(graph, vfilt);
#  comp, hist = gtop.label_components(gsub);
#  comp = np.array(comp.a); 
#  comp[vfilt] += 1;
#  comp = _groupArray(comp);
#  comp = comp[1:];  
#  clique_ids = np.where(hist > 1)[0];
#  ncliques = len(clique_ids);
#  
#  if verbose:
#    timer.printElapsedTime('Graph cleaning: detected %d cliques of branch points' % ncliques);
#    timer.reset();
#  
#  # remove meta cliques
#  vfilt = np.ones(nvertices + ncliques, dtype = bool)
#  g = graph.copy();
#  g.add_vertex(ncliques);
#  
#  if 'x' in g.vertex_properties.keys():
#    with_xyz = True;
#    vx = g.vertex_properties['x']  
#    vy = g.vertex_properties['y'] 
#    vz = g.vertex_properties['z'] 
#  else:
#    with_xyz = False;
#    
#  if 'r' in g.vertex_properties.keys():
#    with_vertex_radii = True;
#    vr = g.vertex_properties['r'];
#  
#  #if 'r' in g.edge_properties.keys():
#  #  with_edge_radii = True;
#  #  er = g.edge_properties['r'];  
#  
#  
#  for i,ci in enumerate(clique_ids):
#    cc = comp[ci];
#    
#    #remove clique vertices 
#    vfilt[cc] = False;
#    
#    # get neighbours 
#    nbs = np.hstack([graph.get_out_neighbours(c) for c in cc])
#    nbs = np.setdiff1d(np.unique(nbs), cc);
#    #print [[ni, i+nvertices] for ni in nbs]
#    g.add_edge_list([[ni, i+nvertices] for ni in nbs]);
#    
#    if with_xyz:
#      vx.a[i+nvertices] = np.mean(vx.a[cc]);
#      vy.a[i+nvertices] = np.mean(vy.a[cc]);
#      vz.a[i+nvertices] = np.mean(vz.a[cc]);
#    
#    if with_vertex_radii:
#      vr.a[i+nvertices] = np.max(vr.a[cc]);
#    
#    if verbose and i % 10000 == 0:
#      timer.printElapsedTime('Graph cleaning: reducing %d / %d' % (i, ncliques));
#  
#  #generate new graph
#  gv = gt.GraphView(g, vfilt = vfilt);
#  g = gt.Graph(gv, prune = True);
#    
#  if verbose:
#    timer_all.printElapsedTime('Graph cleaning: removed %d cliques of branch points from %d to %d nodes and %d to %d edges' % (ncliques, graph.num_vertices(), g.num_vertices(), graph.num_edges(), g.num_edges()));
#  
#  return g;
#
#
#
#### Skeleton to NetwrokX graph
#
#def skeletonToNxGraph(skeleton, convertIndices = True, verbose = False):
#  """Converts a binary skeleton image to a networkx graph
#  
#  Arguments:
#    skeleton (array): 2d/3d binary skeleton image
#    values (array): 2d/3d array of same size as skeleton with values for the skeleton points (e.g. radius data)
#    convertIndices (bool): if True convert the nodes to have integer labels adn xyz attributes
#    verbose (bool): if True print info during proecssing
#    
#  Returns:
#    graph: networkx graph object
#  """
#  
#  x,y,z,nh = top.neighbourhoodList(skeleton);
#  ids = np.transpose([x,y,z]);
#  nh = np.reshape(nh,(nh.shape[0], 3,3,3), order = 'F')
#  nh[:, 1, 1, 1] = -1;
#  if verbose:
#    print('Graph: information extracted from skeleton...')
# 
#  if len(ids) == 0:
#     return nx.Graph();
#  elif len(ids) == 1:
#    adj = {};
#    adj[tuple(ids[0])] = [];
#    g = nx.from_dict_of_lists(adj);
#  else:
#    g = nx.Graph();
#    for i,pos in enumerate(ids):
#      if verbose and i %1000 == 0:
#          print('Graph: %d/%d nodes constructed...' % (i, len(ids)));
#      p = tuple(pos);
#      g.add_node(p);
#      posnh = np.where(nh[i]>=0);
#      for pp in np.transpose(posnh):
#          g.add_edge(p, tuple(pp+pos-1));
#  
#  if convertIndices is True:
#    g = nx.convert_node_labels_to_integers(g, label_attribute = 'xyz');
#  
#  return g;
#
#
#
#def plotNxGraph3d(graph, radii = None,  cmap = 'jet', line_width = 2, opacity = .9, zfactor = 1.0, colorbar = True):
#  """Plot a 3d graph of the skeleton
#  
#  Arguments:
#    radii: radii of the edges used in color code, if None uniform color
#    
#  Returns:
#    mayavi scence
#  """
#  # get graph positions
#  xyz = nx.get_node_attributes(graph, 'xyz');
#  xyz = np.array([xyz[n] for n in graph.nodes()], dtype = 'int32');
#  #g2 = nx.convert_node_labels_to_integers(graph);
#
#  # scalar colors
#  if radii is not None:
#    scalars = np.array([radii[tuple(x)] for x in xyz], dtype = 'float32');
#  else:
#    #scalars = np.arange(5, xyz.shape[0]+5);
#    scalars = np.ones(xyz.shape[0], dtype = 'float32');
#  
#  #pts = mlab.points3d(xyz[:,0], xyz[:,1], xyz[:,2],
#  #                    scalars,
#  #                    scale_factor=node_size,
#  #                    scale_mode='none',
#  #                    colormap=graph_colormap,
#  #                    resolution=20)
#
#  pts = mlab.pipeline.scalar_scatter(xyz[:,0], xyz[:,1], zfactor * xyz[:,2], scalars)
#  
#  pts.mlab_source.dataset.lines = np.array(graph.edges(), dtype = 'int32')
#  pts.update()    
#  
#  #tube = mlab.pipeline.tube(pts, tube_radius=edge_size)
#  #lab.pipeline.surface(tube, color=edge_color)
#  
#  lines = mlab.pipeline.stripper(pts);
#  mlab.pipeline.surface(lines, colormap = cmap, line_width = line_width, opacity = opacity)
#  
#  if radii is not None and colorbar is True:
#      mlab.colorbar(orientation = 'vertical', title='Radius [pix]');    
#      
#  mlab.axes()
#  
#  return lines
#
#def reduceNxGraph(graph):
#  """Reduce graph by replacing all vertices with degree two"""
#  
#  #TODO: include path length and width estimate on each link
#  g = graph.copy();
#  flag = 0
#  while flag == 0:
#    flag = 1
#    for i in g.copy():
#      #if g.degree(i) < 2:
#      #  g.remove_node(i)
#      #  flag = 0
#      if g.degree(i) == 2:
#        u,v = g[i].keys();
#        g.add_edge(u,v)
#        g.remove_node(i)
#        flag = 0;
#  return g;
#
#
#def cleanNxGraph(graph, verbose = True):
#  """Remove all 3 cliques to get pure branch structure"""
#  g = graph.copy();
#  deg = np.array([g.degree(i) for i in g]);
#  ids = np.where(deg >= 3)[0];
#  
#  if verbose:
#    print('Graph cleaning: found %d branch points' % len(ids));
#
#  
#  if len(ids) ==0 :
#    return g;
#
#  gsub = nx.subgraph(g, ids);
#  cc = list(nx.connected_components(gsub));
#  l = np.array([len(c) for c in cc]);
#  cids = np.where(l > 1)[0];
#  cadj = [cc[i] for i in cids];
#  
#  if verbose:
#    print('Graph cleaning: found %d meta branch points' % len(cids));
#
#  
#  # remove meta cliques
#  nnodes = g.number_of_nodes();
#  for i,ci in enumerate(cadj):
#    if verbose and i % 1000 == 0:
#      print('Graph cleaning: reducing %d / %d' % (i, len(cadj)));
#      
#    # get nodes coordinates
#    xyz = np.array([g.node[i]['xyz'] for i in ci]);
#    
#    # get neighbours 
#    nbs = np.hstack([g.adj[n].keys() for n in ci]);
#    nbs = np.setdiff1d(np.unique(nbs), ci);
#    
#    # remove all meta clique nodes
#    g.remove_nodes_from(ci);
#    
#    # generate new center node
#    newxyz = tuple(np.array(np.mean(xyz, axis = 0), dtype = int));
#    g.add_node(nnodes, xyz = newxyz);   
#    nnodes+=1;
#    
#    # reconnect
#    g.add_edges_from([(c, nnodes) for c in nbs]);
#  
#  return g;
#
#
#
#### Tests
#
#if __name__ == "__main__":
#  #from importlib import reload
#  import numpy as np
#  import matplotlib.pyplot as plt
#  from mayavi import mlab
#  
#  
#  import ClearMap.ImageProcessing.Skeletonization.SkeletonGraph as sg
#  reload(sg)
#  
#  data = np.load('./ClearMap/Test/Skeletonization/test_bin.npy');
#  skel = np.load('./ClearMap/Test/Skeletonization/test_skel.npy');
#  radii = np.load('./ClearMap/Test/Skeletonization/test_radii.npy');
#  pts = np.load('./ClearMap/Test/Skeletonization/test_pts.npy');
#  
#  #mlab.figure(bgcolor=(1,1,1));
#  #mlab.contour3d(np.array(data, dtype = int), color = (0.5, 0.5, 0.5), opacity = 0.2);
#  #mlab.contour3d(np.array(skel, dtype = int), color = (1,0,0), opacity = 0.2);  
#  
#  radii[radii > 1000] = 1000;
#  
#  reload(sg)
#  g = sg.skeletonToGtGraph(skel, radii = radii, verbose = True);
#  
#  
#  #get branch points
#  vs = [v for v in g.vertices() if v.out_degree() > 2]
#  branches = np.zeros_like(skel);
#  for v in vs:
#    x,y,z = g.vertex_properties['x'][v],g.vertex_properties['y'][v],g.vertex_properties['z'][v]
#    branches[x,y,z] = 1;
#  
#  
#  dv.plot([data, skel, branches])
#  
#  sg.plotGtGraph3d(g)
#  
#  #mlab.figure(bgcolor = (1,1,1))
#  #sg.plotGtGraph3d(g)
#  #mlab.contour3d(np.array(data, dtype = int), color = (0.5, 0.5, 0.5), opacity = 0.2);
#  
#  
#  gc = sg.cleanGtGraph(g)
#  
#  gc = sg.addEdgeRadiusGtGraph(gc);
#  
#  gr = sg.reduceGtGraph(gc, verbose = True, length = True)
#  
#  sg.plotGtGraph3d(gr)
#  sg.plotGtGraph3d(gc)
#  
#  mlab.figure(bgcolor = (1,1,1))
#  sg.plotGtGraph3d(gr, feature = 'l')
