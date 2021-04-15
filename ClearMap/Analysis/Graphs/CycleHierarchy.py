# -*- coding: utf-8 -*-
"""
Module to analyse graphs using cycle-hierarchies.

Code based on: Modes et al PRX 2016

Note:
  Matrix versoin for now, can be highly optimized
"""

#import ClearMap.Analysis.Graphs.Graph
import numpy as np

import graph_tool as gt
import graph_tool.search as gts

class BFSSearch(gts.BFSVisitor):
  def __init__(self, graph, root):
    self.pre = graph.new_vertex_property('int64');
    self.dist = graph.new_vertex_property('int64');
    self.edges = [];
    self.cycles = None;
    self.tree = None;
  
  def tree_edge(self, e):
    self.pre[e.target()] = e.source();
    self.dist[e.target()] += self.dist[e.source()] + 1;
    self.edges.append(e);


def spanning_tree(graph, root):
  span = BFSSearch(graph, root);
  gts.bfs_search(graph, graph.vertex(root), bfs);
  
  c = graph.new_edge_property('bool');
  c.a = True;
  for e in span.edges:
    c[e] = False;
  gv = gt.GraphView(graph, efilt = c); 
  
  
  return bfs;
  

def cycle_basis(graph, root):
  """Returns the minimal cycle basis of the graph from the minimal spanning tree at the root node
  
  Arguments
  ---------
    graph : Graph
      graph
    root : int
      root node from which to start
      
  Returns
  -------
    spanning_tree : array
    cycles : array
      
  Note
  ----
    Assumes graph is connected
  """

  return s,c


def cycle_basis_vectors(graph, root):
  n_edges = graph.num_edges();
  n_nodes = graph.num_vertices();
  
  s,c = cycle_basis(graph, root);
  pre = s.pre;
  dist = s.dist;
  
  cb = np.zeros((c.a.sum(), n_edges), dtype = bool);
  i = 0;
  for e in graph.edges():
    if c[e]:
      cv = cb[i];
      v1 = e.source();
      v2 = e.target();    
    
      if dist[v1] < dist[v2]:
        v1,v2 = v2,v1;
        
      while dist[v1] > dist[v2]:
        
    while 

    
    while e1.target() != root:
      cv[graph.edge_index[e]] = np.logical_xor(;
      
    s.pred[int(e.target())




def is_singular_matrix(mat):
  

A = Vectors';
% cols are vectors

Nrow=size(A,1);
Ncol=size(A,2);

if any(sum(A,2) == 0)
    isZeroRow = 1;
else
    isZeroRow = 0;
end

%check if there are same vectors with single edges (?) -> not useful for cycles ? -> skip?
AsingleEdges = A(:,sum(full(A),1)==1);
AsingleEdges = unique(AsingleEdges','rows')';
NsingleEdges = size(AsingleEdges,2);
if NsingleEdges<Nrow
    isNotSingularForSure=0;
else if NsingleEdges==Nrow
        isNotSingularForSure=1;
    else
        disp('Impossible!');
        isNotSingularForSure=2;
    end
end

%sort to upper triangular like matrix
jrow=1;
while jrow <= Nrow && isZeroRow==0 && isNotSingularForSure==0
    
    [dum i_max]  = max(abs(A(jrow:Nrow, jrow)));
    i_max = i_max+jrow-1;
    A([i_max jrow], :) = A([jrow i_max ], :);
    
    %make 'upper triangular'
    for i = jrow + 1 : Nrow
        
        if A(i, jrow)~=0
            %B  = xor(A(i, :) , A(jrow, :));
            A(i, :) = xor(A(i, :) , A(jrow, :));
        else
            %B  =A(i, :);
        end
        
        %A(i, :) = B;
        
    end
    
    if any(sum(A,2) == 0)
        isZeroRow = 1;
    end
    
    jrow = jrow+1;
end

% check if matrix is diagonal
AsingleEdges = A(:,sum(full(A),1)==1);
AsingleEdges = unique(AsingleEdges','rows')';
NsingleEdges = size(AsingleEdges,2);
if NsingleEdges<Nrow
    isNotSingularForSure=0;
else if NsingleEdges==Nrow
        isNotSingularForSure=1;
    else
        disp('Impossible!');
        isNotSingularForSure=2;
    end
end




% back substitution
jcount = 1;
while isZeroRow==0 && jcount<500000 && isNotSingularForSure==0

    % find cols that have non single entry
    KcolSingleOnes = find(sum(full(A),1)==1);
    %[dum KcolFirstOnes] = max(full(A'));
    [dum KcolFirstOnes] = max(A,[],2);
    KcolFirstOnes = full(KcolFirstOnes);
    KcolFirstOnesNotSingle = my_setdiff(KcolFirstOnes, KcolSingleOnes);
    
    if ~isempty(KcolFirstOnesNotSingle)
        columnLastOne =max(KcolFirstOnesNotSingle);
        
        rowLastOne = find(KcolFirstOnes==columnLastOne);
        
        rowsToConsid = find(A(:,columnLastOne));
        
        rowsForSubtraction = rowsToConsid(rowsToConsid~=rowLastOne(1));
        
        for js=1:length(rowsForSubtraction)
                A(rowsForSubtraction(js),:) = xor(A(rowsForSubtraction(js),:),A(rowLastOne(1),:));
        end
 
    else
        isNotSingularForSure=1;
    end
    
    isZeroRow = any(sum(full(A),2)==0);
    jcount = jcount+1;
end

if jcount>400000
    disp(['Problem: too many total iter: ' int2str(jcount)])
end


if isZeroRow==1 && isNotSingularForSure==0
    isSingular = 1;
    %disp('I think Matrix is singular')
else if (isZeroRow==1 && isNotSingularForSure==1) || (isZeroRow==0 && isNotSingularForSure==0)
        isSingular = 2;
        disp('Problem')
    else if isZeroRow==0 && isNotSingularForSure==1
            isSingular = 0;
            %disp('Not singular')
        end
    end
end

VectorsF = A';



if __name__ == "__main__":
  import graph_tool as gt
  import graph_tool.draw as gtd
  import ClearMap.Analysis.Graphs.CycleHierarchy as cyc
  reload(cyc)
  
  g = gt.Graph(directed = False);
  g.add_vertex(n = 5);
  g.add_edge_list([(0,1), (1,2), (1,3), (2,4), (4,0)])
  g.add_vertex(n=2);
  g.add_edge_list([(5,6),(0,5),(6,0)]);
  pos = gtd.sfdp_layout(g)
  gtd.graph_draw(g, pos = pos)
  
  b = cyc.spanning_tree(g, 1)
  
  efilt = g.new_edge_property('bool');
  for e in b.edges:
    efilt[e] = True;
  gv = gt.GraphView(g, efilt = efilt);
  gtd.graph_draw(gv, pos = pos);
  
  
  s,c = cyc.cycle_basis(g, 1);
  
  gv = gt.GraphView(g, efilt = c);
  gtd.graph_draw(gv, pos = pos);
  
  cyc.cycle_length(s,c);
  
  