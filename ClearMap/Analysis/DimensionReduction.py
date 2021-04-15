# -*- coding: utf-8 -*-
"""
Module hosting dimension reduction and embedding routines
"""

import matplotlib.pyplot as plt;

import sklearn.manifold as sl;


def tsne(data, n_components = 2, precomputed = False, 
               perplexity=30.0, random_state = 0, init = 'pca',
               early_exaggeration=12.0, learning_rate=200.0, n_iter=1000, **kwargs):
  """Perform t-SNE"""
  
  if precomputed:
    metric = 'precomputed';
  else:
    metric = None;
  
  tsne = sl.TSNE(n_components = n_components, init = init, random_state = random_state, 
                 metric = metric, perplexity = perplexity, early_exaggeration = early_exaggeration,
                 learning_rate = learning_rate, n_iter = n_iter, **kwargs);
  return tsne.fit_transform(data);


#TODO: for million of points -> use vispy
def plot_tsne(data, cmap = plt.cm.Spectral, colors = None, title = 't-SNE'):
  if data.ndim == 1:
    n_components = 1;
  else:
    n_components = data.shape[1];
    
  if colors is None:
    colors = range(len(data[:,0]));
  
  if n_components == 1:
    plt.plot(data);
  elif n_components == 2:
    plt.scatter(data[:,0], data[:,1], c = colors, cmap = cmap);
  else:
    ax = plt.gcf().add_subplot(1, 1, 1, projection = '3d');
    ax.scatter(data[:, 0], data[:, 1], data[:,2], c = colors, cmap = cmap)
  plt.title(title)
  plt.tight_layout();
