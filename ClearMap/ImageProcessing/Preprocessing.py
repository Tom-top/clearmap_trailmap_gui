# -*- coding: utf-8 -*-
"""
Preprocessing Pipelines
"""
__author__    = 'Christoph Kirst <ckirst@rockefeller.edu>'
__license__   = 'MIT License <http://www.opensource.org/licenses/mit-license.php>'
__copyright__ = 'Copyright 2017 by Christoph Kirst, The Rockefeller University, New York City'

import numpy as np;

from scipy import ndimage as ndi  

import skimage.filters.rank as rnk;

import ClearMap.ImageProcessing.Filter.Rank as rnk3d

from ClearMap.ImageProcessing.Filter.StructureElement import structureElement;


def percentile(image, selem, out = None, p0 = 0, max_bin = None,
               out_dtype = None, mask = None, shift_x = False, shift_y = False):
    """Percentile rank filer with normalized max_bin parameter"""
  
    rnk.generic.assert_nD(image, 2);
    
    image, selem, out, mask, m_bin = rnk.generic._handle_input(image, selem, out, mask, out_dtype);
    
    if max_bin is None:
      max_bin = m_bin; 
    #print max_bin
    
    rnk.percentile_cy._percentile(image, selem, out=out, mask=mask, shift_x=shift_x, shift_y=shift_y, max_bin=max_bin, p0=p0, p1=0.);

    return out.reshape(out.shape[:2])


def autolevel(image, selem, out = None, p0 = 0, p1 = 1, max_bin = None,
               out_dtype = None, mask = None, shift_x = False, shift_y = False):
    """Autolevel rank filer with normalized max_bin parameter"""
  
    rnk.generic.assert_nD(image, 2);
    
    image, selem, out, mask, m_bin = rnk.generic._handle_input(image, selem, out, mask, out_dtype);
    
    if max_bin is None:
      max_bin = m_bin;
    
    return rnk.percentile_cy._autolevel(image, selem, out=out, mask=mask, shift_x=shift_x, shift_y=shift_y, max_bin=max_bin, p0=p0, p1=p1);


def threshold_local(image, selem, out = None, p0 = 0, p1 = 1, max_bin = None,
                    out_dtype = None, mask = None, shift_x = False, shift_y = False):
    """Autolevel rank filer with normalized max_bin parameter"""
  
    rnk.generic.assert_nD(image, 2);
    
    image, selem, out, mask, m_bin = rnk.generic._handle_input(image, selem, out, mask, out_dtype);
    
    if max_bin is None:
      max_bin = m_bin;
    
    return rnk.percentile_cy._threshold(image, selem, out=out, mask=mask, shift_x=shift_x, shift_y=shift_y, max_bin=max_bin, p0=p0, p1=0);




def preprocess_old(source, sink = None, clipmin = 700, clipmax = 20000, norm = 255, nbins = 255,
                                    backgroundSize = 100, backgroundPercentile = 0.2,
                                    autolevelSize = 100, autolevelPercentile0 = 0.4, autolevelPercentile1 = 1.0,
                                    thresholdSize = 100, thresholdPercentile = 0.75 ):
  
  if sink is None:
    sink = source.copy(order = 'A');
  else:
    sink[:] = source;
    
  if sink.flags['F_CONTIGUOUS']:
    sinkTranspose = [2,1,0];
    sink = sink.transpose([2,1,0]);
  else:
    sinkTranspose = None;
    
  if not sink.flags['C_CONTIGUOUS']:
    print('warning the array is not in contiguous form, processing might be slow !');
    print(sink.flags);
  
  # clip data and normalize
  sink[sink > clipmax] = clipmax;
  sink[sink < clipmin] = clipmin;
  sink[:] = norm * (np.array(sink, dtype = float)-clipmin)/(clipmax-clipmin);   
  
  # remove background via percentiled bottomhat
  tmp = sink.copy(order = 'A');  
  struct = structureElement(setype = 'Disk', sesize = backgroundSize);
  for z in range(sink.shape[0]):
    tmp2 = sink[z].copy(order = 'A');
    res = percentile(tmp2, selem = struct, out = tmp[z], p0 = backgroundPercentile, max_bin = nbins);
  tmp = np.array(sink, dtype = int, order = 'A') - tmp;
  tmp[tmp < 0] = 0;
  sink[:] = tmp;

  # autolevel
  struct = structureElement(setype = 'Disk', sesize = autolevelSize);
  for z in range(sink.shape[0]):
    tmp = sink[z].copy();
    res = autolevel(tmp, selem = struct, out = sink[z], p0 = autolevelPercentile0, p1= autolevelPercentile1, max_bin = nbins);

  # threshold
  struct = structureElement(setype = 'Disk', sesize = thresholdSize);
  for z in range(sink.shape[0]):
    tmp = sink[z].copy();
    res = threshold_local(tmp, selem = struct, out = sink[z], p0 = thresholdPercentile, max_bin = nbins);
  
  if sinkTranspose is not None:
    sink = sink.transpose(sinkTranspose);

  sink[source < clipmin] = 0;
  
  return sink;
  
  
import scipy.ndimage.morphology as morph
  
def cleanBinary(source, sink = None, erosionSize = 3, dilationSize = 5, fill = True):
  
  # morphological erosion -> get rid of small individual pizel
  if erosionSize is not None:
    tmp = np.empty_like(source, order = 'A');
    struct = structureElement(setype = 'Disk', sesize = [erosionSize]*3);
    morph.binary_erosion(source, output = tmp, structure = struct);
  else:
    tmp = source;
  
  # morphological dilation -> fill gaps
  if sink is None:
    sink = np.empty_like(source, order = 'A');
  
  if dilationSize is not None:
    struct = structureElement(setype = 'Disk', sesize = [dilationSize]*3);
    morph.binary_dilation(tmp, output = sink, structure = struct);
  else:
    sink[:] = tmp;
  tmp = [];
    
  #fill holes
  if fill:
    morph.binary_fill_holes(sink, output = sink);
  
  return sink;
      


import ClearMap.ImageProcessing.Filter.Clip.Clip as clp

def preprocess(source, sink = None, clipmin = 500, clipmax = 3000, norm = 255, median = 4):
  
  if sink is None:
    sink_dtype = 'uint16';
  else:
    sink_dtype = sink.dtype;
  
  sink = clp.clip(source, sink = sink, clipmin =clipmin, clipmax = clipmax, norm = norm, sink_dtype = sink_dtype);
  
  # clip data and normalize
  #sink[sink > clipmax] = clipmax;
  #sink[sink < clipmin] = clipmin;
  #sink[:] = norm * (np.array(sink, dtype = float)-clipmin)/(clipmax-clipmin);   
  
  print('done normalizing');
  
  
  if sink.flags['F_CONTIGUOUS']:
    sinkTranspose = [2,1,0];
    sink = sink.transpose([2,1,0]);
  else:
    sinkTranspose = None;
  
  #3d median filtering
  sink[:] = rnk3d.median(sink, structureElement('Disk', (median,)*3));
  
  print('done median');
  
  if sinkTranspose is not None:
    sink = sink.transpose(sinkTranspose);
  
  return sink;



import ClearMap.ImageProcessing.Filter.Curvature as cur

def binarize(source, sink = None, threshold_main = 10, threshold_tube = [0.27, 5], sigma_tube = [2.5, 1.0], fill = True):
  
  if sink is None:
    sink = np.zeros(source.shape, dtype = bool, order = 'F');
  sink[:] = (source > threshold_main);
  
  print('done main thresholding')
  
  #if sink.flags['F_CONTIGUOUS']:
  #  sinkTranspose = [2,1,0];
  #  sink = sink.transpose([2,1,0]);
  #else:
  #  sinkTranspose = None;
  
  for t,s in zip(threshold_tube, sigma_tube):
    smooth = ndi.gaussian_filter(np.asarray(source, dtype = 'float32'), sigma = s);
    print('done gaussian sigma = %f' % s);
    sink[:] = np.logical_or(sink, cur.tubeness(smooth, threshold = t));
    print('done tubness thresholding sigma=%f threshold=%r' % (s,t));
    
  #if sinkTranspose is not None:
  #  sink = sink.transpose(sinkTranspose);
    
  #fill holes
  #if fill:
    #print sink.shape
    #print sink.flags;
  #  morph.binary_fill_holes(sink, output = sink);
    #print sink.flags;
  #  print('done filling');
    
  if fill:
    for z in range(sink.shape[2]):
      morph.binary_fill_holes(sink[:,:,z], output = sink[:,:,z]);
    print('done filling');
  
  return sink;



def binarize_arteries(source, sink = None, threshold = 10, fill = False, dilate = None, opening = 3):
  
  if sink is None:
    sink = np.zeros(source.shape, dtype = bool, order = 'F');
  sink[:] = (source > threshold);
  print('done main thresholding')
  print(sink.sum())
  
  #if sink.flags['F_CONTIGUOUS']:
  #  sinkTranspose = [2,1,0];
  #  sink = sink.transpose([2,1,0]);
  #else:
  #  sinkTranspose = None;
  
  if fill:
    for z in range(sink.shape[2]):
      morph.binary_fill_holes(sink[:,:,z], output = sink[:,:,z]);
    print('done filling');
    
  if opening is not None:
    se = structureElement('Disk', opening, ndim = 3).astype('uint8');
    #sink[:] = rnk3d.maximum(rnk3d.minimum(sink, se), se);
    sink[:] = ndi.morphology.binary_opening(sink, se);
    print('done opening');
    print(sink.sum())
    
  if dilate is not None:
    se = structureElement('Disk', dilate, ndim = 3).astype('uint8');
    sink[:] = rnk3d.maximum(sink, se);
    print('done main dilating')
    print(sink.sum())
  
  return sink;





def fill(source, sink = None, axis = [2]):
  
  if sink is None:
    sink = np.zeros(source.shape, dtype = bool, order = 'F');
  
  if axis is None:
    morph.binary_fill_holes(source, output = sink);
  else:
    #fill holes
    if axis is all:
      axis = [0,1,2];
    
    sp = source;
    if 2 in axis:
      for z in range(sink.shape[2]):
        sink[:,:,z] = morph.binary_fill_holes(sp[:,:,z]);
      sp = sink;
      print('done filling z');
    
    if 1 in axis:
      for y in range(sink.shape[1]):
        sink[:,y,:] = morph.binary_fill_holes(sp[:,y,:]);
      sp = sink;
      print('done filling y');        
        
    if 0 in axis:
      for x in range(sink.shape[0]):
        sink[x,:,:] = morph.binary_fill_holes(sp[x,:,:]);
      print('done filling x');
  
  print('done filling');

  return sink;



import ClearMap.ParallelProcessing.SharedMemoryProcessing as smm 

def multiFill(source, sink, axis = [2,1,0], **parameter):
  for a in axis:
    smm.process(source, sink, function = fill, chunkAxis = a, axis = [a], **parameter)
    source = sink;



import sys

from ClearMap.Utils.Timer import Timer;

import cv2

def removeBackgroundCV2(source, sink = None, size = None, verbose = False, out = sys.stdout):
    """Remove background via subtracting a morphological opening from the original image"""
    
    out.write('Background routine: size = %r' % (size,));
    
    if size is None:    
        return source;
    
    # change type to float in order to prevent 
    img = np.array(source, dtype = float);
    
    out.write('Background source size = %r' % (img.shape,));

    timer = Timer();
    # background subtraction in each slice
    se = structureElement('Disk', size).astype('uint8');
    for z in range(source.shape[2]):
         #img[:,:,z] = img[:,:,z] - grey_opening(img[:,:,z], structure = structureElement('Disk', (30,30)));
         #img[:,:,z] = img[:,:,z] - morph.grey_opening(img[:,:,z], structure = self.structureELement('Disk', (150,150)));
         img[:,:,z] = img[:,:,z] - cv2.morphologyEx(img[:,:,z], cv2.MORPH_OPEN, se)
         #img[:,:,z] = img[:,:,z] - cv2.morphologyEx(img[:,:,z], cv2.MORPH_CLOSE, se)
    
    img[img < 0] = 0;
    img = np.array(img, dtype = source.dtype);
    

    #if verbose:
    out.write(timer.elapsedTime(head = 'Background') + '\n');
    
    return img


def preprocessBC(source, sink = None, 
                 clipmin = 300, clipmax = 9000, norm = 255,
                 keep_above = 5000,
                 backgroundSize = (50,50,1), backgroundPercentile = 0.8,
                 median = 4,
                 verbose = False, out = sys.stdout):
  """Remove background via subtracting a morphological opening from the original image"""
  
  
  #print 'function...';
  #print source.flags
  
  timer = Timer();
  
  if sink is None:
    sink_dtype = 'uint16';
  else:
    sink_dtype = sink.dtype;
  
  sink = clp.clip(source, sink = sink, clipmin =clipmin, clipmax = clipmax, norm = norm, sink_dtype = sink_dtype);
    
  out.write(timer.elapsedTime(head = 'Clipping') + '\n');
  
  
  #3d median filtering
  if median is not None:
    
    #TODO: fix fortran vs c order in rank filter !
    #if source.flags['F_CONTIGUOUS']:
    #  source = source.transpose([2,1,0]);
    #  transpose = [2,1,0];
    #else:
    # transpose = None;
      
    sink[:] = rnk3d.median(sink, structureElement('Disk', (median,)*3));
    
    #if transpose is not None:
    #  sink = sink.transpose(transpose);
    #  source = source.transpose(transpose);
      
    
    out.write(timer.elapsedTime(head = 'Median') + '\n');
  
  
  # remove background via percentiled bottomhat 
  if backgroundSize is not None:
#    struct = structureElement(setype = 'Disk', sesize = backgroundSize); 
#    tmpout = sink[0].copy(order = 'A');
#    for z in range(sink.shape[0]):
#      tmp = sink[z].copy(order = 'A');
#      percentile(tmp, selem = struct, out = tmpout, p0 = backgroundPercentile, max_bin = norm);
#      tmp = np.array(tmp, dtype = int, order = 'A') - tmpout;
#      tmp[tmp < 0] = 0;
#      sink[z] = tmp;

      #TODO: see if we need to transpose / or not ?
      struct = structureElement('Disk', backgroundSize);
      tmp = rnk3d.percentile(sink, struct, p0 = backgroundPercentile);
      tmp = np.array(sink, dtype = int, order = 'A') - tmp;
      tmp[tmp < 0] = 0;
      tmp[sink >= norm] = norm;
      sink[:] = tmp;
      tmp = [];
      
      out.write(timer.elapsedTime(head = 'Background correction') + '\n');
  
  
  if keep_above is not None:
    #tmp = sink.copy();
    #tmp[source > keep_above] = norm;
    #sink[:] = tmp;
    #tmp = [];
    sink[source > keep_above] = norm;
    
    out.write(timer.elapsedTime(head = 'Restoring high intensity pixels') + '\n');
    
  
  return sink;

