# -*- coding: utf-8 -*-
"""
Created on Fri Sep  1 11:34:09 2017

@author: ckirst
"""

import numpy as np
from scipy import ndimage as ndi

def hysteresis_threshold(image, low, high):
    """Apply hysteresis thresholding to `image`.
    This algorithm finds regions where `image` is greater than `high`
    OR (`image` is greater than `low` *and* that region is connected to
    a region greater than `high`.
    Parameters
    ----------
    image : array, shape (M,[ N, ..., P])
        Grayscale input image.
    low : float, or array of same shape as `image`
        Lower threshold.
    high : float, or array of same shape as `image`
        Higher threshold.
    Returns
    -------
    thresholded : array of bool, same shape as `image`
        Array in which `True` indicates the locations where `image`
        was above the hysteresis threshold.
    Examples
    --------
    >>> image = np.array([1, 2, 3, 2, 1, 2, 1, 3, 2])
    >>> apply_hysteresis_threshold(image, 1.5, 2.5).astype(int)
    array([0, 1, 1, 1, 0, 0, 0, 1, 1])
    References
    ----------
    .. [1] J. Canny. A computational approach to edge detection.
           IEEE Transactions on Pattern Analysis and Machine Intelligence.
           1986; vol. 8, pp.679-698.
           DOI: 10.1109/TPAMI.1986.4767851
    """
    low = np.clip(low, a_min=None, a_max=high)  # ensure low always below high
    mask_low = image > low
    mask_high = image > high
    # Connected components of mask_low
    labels_low, num_labels = ndi.label(mask_low)
    # Check which connected components contain pixels from mask_high
    sums = ndi.sum(mask_high, labels_low, np.arange(num_labels + 1))
    connected_to_high = sums > 0
    thresholded = connected_to_high[labels_low]
    
    print('done')
    return thresholded
 

import ClearMap.ImageProcessing.Filter.StructureElement as sel   
    
def local_max_hysteresis_thresholding(image, size, high, height, low, percentile):
    
    #find local maxima
    se = sel.structureElement('Disk', size);
    data_max = ndi.filters.maximum_filter(image, footprint = se);
    data_max[data_max != image] = 0;
    #data_max[np.logical_not(mask_max)] = 0;
    
    #data_min = ndi.filters.minimum_filter(image, footprint = se);
    #mask_max = ((data_max - data_min) > height)
    #mask_max = np.logical_and(mask_max, image > high);    
    mask_max = data_max > high;    
    
    res = image.copy();
    res[mask_max == 0] = 0
    
    labels_max, num_labels = ndi.label(mask_max, structure = np.ones((3,3), dtype = bool))


    maxs = ndi.median(image, labels_max, np.arange(num_labels + 1));
    
    #this is slow for testing
    thresholded = np.zeros(image.shape, dtype = bool);
    
    data_min = image > low    
    for l in range(1, num_labels+1):
        print("%d / %d" % (l, num_labels))
        mask_low = np.logical_and(image > (maxs[l] * percentile), data_min);
        labels_low, num_labels2 = ndi.label(mask_low)
        
        # Check which connected components contain pixels from mask_high
        sums = ndi.sum(labels_max == l, labels_low, np.arange(num_labels2 + 1))
        connected_to_high = sums > 0
        thresholded = np.logical_or(thresholded, connected_to_high[labels_low]);
    
    print('done')
    return thresholded
    
if __name__ == "__main__":
    from init_tests import *;
    import ClearMap.ImageProcessing.Filter.Thresholding as ths
    reload(ths);
    
    image = img[:,:,300];    
    res = ths.local_max_hysteresis_thresholding(image, size = 10, high = 1000, height = 10, low = 500, percentile = 0.7)
    
  
    plt.figure(4, figsize = (32,7)); plt.clf();
    ax = plt.subplot(1,2,1);
    pim1 = plt.imshow(image);
    plt.subplot(1,2,2, sharex = ax, sharey = ax);
    pim2 = plt.imshow(res);
    plt.tight_layout();
    