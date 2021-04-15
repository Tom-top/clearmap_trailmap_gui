#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 13:10:44 2019

@author: thomas.topilko
"""

import cv2;
import os;
import ClearMap.IO.IO as io;
import numpy as np;
import matplotlib.pyplot as plt;

#filePath = "/home/thomas.topilko/Documents/ClearMap_GUI/ClearMap_Ressources/Regions_annotations/annotation_full_new_EW_symetric.tif"
filePath = "/raid/Thomas_TOPILKO/190225-665/elastix_annotation_transformed_to_raw/annotation_transformed_template_to_cfos_2.tif"

#annotation = io.read(filePath);
annotation = io.readData(filePath);
#annotationCoronal = np.transpose(annotation, (0, 3, 2, 1));
#annotationCoronal = np.transpose(annotation, (1, 2, 0));
#plt.imshow(annotationCoronal[100]);
#annotationCoronal.shape

sink = [];

for i in annotationCoronal :
    
    i = np.uint8(i);
    edges = cv2.Canny(i,0,1);
    sink.append(edges);
    
sink = np.array(sink);
#io.write(os.path.join("/home/thomas.topilko/Desktop","edges2.tif"), sink)
io.writeData(os.path.join("/home/thomas.topilko/Desktop","edges2.tif"), sink)

#import matplotlib.pyplot as plt
#
#plt.subplot(111),plt.imshow(edges,cmap = 'gray');
