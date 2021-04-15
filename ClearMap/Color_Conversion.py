#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May  2 09:59:01 2019

@author: thomas.topilko
"""

import os
import matplotlib.pyplot as plt;
import numpy as np;
from multiprocessing import Process;
import ClearMap.IO as io;

path2Folder = "/network/lustre/iss01/renier/Thomas/thomas.topilko/Experiments/Nesting_Project/Tracing/Anterograde_Tracing_EWcp/190305-04_EW_Tracing_Pregnant_Virgin_215-216-217-226-247-248/LowRes/190225-215"

source = io.readData(os.path.join(path2Folder,"transform_cfos_to_template/result.tiff"));
transposed = source.transpose();

imageBuffer = np.zeros((transposed.shape[0],transposed.shape[1],transposed.shape[2],3),dtype="uint16");

def TransformColor(index,plane) :
        
    print("Reading plane {0}/{1}".format(index,transposed.shape[2]))
    
    planeBuffer = np.zeros((plane.shape[0],plane.shape[1],3),dtype="uint16");
    
    for x,line in enumerate(plane) :
        
        for y,pixel in enumerate(line) :
            
            planeBuffer[x][y] = np.array([0,int(pixel),int(pixel)]) #cyan
            
    imageBuffer[index] = planeBuffer;
        
procs = [];
 
for index ,plane in enumerate(transposed) :
    
    proc = Process(target=TransformColor, args=(index,plane,));
    procs.append(proc);
    proc.start();
 
for proc in procs:
    proc.join();

io.writeData(os.path.join(path2Folder,"Colored_Tracing.tiff"),imageBuffer)
