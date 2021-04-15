# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 15:19:22 2017

@author: ckirst
"""

import numpy as np;
data = np.load('../test_bin.npy');
data[[0,-1],:,:] = 0;
data[:,[0,-1],:] = 0;
data[:,:,[0,-1]] = 0;


import skeleton.thinVolume as s;

skel = s.get_thinned(data);

#plot the result
import os;
cmpath = '/home/ckirst/Science/Projects/WholeBrainClearing/Vasculature/Analysis/ClearVessel'
skpath = '/home/ckirst/Science/Projects/WholeBrainClearing/Vasculature/Analysis/ClearVessel/ClearMap/ImageProcessing/Skeletonization'
os.chdir(cmpath);
import ClearMap.GUI.DataViewer as dv;
dv.DataViewer([data, skel])
os.chdir(skpath)