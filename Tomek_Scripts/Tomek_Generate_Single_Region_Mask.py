#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 14:12:00 2019

@author: thomas.topilko
"""

import os;
import numpy as np;

import ClearMap.IO as io;

testFile = "/mnt/vol00-renier/Catarina/190713/Catarina/190701-Patricia-01/upsampled_annotation.tif";

data = io.readData(testFile);

#mask = data != 689; #689
mask = data != 263; #689

buf = np.full_like(data, 0);

buf[~mask] = 1;

io.writeData( os.path.join(os.path.dirname(testFile), "AVP.tif"), buf )
