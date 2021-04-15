#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 27 14:16:03 2020

@author: thomas.topilko
"""

import os;

folder = "/raid/200217-700/200515_cfos_09-40-38";

for f in os.listdir(folder) :
    
    path2File = os.path.join(folder, f);
    new = [f.split("_")[0]] + ["cfos"] + f.split("_")[2:];
    newPath2File = os.path.join(folder, "_".join(new));
    
    print(new)
    
    os.rename(path2File, newPath2File);