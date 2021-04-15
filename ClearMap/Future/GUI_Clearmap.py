#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 18:22:15 2019

@author: thomas.topilko
"""
             
with open(os.path.join(clearmapDir,'Parameters.py'), 'r') as f :
                
    for n,line in enumerate(f) :
        
        print(line.split(" "))
        
        if line.split(" ") == ['detectCellShapeParameter', '=', '{\n'] :
            pos = n+1;

backgroundThreshPos = None; 
cellSizeThreshPos = None;     

with open(os.path.join(p.clearmapDir,"Temp/Parameters.py"),'w') as new_file:
                
    with open(os.path.join(p.clearmapDir,'Parameters.py')) as old_file:
        
        for n,line in enumerate(old_file) :
            
            if line.split(" ") == ['detectCellShapeParameter', '=', '{\n'] :
                
                backgroundThreshPos = n+1;
                
            elif line.split(" ") == ['thresholdPointParameter', '=', '{\n'] :
                
                cellSizeThreshPos = n+1;
            
            if n == backgroundThreshPos :
                
                new_file.write("    'threshold' : {0}, \n".format(600));
                
            elif n == cellSizeThreshPos :
                
                new_file.write("    'threshold' : {0}, \n".format((50,933)));
                
            else :
            
                new_file.write(line)