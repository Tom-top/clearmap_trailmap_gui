#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 10:34:04 2019

@author: thomas.topilko
"""

#%%###########################################################################
### Analysis
##############################################################################

import os;

clearmapDir = '/home/thomas.topilko/Documents/clearmap_gui_tomek/';

if not os.getcwd() == clearmapDir :
    os.chdir(clearmapDir)

import ClearMap.IO.IO as io;
import ClearMap.Alignment.Resampling as rs;
import numpy as np;
from natsort import natsorted;

import ClearMap.Tomek_Utilities as ut;
import ClearMap.GUI_Analysis as GUI_Analysis;
import Parameters;

clearmapRessourcesDir = os.path.join(clearmapDir,"ClearMap_Ressources");
homedir = '/raid/Tomek'; #Directory where the data folders will be placed '/mnt/raid/Thomas_TOPILKO/'
experiment = '210225'; #Prefix of the name of the folders to analyze. Format : "name-XXX"
Parameters.analysisDir = os.path.join(homedir, "Analysis");
#Parameters.analysisDir = None;

group1 = 'SS';
group2 = 'SN';
group3 = 'NS';
group4 = 'NN';

allGroups = [group1,group2,group3,group4];

'''

Operations format : 'XXXXXXX' (len == 7) X = T or F 
T -> True : the animal will be included in the group
F -> False : the animal will not be included in the group
all -> will return 'TTTTTTT'
none -> will return 'FFFFFFF'

X1 = Stiching
X2 = Resampling auto
X3 = Resampling cfos
X4 = Align cFOS to Auto
X5 = Align Template to Auto
X6 = Cell Detection
X7 = Heatmap generation

'''

sep = "-";
folders = [];

for folder in natsorted(os.listdir(homedir)) :
    if os.path.isdir(os.path.join(homedir,folder)) :
        if folder.split(sep)[0] == experiment :
            folders.append(folder);
      
_,operations = GUI_Analysis.main(folders, experiment, allGroups);

#%%###########################################################################
### Lauching Analysis Pipeline
##############################################################################

print(ut.coloredMessage(ut.titleMessage('Starting analysis pipeline'),'darkgreen'));

sampleNames = ut.getSampleNames(homedir, experiment, sep);

AnnotationFile = os.path.join(clearmapRessourcesDir+"/Regions_annotations/","annotation_279_new.tif");    

if Parameters.analysisDir != None : 
    saveDir = os.path.join(homedir, Parameters.analysisDir);
else :
    saveDir = homedir;

nGroups = len(allGroups);
#nCombinations = np.math.factorial(nGroups) / ( np.math.factorial(2) * np.math.factorial(nGroups-2) )
combinations = ut.combinlisterep(np.arange(0,nGroups), 2);

groups = {};

for group, operation in zip(operations.keys(), operations.values()) :
    
    args = {"group" : group,
        "operation" : operation,
        "homedir" : homedir,
        "savingDir" : saveDir,
        "analysisDir" : Parameters.analysisDir,
        "experiment" : experiment,
        "sep" : sep,
        "sampleNames" : sampleNames,
        }
  
    groups[group] = ut.Group(**args);
    
    io.writeData(os.path.join(saveDir, str(group)+'.raw'), rs.sagittalToCoronalData(groups[group].mean));
  
    if Parameters.analysisDir == None : 
      
        groups[group].paths2 = [homedir+os.sep+experiment+sep+str(y)+os.sep+\
                                'cells_transformed_to_Atlas.npy' for x,y in zip(operation,sampleNames) if x == 'T'];

    else :

        groups[group].paths2 = [os.path.join(saveDir, 'cells_transformed_to_Atlas_{0}_{1}.npy'.format(experiment, y))\
                                   for x,y in zip(operation,sampleNames) if x == 'T'];
              
    groups[group].i = [fn.replace('cells_transformed_to_Atlas', 'intensities') for fn in groups[group].paths2];

#%%
        
ut.launchPairAnalysis(groups, combinations, saveDir, AnnotationFile, cutoff=0.05);
    