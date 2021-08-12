# %%###########################################################################
### Setting up operations to perform
##############################################################################

# !/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 17:47:32 2017

@author: thomas.topilko
"""

import os

clearmapDir = '/home/thomas.topilko/PycharmProjects/clearmap_gui_tomek'

if not os.getcwd() == clearmapDir:
    os.chdir(clearmapDir)

clearmapRessourcesDir = os.path.join(clearmapDir, "ClearMap_Ressources")
homedir = '/raid/thomas.topilko/Youenn_Psilocybin'  # Directory where the data folders will be placed '/mnt/raid/Thomas_TOPILKO/'
experiment = '210602'  # Prefix of the name of the folders to analyze. Format : "name-XXX"
# If the format is "name" make sure that the name of the folder exactly matches the variable "experiment"

import ClearMap.Tomek_Utilities as ut
import ClearMap.GUI_Clearmap as GUI
import Run as run

"""

Operations format : 'XXXXXXXXX' (len == 10) X = T or F 
T -> True : the function will be executed
F -> False : the function will not be executed
all -> will return 'TTTTTTTTTT'
none -> will return 'FFFFFFFFFF'

X1 = Stiching
X2 = Resampling auto
X3 = Resampling cfos/the data that was stiched
X4 = Align cFOS to Auto
X5 = Align Template to Auto
X6 = Align Auto to Template
X7 = Transform to Template
X8 = Cell Detection
X9 = Heatmap generation
X10 = Transform data from template alignement

"""

# %%###########################################################################
### Launch Manual/GUI mode to select future operations
##############################################################################

operations = GUI.launchGUI(homedir, experiment)

# %%###########################################################################
### Launching pipeline
##############################################################################

ClearMapParams = {
    "clearmapDir": clearmapDir,
    "clearmapRessourcesDir": clearmapRessourcesDir,
    "homedir": homedir,
    "experiment": experiment,
    "playSound": True,
    "operations": operations["operations"],
    "orient": operations["orients"],
    "bgtParams": operations["bgthresholds"],
    "cstParams": operations["csthresholds"],
    "cfosX": operations["cfosX"],
    "cfosY": operations["cfosY"],
    "cfosZ": operations["cfosZ"],
}

if os.path.exists(os.path.join(clearmapDir, 'Run.py')):
    print(ut.coloredMessage('[INFO] Run file has been successfully detected', 'darkgreen'))
else:
    raise RuntimeError('[INFO] No Run file detected in homedir!')

run.LaunchClearMap(**ClearMapParams)
