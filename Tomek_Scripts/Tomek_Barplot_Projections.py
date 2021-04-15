#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 17:32:50 2019

@author: thomas.topilko
"""

import os;
import matplotlib.pyplot as plt;
import numpy as np;
#import seaborn as sns;
import pandas as pd;

workingDir = "/raid/Thomas_TOPILKO";
regionIdPath = "/home/thomas.topilko/Documents/ClearMap_GUI/ClearMap_Ressources/Regions_annotations/regions IDs.csv";
voxelNumbersPath = "/home/thomas.topilko/Documents/ClearMap_GUI/Tomek_Scripts/Voxel_Numbers_Annotation.xlsx";
projectionTablePath = os.path.join(workingDir, "counts-cells_table_Pregnant_vs_Virgin.csv");
matedActivityTablePath = os.path.join(workingDir, "counts-cells_table_mated_cotton_vs_mated_cup.csv");
virginActivityTablePath = os.path.join(workingDir, "counts-cells_table_virgin_cotton_vs_virgin_cup.csv");

IdsMetaRegions = ["OLF", "CTXsp", "TH", "HY", "MB", "MY", "CB"];

IdTableDF = pd.read_csv(regionIdPath, header=0, usecols=["id", "red", "green", "blue", "acronym", "structure_order", "parent_id"]);
print("{0} Regions detected in voxel table".format(len(IdTableDF["id"])))

voxelTableDF = pd.read_excel(voxelNumbersPath, header=0, usecols=["id", "voxel"]);
print("{0} Regions detected in voxel table".format(len(voxelTableDF["id"])));

def acronymToID(name, table) :
    
    pos = np.where( table["acronym"] == name );
    Id = table["id"][pos[0]];
    
    return Id;

#testName = "AN";
#np.where( voxelTableDF["id"] == acronymToID(testName, IdTableDF).tolist()[0] )

#nColumns = pd.read_csv(projectionTablePath, header=0,nrows=0).columns;
projectionTableDF = pd.read_csv(projectionTablePath, header=0, usecols=["id", " mean1", " mean2"]);

orderedIds = set(projectionTableDF["id"]);

#pregnantCotton = ["6", "13", "14", "24", "25", "28"];
matedActivityTableDF = pd.read_csv(matedActivityTablePath, header=0, usecols=["id", " mean1", " mean2"]);
virginActivityTableDF = pd.read_csv(virginActivityTablePath, header=0, usecols=["id", " mean1", " mean2"]);

sink = plotMatrix2 = [];
testedGroup = matedActivityTableDF;

for region in orderedIds :
    
    posRegionInActivityTable = testedGroup["id"][testedGroup["id"] == region].index.tolist();
#    posRegionInActivityTable2 = virginActivityTableDF["id"][virginActivityTableDF["id"] == region].index.tolist();
    activityInRegion0 = testedGroup[" mean1"].loc[posRegionInActivityTable].values.tolist();
    activityInRegion1 = testedGroup[" mean2"].loc[posRegionInActivityTable].values.tolist();
    
    posRegionInDensityTable = projectionTableDF["id"][projectionTableDF["id"] == region].index.tolist();
    densityInRegion = projectionTableDF[" mean1"].loc[posRegionInDensityTable].values.tolist();
    
    posRegionInVoxelTable = voxelTableDF["id"][voxelTableDF["id"] == region].index.tolist();
    voxelsInRegion = voxelTableDF["voxel"].loc[posRegionInVoxelTable].values.tolist();
    
    posRegionInIdTable = IdTableDF["id"][IdTableDF["id"] == region].index.tolist();
    redInRegion = IdTableDF["red"].loc[posRegionInIdTable].values.tolist();
    greenInRegion = IdTableDF["green"].loc[posRegionInIdTable].values.tolist();
    blueInRegion = IdTableDF["blue"].loc[posRegionInIdTable].values.tolist();
    acronymInRegion = IdTableDF["acronym"].loc[posRegionInIdTable].values.tolist();
    structureOrderInRegion = IdTableDF["structure_order"].loc[posRegionInIdTable].values.tolist();
    parentIdInRegion = IdTableDF["parent_id"].loc[posRegionInIdTable].values.tolist();
    
#    print("\n")
#    print(region, activityInRegion, densityInRegion, redInRegion, greenInRegion, blueInRegion, acronymInRegion)
    availableData = [] not in (activityInRegion0, activityInRegion1, densityInRegion, redInRegion,\
                                greenInRegion, blueInRegion, acronymInRegion, voxelsInRegion,\
                                structureOrderInRegion, parentIdInRegion);
    
    if availableData and acronymInRegion[0] not in IdsMetaRegions :
        
        activityInRegion = [activityInRegion0[0] - activityInRegion1[0]];
#        print(activityInRegion0 , activityInRegion1)
        activityInRegionNorm = [activityInRegion[0]/voxelsInRegion[0]];
        densityInRegionNorm = [densityInRegion[0]/voxelsInRegion[0]];
        
        if redInRegion[0] != greenInRegion[0] != blueInRegion[0]  and\
        (redInRegion[0], greenInRegion[0], blueInRegion[0]) not in [(176,255,184)] :
        
            sink.append(np.array([region, densityInRegionNorm[0], activityInRegionNorm[0], redInRegion[0],\
                               greenInRegion[0], blueInRegion[0], acronymInRegion[0], structureOrderInRegion[0],\
                               parentIdInRegion[0], voxelsInRegion[0], densityInRegion[0]]));
    
    
regionsStructureOrder = np.array(np.array(sink)[:,7], dtype="int64");
ordered = np.argsort(regionsStructureOrder);
ordered = np.argsort(ordered);

def orderList(source, mask) :
    
    sink = np.full_like( source, np.full_like(source[0], 0.) );
    
    for n, x in zip(mask, source) :
        
        sink[n] = x;
        
    return sink;

colors = [(r,g,b) for r,g,b in zip(np.array(np.array(sink)[:, 3], dtype="float64")/255.,\
           np.array(np.array(sink)[:, 4], dtype="float64")/255., np.array(np.array(sink)[:, 5], dtype="float64")/255.)];

orderedDensity = orderList(np.array(np.array(sink)[:, 1], dtype="float64"), ordered);
orderedDensityRaw = orderList(np.array(np.array(sink)[:, 10], dtype="float64"), ordered);
totalDensity = np.sum(orderedDensityRaw);
orderedColors = orderList(colors, ordered);
orderedNames = orderList(np.array(np.array(sink)[:, 6], dtype="str"), ordered);

fig = plt.figure(figsize=(10, 7));
ax = plt.subplot(111);
#threshYHigh = 0.2;
threshYHigh = 10000;
#threshYLow = 0.045;
threshYLow = 2000;
offset = 0.003;

ax.bar(np.arange(0, len(orderedDensityRaw)),\
       orderedDensityRaw,\
       color=orderedColors,\
       edgecolor="black",\
       linewidth=0.4,\
       );

for x,y,n in zip(np.arange(0, len(orderedDensityRaw)), orderedDensityRaw, orderedNames) :
    
    if y > threshYLow and y < threshYHigh : 
        
        ax.text(x, y+offset, n, fontsize=15., rotation=90., ha="center", va="bottom", weight='bold');
        
    elif y > threshYHigh :
        
        ax.text(x, threshYHigh+offset, n+" : {0}".format(round(y,2)), fontsize=15., rotation=90., ha="center", va="bottom", weight='bold');
        
#    if n == "LPO" :
#        
#        ax.text(x, y+0.005, n, fontsize=6., rotation=90., ha="center", va="center");

ax.set_ylim(0,threshY);
ax.set_xlim(-0.5,len(YNewNorm)+0.5);
ax.set_ylabel("Density of detected puncta per brain area (n points / n Pixels in area)", fontsize=15.);

labels = ["Cerebral cortex", "Striatum", "Pallidum", "Thalamus", "Hypothalamus", "Midbrain", "Hindbrain", "Cerebellum"];
handlesColors = ["#70FF70", "#98D6F9", "#8599CC", "#FF7080", "#E64438", "#FF64FF", "#FF9B88", "#F0F080"];
handles =  [plt.Rectangle((0, 0), 1, 1, fc=c, edgecolor="black") for c in handlesColors]
plt.legend(handles=handles, labels=labels, shadow=True);
ax.set_xticks([]);
plt.savefig( os.path.join("/home/thomas.topilko/Desktop", "Density_Projections_Per_Area.png"), dpi=200. )
        
        
        
        
        
        
        
        
        
        
###############################################################################      
#Checking missing ids detected from annotation (parent structures)
###############################################################################    
        
parentIds = fcr["parent_id"];
rawIds = fcr["id"];
added = 0;

matchingIDS = [];

for id1 in idsVoxels.keys() :
    
    if id1 in ids.tolist() :
        
        matchingIDS.append(id1)

while len(matchingIDS)+added != len(ids) :
    
    print("\n");

    for ID in ids :
        
        try :
            
            idsVoxels[ID];
        
        except :
            
            parents = [];
            
            for rawId, parent in zip(rawIds, parentIds) :
                
                if parent == ID :
                    
                    parents.append(rawId);
                    
            allOk = True;
            
            print(parents)
                    
            for parent in parents :
                
                try :
                    
                    idsVoxels[parent];
                    
                except :
                    
                    allOk = False;
                    
            if allOk :
                
                idsVoxels[ID] = 0;
                
                added += 1;
                
                print("ID n_{0} added".format(ID));
                print(len(matchingIDS)+added, len(ids))
                
                for parent in parents :
                
                    idsVoxels[ID] += idsVoxels[parent];
                    
YNormalized = [];
idsNormalized = [];

for ID, y in zip(idsFiltered,Y) :
    
    try :
    
        YNormalized.append( y/idsVoxels[ID] );
        idsNormalized.append( ID );
    
    except :
        
        pass;




















































