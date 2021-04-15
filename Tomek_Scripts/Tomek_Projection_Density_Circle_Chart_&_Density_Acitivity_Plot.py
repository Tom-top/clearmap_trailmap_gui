#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 13:47:46 2020

@author: thomas.topilko
"""

import os;
import numpy as np;
import pandas as pd;
import matplotlib.pyplot as plt;
from matplotlib.lines import Line2D;
import matplotlib;
import matplotlib.cm as cm;
from matplotlib.patches import Rectangle;
from matplotlib.colors import ListedColormap, LinearSegmentedColormap;
from mpl_toolkits.axes_grid1 import make_axes_locatable;

workingDir = "/raid/Thomas_TOPILKO";
regionIdPath = "/home/thomas.topilko/Documents/ClearMap_GUI/ClearMap_Ressources/Regions_annotations/regions IDs.csv";
voxelNumbersPath = "/home/thomas.topilko/Documents/ClearMap_GUI/Tomek_Scripts/Voxel_Numbers_Annotation.xlsx";
projectionTablePath = os.path.join(workingDir, "counts-cells_table_Pregnant_vs_Virgin.csv");
matedActivityTablePath = os.path.join(workingDir, "counts-cells_table_mated_cotton_vs_mated_cup.csv");
virginActivityTablePath = os.path.join(workingDir, "counts-cells_table_virgin_cotton_vs_virgin_cup.csv");

IdTableDF = pd.read_csv(regionIdPath, header=0, usecols=["id", "red", "green", "blue", "acronym", "structure_order", "parent_id"]);

voxelTableDF = pd.read_excel(voxelNumbersPath, header=0, usecols=["id", "voxel"]);

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
    activityInRegion1 = testedGroup[" mean1"].loc[posRegionInActivityTable].values.tolist();
    
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
    
    if availableData :
        
        activityInRegion = [activityInRegion0[0] - activityInRegion1[0]];
        
        activityInRegionNorm = [activityInRegion[0]/voxelsInRegion[0]];
        densityInRegionNorm = [densityInRegion[0]/voxelsInRegion[0]];
        
        if redInRegion[0] != greenInRegion[0] != blueInRegion[0]  and\
        (redInRegion[0], greenInRegion[0], blueInRegion[0]) not in [(176,255,184)] :
        
            sink.append([region, densityInRegionNorm[0], activityInRegionNorm[0], redInRegion[0],\
                               greenInRegion[0], blueInRegion[0], acronymInRegion[0], structureOrderInRegion[0],\
                               parentIdInRegion[0]]);
    
###############################################################################
###############################################################################
colors = np.array(np.array(sink)[:,3:6], dtype=np.int64)/255.

fig = plt.figure(figsize=(10,7), dpi=200.);
ax0 = plt.subplot(1,1,1);

ax0.scatter(np.array(np.array(sink)[:,1], dtype=np.float64),\
            np.array(np.array(sink)[:,2], dtype=np.float64),\
            color = colors);
            
for x, y, t in zip(np.array(np.array(sink)[:,1], dtype=np.float64),\
                   np.array(np.array(sink)[:,2], dtype=np.float64),\
                   np.array(sink)[:,6]) :
    
    ax0.text(x,y,t,fontsize=10);

ax0.set_xlabel("(n puncta from Projections in region) / (n Voxels in region)");
ax0.set_ylabel("(n cFOS+ cells in region) / (n Voxels in region)");
ax0.set_title("Pregnant Cotton vs Pregnant Cup");
ax0.set_ylim(-0.05,0.15);
ax0.set_xlim(-0.01,0.21);
ax0.axhline(0,0, np.inf)

###############################################################################
###############################################################################

from natsort import natsorted;

common = np.intersect1d(np.array(plotMatrix1)[:,0], np.array(plotMatrix2)[:,0]);

mask1 = [True if x in common else False for x in np.array(plotMatrix1)[:,0] ];
mask2 = [True if x in common else False for x in np.array(plotMatrix2)[:,0] ];

plotMatrix1New = np.array(plotMatrix1)[mask1];
plotMatrix2New = np.array(plotMatrix2)[mask2];

colors = np.array(np.array(plotMatrix1)[:,3:6][mask1], dtype=np.int64)/255.

fig = plt.figure(figsize=(10,7), dpi=200.);
ax0 = plt.subplot(1,1,1);

for r1, r2, c in zip(plotMatrix1New, plotMatrix2New, colors) :
    
    dx = float(r2[1])-float(r1[1])
    dy = float(r2[2])-float(r1[2])
    
    ax0.arrow(float(r1[1]), float(r1[2]), dx, dy,  color=c, width=0.0005);
    
ax0.set_ylim(-0.1,0.2);
ax0.set_xlim(-0.01, 0.3);
            
for x, y, t in zip(np.array(np.array(plotMatrix1)[:,1], dtype=np.float64),\
                   np.array(np.array(plotMatrix1)[:,2], dtype=np.float64),\
                   np.array(plotMatrix1)[:,6]) :
    
    ax0.text(x,y,t,fontsize=5);

###############################################################################
###############################################################################
#Connectivity plot

center = (0, 0);
brainRegions = np.array(sink)[:,0];
nBrainRegions = len(brainRegions);

regionsStructureOrder = np.array(np.array(sink)[:,7], dtype="int64");
ordered = np.argsort(regionsStructureOrder);
ordered = np.argsort(ordered);

orderedColors = np.full_like(colors, np.array([0.,0.,0.]));
orderedStructureOrder = np.full_like(regionsStructureOrder, 0);
orderedRegionNames = np.full_like(np.array(sink)[:,6], "");
orderedDensity = np.full_like(np.array(sink)[:,1], 0);
orderedActivity = np.full_like(np.array(sink)[:,2], 0);

for x, c, n, o, d, a in zip(ordered, colors, np.array(sink)[:,6],\
                      regionsStructureOrder, np.array(sink)[:,1],\
                      np.array(sink)[:,2]) :
    
    orderedStructureOrder[x] = o;
    orderedColors[x] = c;
    orderedRegionNames[x] = n;
    orderedDensity[x] = d;
    orderedActivity[x] = a;

orderedActivity = np.array(orderedActivity, dtype="float64");
#minima = min(orderedActivity);
#orderedActivityNew = orderedActivity+abs(minima);
#minima = min(orderedActivity);
#maxima = max(orderedActivity);
norm = matplotlib.colors.Normalize(vmin=-0.1, vmax=0.1, clip=True);
#norm = matplotlib.colors.LogNorm(vmin=0, vmax=0.1, clip=True);
mapper = cm.ScalarMappable(norm=norm, cmap=cm.coolwarm);

fig = plt.figure(figsize=(7,7), dpi=200.);
ax0 = plt.subplot(1, 1, 1);
delta = 0.07;
factor = 5;
#ax0.scatter(0, 0, c=(255./255., 178./255., 253./255.), s=100)

#pos = 0;
radius = 5.;
startAngle = 0;
angle = (360*(np.pi/180)) /nBrainRegions;

cmapBrainRegions = LinearSegmentedColormap.from_list("BR", orderedColors);

for n in np.arange(0, len(orderedColors), 1) :
    
    color = cmapBrainRegions(n/len(orderedColors));
    xi, yi = (radius+0.7)*np.cos(startAngle), (radius+0.7)*np.sin(startAngle);
    xf, yf = (radius+0.7)*np.cos(startAngle+angle), (radius+0.7)*np.sin(startAngle+angle);
    
    ax0.plot((xi,xf), (yi,yf), c=color, lw=3);
    
    startAngle += angle;
    
startAngle = 0;
scatter = [];
    
for c, n, d, a in zip(orderedColors, orderedRegionNames, orderedDensity, orderedActivity) :
    
    mappedColor = mapper.to_rgba(float(a))
    x = radius*np.cos(startAngle);
    y = radius*np.sin(startAngle);
    
    xLine = (radius-delta)*np.cos(startAngle);
    yLine = (radius-delta)*np.sin(startAngle);
    
    plt.plot([0,xLine], [0,yLine], 'k-', color="black", lw=float(d)*factor)
    sc = ax0.scatter(x, y, c=mappedColor, s=10);
    scatter.append(sc)
    
    xn, yn = ((radius+0.35)*np.cos(startAngle), (radius+0.35)*np.sin(startAngle))
    ax0.text(xn, yn, n, rotation = np.rad2deg(startAngle), ha="center", va="center",\
             fontsize=5);
    
#    xo, yo = (11.3*np.cos(startAngle), 11.3*np.sin(startAngle))
    
#    patch = Rectangle((xo, yo), 0.2, 0.2, alpha=1, color=c);
#    transform = matplotlib.transforms.Affine2D().rotate_deg(np.rad2deg(startAngle)) + ax0.transData;
#    patch.set_transform(transform);
#    ax0.add_patch(patch);
#    ax0.scatter(xo, yo, c=c, s=30, marker="8");
             
    startAngle += angle;
             
#    pos += 1;
    
#    if pos == 10 :
#        break

ax0.set_xticks([]);
ax0.set_yticks([]);

handles = ["Density of projection",
           Line2D([0], [0], color="black", lw=0.2*factor),
           Line2D([0], [0], color="black", lw=0.1*factor),
           Line2D([0], [0], color="black", lw=0.05*factor)]

ax0.legend(handles=handles, labels=["Density of projection","20%","10%","5%"], loc=2, fontsize=5);
#divider = make_axes_locatable(ax0)
#cax = divider.append_axes('right', size='5%', pad=0.05);
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
axins = inset_axes(ax0,
                   width="2.5%",  # width = 5% of parent_bbox width
                   height="50%",  # height : 50%
                   loc='upper left',
                   bbox_to_anchor=(1.01, 0., 1, 1),
                   bbox_transform=ax0.transAxes,
                   borderpad=0,
                   )

cbar = plt.colorbar(mapper, cax=axins);
cbar.ax.tick_params(labelsize=5);

#ax0.set_title("Density of efferent projections form the EWcp neurons \n\
#              & the fold change of cFOS activity between Pregnant Females")
#plt.tight_layout();
#fig.colorbar(sc, cax=ax0, orientation='vertical');

plt.savefig(os.path.join(workingDir, "Projection_Density_&_Activity_Pregnant_Cotton_Pregnant_Cup.png"), dpi=200.);
plt.savefig(os.path.join(workingDir, "Projection_Density_&_Activity_Pregnant_Cotton_Pregnant_Cup.svg"), dpi=200.);

















