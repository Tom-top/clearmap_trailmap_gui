#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 14:57:14 2019

@author: thomas.topilko
"""

import os;
import pandas as pd;
import matplotlib.pyplot as plt;
import math;

File = "/mnt/raid/Test/counts-intensity_table_virgin_cotton_vs_virgin_cup.csv";
FileID = "/home/thomas.topilko/Documents/ClearMap_GUI/ClearMap_Ressources/Regions_annotations/regions IDs.csv";

grp1 = "virgin_cotton"
grp2 = "virgin_cup"

f = pd.read_csv(File);
fID = pd.read_csv(FileID);

threshVal = 0.08;
thresh = -math.log(0.05, 10);
threshTendance = -math.log(threshVal, 10);

fig = plt.figure(figsize=(8,8));
ax = plt.subplot(111);

maxXVal = 0;
maxYVal = 0;

highlightID = [975,226,867,1048,661,286,307,938,452,263,689]

for m1, m2, pval, ID, name in zip(f[" mean1"], f[" mean2"], f[" pvalue"], f["id"], f[" name"]) :
    
    foldChange = abs(m1-m2);
    name = unicode(name, "utf-8");
    
    try :
        foldChange = math.log(foldChange, 2);
    except :
        foldChange = 0;
    
    logPVal = -math.log(pval,10);
    
    if logPVal > maxXVal :
        
        maxXVal = logPVal;
        
    if foldChange > maxYVal :
        
        maxYVal = foldChange;
    
    if pval < threshVal :    
        
        pos = fID["id"] == ID;
        index = [i for i, x in enumerate(pos) if x];
#        print(index, ID)
        print(name)
        
        r, g, b = float(fID["red"][index[0]])/255, float(fID["green"][index[0]])/255, float(fID["blue"][index[0]])/255;
        
        if foldChange > 0 :
            
            ax.scatter(foldChange, logPVal, color=[(r,g,b)], alpha=0.5, marker="o", edgecolor="red");
            
            if ID in highlightID :
                ax.text(foldChange, logPVal, name, fontsize=6.)
        
        elif foldChange < 0 :
            
            ax.scatter(foldChange, logPVal, color=[(r,g,b)], alpha=0.5, marker="o", edgecolor="blue");
            if ID in highlightID :
                ax.text(foldChange, logPVal, name, fontsize=6.)
            
    else :
        
        if foldChange > 0 :
        
            ax.scatter(foldChange, logPVal, color="gray", alpha=0.5, marker="o");
        
        elif foldChange < 0 :
            
            ax.scatter(foldChange, logPVal, color="gray", alpha=0.5, marker="o");

from matplotlib.lines import Line2D;
      
customLeg = [Line2D([0], [0], color="white", markeredgecolor="red", marker="o", alpha=0.5, linewidth=0),
             Line2D([0], [0], color="white", markeredgecolor="blue", marker="o", alpha=0.5,linewidth=0),
             Line2D([0], [0], color="gray", markeredgecolor="white", marker="o", alpha=0.5,linewidth=0)]       

plt.xlim(-maxYVal-maxYVal*0.1,maxYVal+maxYVal*0.1);
plt.ylim(-0.1,maxXVal+1);
plt.plot([-maxYVal-maxYVal*0.1,maxYVal+maxYVal*0.1], [thresh, thresh], "--", color='red', linewidth=1.);
plt.plot([-maxYVal-maxYVal*0.1,maxYVal+maxYVal*0.1], [threshTendance, threshTendance], "--", color='red', alpha=0.5, linewidth=1.);
plt.plot([0,0], [-0.1, maxXVal+1], "-", color='gray', alpha=0.5, linewidth=1.);
plt.legend(handles=customLeg, labels=["Up","Down","N.S"], loc="upper left", shadow=True, fontsize=8.);
plt.xlabel("log2(fold-change)");
plt.ylabel("-log10(p-value)");
plt.title("{0} & {1}\np-value in function of cFOS fold-change".format(grp1,grp2));

plt.savefig(os.path.join('/home/thomas.topilko/Desktop',"{0} & {1}\np-value in function of cFOS fold-change.png".format(grp1,grp2)), dpi=100.)
    
    