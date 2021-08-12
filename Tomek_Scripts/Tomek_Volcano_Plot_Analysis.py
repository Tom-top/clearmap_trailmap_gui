#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 14:57:14 2019

@author: thomas.topilko
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import math

File = "/network/lustre/dtlake01/renier/Thomas/thomas.topilko/Experiments/Nesting_Project/Behavior/Global_Analysis_cFOS/Thumbnail_pvalues/Raw_Data/Females/counts-cells_table_pregnant_cup_vs_virgin_cup.csv"
FileID = "/home/thomas.topilko/PycharmProjects/clearmap_gui_tomek/ClearMap_Ressources/Regions_annotations/regions_IDs_left_right.csv"

grp1 = "pregnant_cup"
grp2 = "virgin_cup"

f = pd.read_csv(File)
fID = pd.read_csv(FileID)

threshVal = 0.08
thresh = -math.log(0.05, 10)
threshTendance = -math.log(threshVal, 10)

# cm = 1 / 2.53
fig = plt.figure(figsize=(5,5), dpi=300)
ax = plt.subplot(111)

maxXVal = 0
maxYVal = 0
marker_size = 30

# highlightID = [975,226,867,1048,661,286,307,938,452,263,689]

for m1, m2, pval, ID, name in zip(f[" mean1"], f[" mean2"], f[" pvalue"], f["id"], f[" name"]) :
    try:
        foldChange = m1/m2
    except:
        foldChange = 0
    name = str(name)
    # name = name.decode("utf-8", "strict")
    try:
        foldChange = math.log(foldChange, 2)
    except:
        foldChange = 0
    logPVal = -math.log(pval,10)
    if logPVal > maxXVal:
        maxXVal = logPVal
    if foldChange > maxYVal:
        maxYVal = foldChange
    if pval < threshVal:
        pos = fID["id"] == ID
        index = [i for i, x in enumerate(pos) if x]
        name = str(fID["acronym"][index[0]])
        r, g, b = float(fID["red"][index[0]])/255, float(fID["green"][index[0]])/255, float(fID["blue"][index[0]])/255
        if foldChange > 0:
            ax.scatter(foldChange, logPVal, color=[(r, g, b)], alpha=1, marker="o", s=marker_size, linewidth=0.2, edgecolor="black")
            # if ID in highlightID:
            ax.text(foldChange, logPVal, name, fontsize=6., ha="center", va="center")
        elif foldChange < 0:
            ax.scatter(foldChange, logPVal, color=[(r, g, b)], alpha=1, marker="o", s=marker_size, linewidth=0.2, edgecolor="black")
            # if ID in highlightID:
            ax.text(foldChange, logPVal, name, fontsize=6., ha="center", va="center")
    else:
        if foldChange > 0:
            ax.scatter(foldChange, logPVal, color="gray", alpha=0.5, marker="o", s=marker_size, linewidth=0.2, edgecolor="black")
        elif foldChange < 0:
            ax.scatter(foldChange, logPVal, color="gray", alpha=0.5, marker="o", s=marker_size, linewidth=0.2, edgecolor="black")

# customLeg = [Line2D([0], [0], color="white", markeredgecolor="red", marker="o", alpha=0.5, linewidth=0),
#              Line2D([0], [0], color="white", markeredgecolor="blue", marker="o", alpha=0.5,linewidth=0),
#              Line2D([0], [0], color="gray", markeredgecolor="white", marker="o", alpha=0.5,linewidth=0)]
customLeg = [Line2D([0], [0], color="gray", markeredgecolor="white", marker="o", alpha=0.5,linewidth=0)]


plt.xlim(-maxYVal-maxYVal*0.1,maxYVal+maxYVal*0.1)
plt.ylim(-0.1,maxXVal+1)
plt.plot([-maxYVal-maxYVal*0.1,maxYVal+maxYVal*0.1], [thresh, thresh], "--", color='red', linewidth=1.)
plt.plot([-maxYVal-maxYVal*0.1,maxYVal+maxYVal*0.1], [threshTendance, threshTendance], "--", color='red', alpha=0.5, linewidth=1.)
plt.plot([0,0], [-0.1, maxXVal+1], "-", color='gray', alpha=0.5, linewidth=1.)
plt.legend(handles=customLeg, labels=["N.S"], loc="upper left", shadow=True, fontsize=8.)
plt.xlabel("log2(fold-change)")
plt.ylabel("-log10(p-value)")
plt.title("{0} & {1} p-value in function\nof the fold-change in projection density".format(grp1, grp2))

plt.savefig(os.path.join('/home/thomas.topilko/Bureau',
                         "{0} & {1}\np-value in function of the fold-change in projection density.png".
                         format(grp1, grp2)), dpi=300.)
plt.savefig(os.path.join('/home/thomas.topilko/Bureau',
                         "{0} & {1}\np-value in function of the fold-change in projection density.svg".
                         format(grp1, grp2)), dpi=300.)
plt.show()


# file = "/home/thomas.topilko/PycharmProjects/clearmap_gui_tomek/ClearMap/Data/ARA2_annotation_info_collapse.csv"
# import pandas as pd

# df = pd.read_csv(file)
# df_new = df.copy()
#
# df["name"] += " left"
# df_new["name"] += " right"
# df_new["id"] += 20000
#
# df_f = pd.concat([df, df_new])
# df_f.to_csv("/home/thomas.topilko/PycharmProjects/clearmap_gui_tomek/ClearMap_Ressources/Regions_annotations/regions_IDs_left_right.csv")