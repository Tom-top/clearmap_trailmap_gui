#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 15:23:49 2019

@author: ckirst
"""
### Graph construction

#%% init
import os;

clearMapPath = "/home/thomas.topilko/Documents/ClearMap_2";

if os.getcwd() != clearMapPath : 
    os.chdir(clearMapPath);
    
from ClearMap.Scripts.initialize import *;
import ClearMap.ParallelProcessing.DataProcessing.LargeData as ld;
from ClearMap.Utils.Timer import Timer;
    
filesToClean = ["data_layout_aligned.lyt", "data_layout_aligned_axis.lyt", "data_layout_placed.lyt",\
                "data_stitched.npy"];
    
def cleanNpyFiles(directory, files) :
    
    for f in files :
        
        if os.path.exists( os.path.join(directory, f) ) :
            
            os.remove( os.path.join(directory, f) );
            
    for d in os.listdir(directory) :
        
        if len(d.split("_")) > 1 :
        
            if d.split("_")[1] == "cfos" :
                
                path = os.path.join(directory, d);
                
                for f in os.listdir( path ) :
                    
                    if f.split(".")[-1] == "npy" :
                        
                         os.remove( os.path.join(path, f) );
                         
directory = '/raid/Thomas_TOPILKO/190701-21';

for d in os.listdir(directory) :
    
    if len(d.split("_")) > 1 :
        
        if d.split("_")[1] == "cfos" :
            
            timeStamp = d.split("_")[2];
            expression = os.path.join(d, "{0}_cfos_UltraII[<Y,2> x <X,2>]_C00.ome.tif".format(timeStamp))
            
        elif d.split("_")[1] == "auto" :
            
            timeStamp = d.split("_")[2];
            AutoFile = os.path.join( directory, d, "{0}_auto_UltraII_C00_xyz-Table <Z,4>.ome.tif".format(timeStamp) )
                                                                                                                                                      
ws = wsp.Workspace(name = 'data', directory = directory, expression=expression); 
              
ws.info

import numpy as np

import ClearMap.Settings as settings
import ClearMap.IO.IO as io

import ClearMap.ImageProcessing.Skeletonization.Skeletonize as skl
import ClearMap.Alignment.Resampling as res
import ClearMap.Alignment.Elastix as elx  
import ClearMap.Analysis.Measurements.MeasureRadius as mr
import ClearMap.Analysis.Graphs.GraphProcessing as gp
import ClearMap.Visualization.Vispy.GraphVisual as gv
import ClearMap.Analysis.Graphs.GraphGt as ggt

import ClearMap.Visualization.Plot3d as p3d
import ClearMap.ImageProcessing.Vasculature_final_2019_02_02 as vasc
import ClearMap.ImageProcessing.Vasculature_arteries as vasc_arteries
import ClearMap.Analysis.Measurements.Voxelization as vox


import ClearMap.Alignment.Annotation as ano

import ClearMap.Analysis.Measurements.MeasureRadius as mr
import ClearMap.Analysis.Measurements.MeasureExpression as me

import ClearMap.Visualization.Color as col

TemplateFile = '/home/thomas.topilko/Documents/ClearMap_GUI/ClearMap_Ressources/25um_Autofluo_Reference/Template_1-279.tif'
#ano.set_annotation_file('/home/nicolas.renier/Documents/ClearMap_Ressources/annotation_25_1-246R.nrrd')               


############################################################################################################################

#%% Tile conversion (Conversion dos arquivos tif en nummpy, unha vez rematada, deben eliminarse os orixinais)
###############################################################################              
                  
io.convert_files(ws.file_list('expression', extension = 'tif'), extension = 'npy', verbose = True, processes = 'serial')               

############################################################################################################################
# Layout (Posicionamento a grossomodo, baseado nas coordenadas dos tacos, debemos indicar o solapamento, considerando que 
#o microscopio sempre comete un certo erro por exceso, ou sexa, aumenta o solapamento teorico)
###############################################################################

l = layout = stw.WobblyLayout(expression = ws.filename('expression'), tile_axes = ['X', 'Y'], overlaps = (145, 145));  
                                               
############################################################################################################################
# Place tiles rigid (Mesmo posicionamento ca TeraStitcher, vai optimizando a localizacion dos tacos por pares, e ao final
# fai unha optimizacion/media entre todolos tacos)(require que aumentemos a superficie de solapamento un pouco para funcionar 
# ben, e temos que definir tamen o maximo desprazamento nas 3 coordenadas que imos permitir)
###############################################################################                          

stb.align_layout_rigid_mip(layout, depth = [150, 150, None], max_shifts = [(-30,30),(-30,30),(-20,20)],
                           ranges = [None,None,None],
                           background = (10, 100), clip = 25000, verbose = True, processes = 'serial')

stb.place_layout(layout, method = 'optimization', min_quality = -np.inf, lower_to_origin = True, verbose = True)

# Garda a nova posicion dos tacos

stb.save_layout(ws.filename('layout', postfix='aligned_axis'), layout)
#l = layout = stb.load_layout(ws.filename('layout', postfix='aligned_axis'))


# check axis-alignment (Ilustra o desprazamento realizado para acadar a posicion final no paso anterior)

#a = layout.alignment_from_tile_position((3,4),(3,5))
#a = layout.alignments[0]
#
#plt.figure(1); plt.clf()
#a.plot_mip(depth = (55, 155), max_shifts = [(-30,30),(-30,30),(-20,20)])

############################################################################################################################
# Wobbly alignment (Parte sofisticada do posicionamento, no que se busca a posicion ideal de CADA PLANO en cada taco, co vecinho)
############################################################################### 
l = layout = stb.load_layout(ws.filename('layout', postfix='aligned_axis'))

stw.align_layout(layout, axis_range = (None, None, 3), max_shifts =  [(-30,30),(-30,30),(0,0)], axis_mip = None,
                 validate = dict(method='foreground', valid_range = (10, None), size = None),
                 prepare = dict(method='normalization', clip = None, normalize=True),
                 validate_slice = dict(method='foreground', valid_range = (10,60000), size = 1300),
                 prepare_slice = None,
                 find_shifts = dict(method='tracing', cutoff = 3*np.sqrt(2)),
                 processes = 'serial', verbose = True)


# Garda a nova posicion dos tacos
stb.save_layout(ws.filename('layout', postfix='aligned'), layout)
#l = layout = stb.load_layout(ws.filename('layout', postfix='aligned'));

#
#l.alignment_info((0,0), 1000)                            
############################################################################################################################                                                       
# Placement (Optimizacion final que ainda non comprendemos ben, pero non hai que cambiar ningun parametro)
############################################################################### 

l = layout = stb.load_layout(ws.filename('layout', postfix='aligned'));

stw.place_layout(layout, min_quality = -np.inf, 
                 method = 'optimization', 
                 smooth = dict(method = 'window', window = 'bartlett', window_length = 100, binary = None), 
                 smooth_optimized = dict(method = 'window', window = 'bartlett', window_length = 20, binary = 10),                             
                 fix_isolated = False, lower_to_origin = True,
                 processes = 'serial', verbose = True)


# Garda a nova posicion dos tacos
stb.save_layout(ws.filename('layout', postfix='placed'), layout)
#l = layout = stb.load_layout(ws.filename('layout', postfix='placed'));

#
#l.alignment_info((3,4), 2000,use_displacements = False) 

############################################################################################################################                                   
# Stitching (a felicidade :)
############################################################################### 

l = layout = stb.load_layout(ws.filename('layout', postfix='placed'));

stw.stitch_layout(layout, sink = ws.filename('stitched'),
                  method = 'interpolation', processes='!serial', verbose=True)

############################################################################################################################
# Conversion in tif file
###############################################################################

#dv.plot(ws.filename('stitched'))

io.convert_files(ws.filename('stitched'), extension = 'tif', verbose = True);
cleanNpyFiles(directory, filesToClean);     

