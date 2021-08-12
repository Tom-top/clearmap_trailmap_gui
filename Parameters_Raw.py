#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 20 09:36:06 2018
@author: thomas.topilko
"""

import sys;
import re;

import ClearMap.Alignment.Resampling as rs;
import ClearMap.Alignment.Stitching as st;
import ClearMap.Alignment.Elastix as elx;
import ClearMap.Tomek_Utilities as ut;
import os;
from ClearMap.Utils.ParameterTools import joinParameter;
from ClearMap.ImageProcessing.CellDetection import detectCells;
import ClearMap.IO as io;
from ClearMap.Analysis.Statistics import thresholdPoints;
from ClearMap.Alignment.Resampling import resamplePoints;
from ClearMap.Alignment.Elastix import transformPoints;
from ClearMap.Analysis.Voxelization import voxelize;

clearmapDir = "";
clearmapRessourcesDir = "";
workdir = "";
cFosfile = "";
autodir = "";
AutoFluoFile = "";

FinalOrientation = (1,2,3);

DataResolution = (1.625,1.625,6);
# DataResolution = (5,5,6);
cFosDetectionRange = {"x" : all, "y" : all, "z" : all};

AutoResolution = (5,5,6);
TemplateResolution10m = (10,10,10);
TemplateResolution = (25,25,25);

Upsample = True;

TemplateFile10m = os.path.join(clearmapRessourcesDir+"/10um_Autofluo_Reference/","average_template_10_horizontal.tif");
TemplateFile = os.path.join(clearmapRessourcesDir+"/25um_Autofluo_Reference/","template_Youenn_1-288.tif");
#TemplateFile = os.path.join(clearmapRessourcesDir+"/25um_Autofluo_Reference/","Template_1-279_No_Bulb.tif");
TemplateFileNrrd = os.path.join(clearmapRessourcesDir+"/25um_Autofluo_Reference/","template_horizontal_1-289_25.nrrd");
AnnotationFile = os.path.join(clearmapRessourcesDir+"/Regions_annotations/","ABA_25um_annotation_Youenn.tif");
AnnotationFileNrrd = os.path.join(clearmapRessourcesDir+"/Regions_annotations/","annotation_279_new.nrrd");
    
try :
  
  workdir;
  
except :
  
  print(ut.coloredMessage('[/!\ Warning] workdir is not defined !','darkred'));
  workdir = None;
  
try :
  
  cFosfile;
  
except :
  
  print(ut.coloredMessage('[/!\ Warning] cFosfile is not defined !','darkred'));
  cFosfile = None;
  
try :
  
  autodir;
  
except :
  
  print(ut.coloredMessage('[/!\ Warning] autodir is not defined !','darkred'));
  autodir = None;
  
try :
  
  AutoFluoFile;
  
except :
  
  print(ut.coloredMessage('[/!\ Warning] AutoFluoFile is not defined !','darkred'));
  AutoFluoFile = None;
  
outdir = os.path.join(workdir, 'stitched');

analysisDir = None;
#os.path.join(os.path.dirname(workdir), "Analysis_40_900");
    
cFosTransformedDir = os.path.join(workdir,"elastix_cfos_transformed");


ResamplingParametercFos = {
    "processes"  : 24,
    "source"     : os.path.join(outdir, 'Z\d{4}.tif'),
    "sink"       : os.path.join(workdir, 'cfos_resampled.tif'),
    "resolutionSource" : DataResolution,
    "resolutionSink"   : TemplateResolution,
    "orientation"      : FinalOrientation,
    "crop" : None,
};
   
ResamplingParameterAutoFluo = {
    "processes"  : 24,
    "source"     : AutoFluoFile,
    "sink"       : os.path.join(workdir, 'autofluo_resampled.tif'),
    "resolutionSource" : AutoResolution,
    "resolutionSink"   : TemplateResolution,
    "orientation"      : FinalOrientation
};
        
    
UpsamplingParameterAnnotation = {
    "processes"  : 24,
    "source"     : None, #os.path.join(cFosTransformedDir,'cfos_transformed.tif')
    "sink"       : os.path.join(workdir,'upsampled_annotation.tif'), #os.path.join(workdir,'upsampled_cfos.tif')
    "dataSizeSink" : None, #(x,y,z)
    "resolutionSource" : TemplateResolution, #TemplateResolution10m
    "resolutionSink"   : TemplateResolution10m,
    "orientation"      : FinalOrientation,
    "interpolation" : None,
};
   
ChannelsAlignmentParameter = {            
    #moving and reference images
    "movingImage" : ResamplingParameterAutoFluo["sink"], #os.path.join(workdir,'autofluo_resampled.tif')
    "fixedImage"  : ResamplingParametercFos["sink"], #ResamplingParametercFos["sink"]
    
    #elastix parameter files for alignment
    "affineParameterFile"  : os.path.join(clearmapRessourcesDir+"/Parameter_files/","Par0000affine.txt"),
    "bSplineParameterFile" : None,
    
    #directory of the alignment result
    "resultDirectory" :  os.path.join(workdir, 'elastix_cfos_to_auto')
    }; 
 

TemplateAlignmentParameter = {            
    #moving and reference images
    "movingImage" : TemplateFile,#TemplateFile
    "fixedImage"  : os.path.join(workdir,'autofluo_resampled.tif') ,
    
    #elastix parameter files for alignment
    "affineParameterFile"  : os.path.join(clearmapRessourcesDir+"/Parameter_files/","Par0000affine_acquisition.txt"),
    "bSplineParameterFile" : os.path.join(clearmapRessourcesDir+"/Parameter_files/","Par0000bspline.txt"),
    
    #directory of the alignment result
    "resultDirectory" :  os.path.join(workdir, 'elastix_template_to_auto')
    }; 

TransformFOS2Template = {            
        
    "source" : os.path.join(ChannelsAlignmentParameter["resultDirectory"],"result.0.mhd"),
    "sink"  : "cFOS_to_annotation.tif",
    "transformDirectory" : None,
    "resultDirectory" : os.path.join(workdir, "transform_cfos_to_template"),
    }; 

######################### Cell Detection Parameters using custom filters


ImageProcessingMethod = "SpotDetection";


correctIlluminationParameter = {
    "flatfield"  : None,  # (str, True or None)  flat field intensities, if None d onot correct image for illumination, if True the 
    "background" : None, # (str, None or array) background image as file name or array, if None background is assumed to be zero
    "scaling"    : "Mean", # (str or None)        scale the corrected result by this factor, if 'max'/'mean' scale to keep max/mean invariant
    "save"       : None,       # (str or None)        save the corrected image to file
    "verbose"    : False    # (bool or int)        print / plot information about this step 
}

removeBackgroundParameter = {
    "size"    : (7,7),  #(10,10)# size for the structure element of the morphological opening
    "save"    : None, #os.path.join(workdir, 'background/Z\d{4}.tif'),     # file name to save result of this operation
    "verbose" : False  # print / plot information about this step       
}


filterDoGParameter = {
    "size"    : None,        # (tuple or None)      size for the DoG filter if None, do not correct for any background
    "sigma"   : None,        # (tuple or None)      std of outer Guassian, if None autmatically determined from size
    "sigma2"  : None,        # (tuple or None)      std of inner Guassian, if None autmatically determined from size
    "save"    : None,        # (str or None)        file name to save result of this operation if None dont save to file 
    "verbose" : False      # (bool or int)        print / plot information about this step
}

findExtendedMaximaParameter = {
    "hMax"      : None,            # (float or None)     h parameter for the initial h-Max transform, if None, do not perform a h-max transform
    "size"      : 3,             # #5 (tuple)             size for the structure element for the local maxima filter
    "threshold" : 0,        # (float or None)     include only maxima larger than a threshold, if None keep all localmaxima
    "save"      : None,         # (str or None)       file name to save result of this operation if None dont save to file 
    "verbose"   : False       # (bool or int)       print / plot information about this step
}

findIntensityParameter = {
    "method" : 'Max',       # (str, func, None)   method to use to determine intensity (e.g. "Max" or "Mean") if None take intensities at the given pixels
    "size"   : (3,3,3)       #(3,3,3) (tuple)             size of the box on which to perform the *method*
}

detectCellShapeParameter = {
    "threshold" : 200,     # (float or None)      threshold to determine mask, pixel below this are background if None no mask is generated
    "save"      : None, #os.path.join(workdir, 'shape/Z\d{4}.tif'),    #os.path.join(workdir, 'shape/Z\d{4}.tif'), # (str or None)        file name to save result of this operation if None dont save to file 
    "verbose"   : False      # (bool or int)        print / plot information about this step if None take intensities at the given pixels
}

## Parameters for cell detection using spot detection algorithm 
detectSpotsParameter = {
    "correctIlluminationParameter" : correctIlluminationParameter,
    "removeBackgroundParameter"    : removeBackgroundParameter,
    "filterDoGParameter"           : filterDoGParameter,
    "findExtendedMaximaParameter"  : findExtendedMaximaParameter,
    "findIntensityParameter"       : findIntensityParameter,
    "detectCellShapeParameter"     : detectCellShapeParameter
}

thresholdPointParameter = {
        "threshold" : (10,900),
        };


#################### Heat map generation

##Voxelization

VoxelizationFile = os.path.join(workdir, 'points_voxelized.tif');

## Parameter to calculate density voxelization
voxelizeParameter = {
    #Method to voxelize
    "method" : 'Spherical', # Spherical,'Rectangular, Gaussian'
       
    # Define bounds of the volume to be voxelized in pixels
    "size" : (15,15,15),  

    # Voxelization weigths (e/g intensities)
    "weights" : None
};


StackProcessingParameter = {
    #max number of parallel processes
    "processes" : 17,
   
    #chunk sizes
    "chunkSizeMax" : 20, #20 10
    "chunkSizeMin" : 15, #15 6
    "chunkOverlap" : 5, #5 4
    "verbose" : False,

    #optimize chunk size and number to number of processes
    "chunkOptimization" : True,
    
    #increase chunk size for optimization (True, False or all = automatic)
    "chunkOptimizationSize" : all,
   
    "processMethod" : "parallel"   
   };

SpotDetectionParameter = {
    "source" : os.path.join(outdir, 'Z\d{4}.tif'),
    "sink"   : (os.path.join(workdir, 'cells-allpoints.npy'),  os.path.join(workdir,  'intensities-allpoints.npy')),
    "detectSpotsParameter" : detectSpotsParameter,
};
        
#print(SpotDetectionParameter["source"]);
    
SpotDetectionParameter = joinParameter(SpotDetectionParameter, cFosDetectionRange);

ImageProcessingParameter = joinParameter(StackProcessingParameter, SpotDetectionParameter);

cellsCheckFolder = os.path.join(workdir, "cells_check");
FilteredCellsFileCheck = (os.path.join(cellsCheckFolder, 'cells.npy'), os.path.join(cellsCheckFolder,  'intensities.npy'));

mouseNumber = os.path.basename(workdir);

if analysisDir != None :
    FilteredCellsFile = (os.path.join(analysisDir, 'cells_{0}.npy'.format(mouseNumber)), os.path.join(analysisDir,  'intensities_{0}.npy'.format(mouseNumber)));
    TransformedCellsFile = os.path.join(analysisDir, 'cells_transformed_to_Atlas_{0}.npy'.format(mouseNumber));
    
else :
    FilteredCellsFile = (os.path.join(workdir, 'cells.npy'), os.path.join(workdir,  'intensities.npy'));
    TransformedCellsFile = os.path.join(workdir, 'cells_transformed_to_Atlas.npy');

ResamplePointsParameters = {
    "pointSource" : FilteredCellsFile[0],
    "dataSizeSource" : os.path.join(outdir, 'Z\d{4}.tif'),
    "pointSink"   : None,
    "resolutionSource" : DataResolution,
    "resolutionSink"   : TemplateResolution,
    "orientation"      : FinalOrientation
};
