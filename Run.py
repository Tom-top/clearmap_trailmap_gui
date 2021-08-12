# %%
# !/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 17:47:32 2017

@author: thomas.topilko
"""

import ClearMap.Tomek_Utilities as ut;
import sys, os;
from natsort import natsorted;
import time;
from shutil import copyfile;
import fileinput;
import imp;
import shutil;
import numpy as np;

import ClearMap.Alignment.Resampling as rs;
import ClearMap.Alignment.Stitching as st;
import ClearMap.Alignment.Elastix as elx;

from ClearMap.ImageProcessing.CellDetection import detectCells;
import ClearMap.IO as io;
from ClearMap.Analysis.Statistics import thresholdPoints;
from ClearMap.Alignment.Resampling import resamplePoints;
from ClearMap.Alignment.Elastix import transformPoints;
from ClearMap.Analysis.Voxelization import voxelize;

import Parameters;


# %%###########################################################################
### Function that lauches the pipeline
##############################################################################

def LaunchClearMap(**kwargs):
    all_false = 'FFFFFFFF'

    os.remove(os.path.join(kwargs["clearmapDir"], 'Parameters.py'))
    copyfile(os.path.join(kwargs["clearmapDir"], "Parameters_Raw.py"),
             os.path.join(kwargs["clearmapDir"], "Parameters.py"))

    try:

        for folder in natsorted(os.listdir(kwargs["homedir"])):

            auto = False
            cfos = False

            path2Folder = os.path.join(kwargs["homedir"], folder)

            if os.path.isdir(path2Folder):

                if folder.split('_')[0] != folder:
                    exp = folder.split('_')[0]
                elif folder.split('-')[0] != folder:
                    exp = folder.split('-')[0]
                else:
                    exp = folder

                    # raise RuntimeError('No folder(s) with the name {} were found in {} \n'.format(experiment,homedir));
                if exp == kwargs["experiment"]:

                    if kwargs["operations"][folder] != all_false:

                        workdir = path2Folder

                        datadir = None
                        nCFosFiles = None
                        cFosfile = None
                        autodir = None
                        nAutoFiles = None
                        AutoFluoFile = None

                        print(ut.coloredMessage('[INFO] Current Workdir is : %s' % (workdir), 'darkgreen'))

                        for subfolder in os.listdir(os.path.join(kwargs["homedir"], folder)):

                            path2Subfolder = os.path.join(path2Folder, subfolder)

                            if os.path.isdir(path2Subfolder):

                                # timeStamp = subfolder.split('_')[-1]
                                timeStamp = subfolder.split('m')[-1]

                                try:

                                    dataType = subfolder.split('m')[0]

                                except:

                                    dataType = None

                                if dataType in ut.cfosPrefix:

                                    datadir = path2Subfolder
                                    nCFosFiles = len(os.listdir(datadir))
                                    print(ut.coloredMessage('[INFO] Current Datadir is : %s' % (datadir),
                                                            'darkgreen'))
                                    cFosfile = os.path.join(datadir,
                                                            r'' + timeStamp + '_UltraII\[(?P<row>\d{2}) x (?P<col>\d{2})\]_C00.ome.tif')
                                    # cFosfile = os.path.join(datadir, r''+timeStamp+'_%s_UltraII\[(?P<row>\d{2}) x (?P<col>\d{2})\]_C00_xyz-Table Z(?P<z>\d{4}).ome.tif'%(dataType))
                                    # cFosfile = os.path.join(datadir, r''+timeStamp+'_%s_UltraII_C00_UltraII Filter0002.ome.tif'%(dataType));
                                    # cFosfile = os.path.join(datadir, "projections/\d{4}.tif")
                                    print(ut.coloredMessage(
                                        '[INFO] %d %s files have been detected' % (nCFosFiles, dataType),
                                        'darkgreen'))
                                    cfos = True

                                elif dataType in ut.autoPrefix:

                                    autodir = path2Subfolder
                                    nAutoFiles = len(os.listdir(autodir))
                                    print(ut.coloredMessage('[INFO] Current Autodir is : %s' % (autodir),
                                                            'darkgreen'))
                                    AutoFluoFile = os.path.join(autodir,
                                                                timeStamp + '_UltraII_C00_xyz-Table Z\d{4}.ome.tif')
                                    # AutoFluoFile = os.path.join(autodir, timeStamp+'_%s_UltraII_C00_xyz-Table Z\d{4}.ome.tif'%(dataType))
                                    # AutoFluoFile = os.path.join(autodir, timeStamp+'_%s_UltraII_C00_UltraII Filter0000.ome.tif'%(dataType))
                                    # AutoFluoFile = os.path.join(autodir, "autofluo.tif")
                                    print(ut.coloredMessage(
                                        '[INFO] %d %s files have been detected \n' % (nAutoFiles, dataType),
                                        'darkgreen'))
                                    auto = True

                        backgroundThreshPos = None
                        cellSizeThreshPos = None

                        with open(os.path.join(kwargs["clearmapDir"], "Temp/Parameters.py"), 'w') as new_file:

                            with open(os.path.join(kwargs["clearmapDir"], 'Parameters.py')) as old_file:

                                for n, line in enumerate(old_file):

                                    if line.split(" ") == ['detectCellShapeParameter', '=', '{\n']:

                                        backgroundThreshPos = n + 1;

                                    elif line.split(" ") == ['thresholdPointParameter', '=', '{\n']:

                                        cellSizeThreshPos = n + 1;

                                    if n == backgroundThreshPos:

                                        new_file.write("    'threshold' : {0}, \n".format(kwargs["bgtParams"][folder]));

                                    elif n == cellSizeThreshPos:

                                        new_file.write("    'threshold' : {0}, \n".format(kwargs["cstParams"][folder]));

                                    elif line.split(" ")[0] == "cFosDetectionRange":

                                        new_file.write(
                                            "cFosDetectionRange = {{ 'x' : {0}, 'y' : {1}, 'z' : {2} }}; \n".format(
                                                kwargs["cfosX"][folder], \
                                                kwargs["cfosY"][folder], kwargs["cfosZ"][folder]));

                                    elif line.split(" ")[0] == 'clearmapDir':

                                        new_file.write("clearmapDir = '{0}'; \n".format(kwargs["clearmapDir"]));

                                    elif line.split(" ")[0] == 'clearmapRessourcesDir':

                                        new_file.write("clearmapRessourcesDir = '{0}'; \n".format(
                                            kwargs["clearmapRessourcesDir"]));

                                    elif line.split(" ")[0] == 'workdir':

                                        new_file.write("workdir = '{0}'; \n".format(workdir));

                                    elif line.split(" ")[0] == 'cFosfile':
                                        print(cFosfile)
                                        new_file.write("cFosfile = '{0}'; \n".format(cFosfile));

                                    elif line.split(" ")[0] == 'autodir':

                                        new_file.write("autodir = '{0}'; \n".format(autodir));

                                    elif line.split(" ")[0] == 'AutoFluoFile':

                                        new_file.write("AutoFluoFile = '{0}'; \n".format(AutoFluoFile));

                                    elif line.split(" ")[0] == 'FinalOrientation':

                                        new_file.write("FinalOrientation = {0}; \n".format(kwargs["orient"][folder]));

                                    else:

                                        new_file.write(line)

                        os.remove(os.path.join(kwargs["clearmapDir"], 'Parameters.py'));
                        os.rename(os.path.join(kwargs["clearmapDir"], "Temp/Parameters.py"),
                                  os.path.join(kwargs["clearmapDir"], 'Parameters.py'));
                        print(ut.coloredMessage('New Parameter file created', 'darkgreen'));

                        # execfile(os.path.join(clearmapDir,'Parameters.py'));
                        # execfile(os.path.join(clearmapDir,'Run.py'));
                        #                    exec(open(os.path.join(kwargs["clearmapDir"],'Parameters.py')).read());

                        if auto == False:
                            print(ut.coloredMessage('[WARNING] No autofluorescence folder was detected !!', 'darkred'));

                        if cfos == False:
                            print(ut.coloredMessage('[WARNING] No cfos folder was detected !!', 'darkred'));

                        #                    imp.reload(Parameters);

                        launchPipeline(kwargs["operations"][folder], datadir, **kwargs);

    except KeyboardInterrupt:

        print(ut.coloredMessage('[WARNING] The script for mouse {0} was manually aborted !'.format(folder), 'darkred'));
        os.remove(os.path.join(kwargs["clearmapDir"], 'Parameters.py'));
        copyfile(os.path.join(kwargs["clearmapDir"], "Parameters_Raw.py"),
                 os.path.join(kwargs["clearmapDir"], "Parameters.py"))

    else:

        print(ut.coloredMessage('[INFO] The script successfully finished'.format(folder), 'darkgreen'));
        os.remove(os.path.join(kwargs["clearmapDir"], 'Parameters.py'));
        copyfile(os.path.join(kwargs["clearmapDir"], "Parameters_Raw.py"),
                 os.path.join(kwargs["clearmapDir"], "Parameters.py"))


def launchPipeline(folderOperation, datadir, **kwargs):
    exec(open(os.path.join(kwargs["clearmapDir"], 'Parameters.py')).read())
    imp.reload(Parameters);

    print(Parameters.cFosDetectionRange)

    if folderOperation[0] == 'T':

        print("\n");
        print(ut.coloredMessage(ut.titleMessage('Starting stiching'), 'darkgreen'));

        print(Parameters.cFosfile)
        print(os.path.join(Parameters.workdir, 'TeraStitcher_import.xml'))
        importfile = st.xmlImportFile(Parameters.cFosfile,
                                      xmlImportFile=os.path.join(Parameters.workdir, 'TeraStitcher_import.xml'));
        st.importData(importfile);

        try:
            os.remove(os.path.join(datadir, 'mdata.bin'));
        except:
            pass;
        print(datadir)
        alignfile = st.alignData(importfile,
                                 search=(48, 59, 8),
                                 xmlResultFile=os.path.join(datadir, 'TeraStitcher_align.xml'),
                                 subRegion=((all, all), (all, all), (700, 750)))

        projectfile = st.projectDisplacements(os.path.join(datadir, 'TeraStitcher_align.xml') \
                                              , xmlResultFile=os.path.join(datadir, 'TeraStitcher_project.xml'));

        thresfile = st.thresholdDisplacements(projectfile,
                                              xmlResultFile=os.path.join(datadir, 'TeraStitcher_threshold.xml'));

        placefile = st.placeTiles(thresfile, xmlResultFile=os.path.join(datadir, 'TeraStitcher_place.xml'));

        result = st.stitchData(placefile, resultPath=os.path.join(Parameters.outdir, 'Z\d{4}.tif') \
                               , bitDepth=16, algorithm='SINBLEND');

        if kwargs["playSound"]:
            ut.playSound(1);

    # ''' Resampling '''

    if folderOperation[1] == 'T':

        print("\n");
        print(ut.coloredMessage(ut.titleMessage('Starting resampling Autofluorescence'), 'darkgreen'));

        print(Parameters.ResamplingParameterAutoFluo["source"])
        rs.resampleData(**Parameters.ResamplingParameterAutoFluo);

        if kwargs["playSound"]:
            ut.playSound(1);

    if folderOperation[2] == 'T':

        print("\n");
        print(ut.coloredMessage(ut.titleMessage('Starting resampling cFos'), 'darkgreen'));

        rs.resampleData(**Parameters.ResamplingParametercFos);

        if kwargs["playSound"]:
            ut.playSound(1);

    # ''' Aligning '''

    if folderOperation[3] == 'T':

        time.sleep(5)

        print("\n");
        print(ut.coloredMessage(ut.titleMessage('Starting alignement Auto/cFOS'), 'darkgreen'));

        resultDirectory = elx.alignData(**Parameters.ChannelsAlignmentParameter);

        if kwargs["playSound"]:
            ut.playSound(1);

    if folderOperation[4] == 'T':

        time.sleep(5)

        print("\n");
        print(ut.coloredMessage(ut.titleMessage('Starting alignement Auto/Template'), 'darkgreen'));

        resultDirectory = elx.alignData(**Parameters.TemplateAlignmentParameter);

        if kwargs["playSound"]:
            ut.playSound(1);

    # ''' Cell Detection '''

    if folderOperation[5] == 'T':

        print("\n");
        print(ut.coloredMessage(ut.titleMessage('Starting cell detection'), 'darkgreen'));

        detectCells(**Parameters.ImageProcessingParameter);

        if kwargs["playSound"]:
            ut.playSound(1);

    if folderOperation[6] == 'T':

        print("\n")
        print(ut.coloredMessage(ut.titleMessage('Starting Heatmap generation'), 'darkgreen'))

        ###############################################################################
        # cell check
        ###############################################################################

        # cellsCheckFolder = os.path.join(Parameters.workdir, "cells_check")
        # if not os.path.exists(cellsCheckFolder):
        #     os.mkdir(cellsCheckFolder)
        #     os.mkdir(os.path.join(cellsCheckFolder, "background"))
        #     os.mkdir(os.path.join(cellsCheckFolder, "cell_shape"))
        # else:
        #     shutil.rmtree(cellsCheckFolder)
        #     os.mkdir(cellsCheckFolder)
        #     os.mkdir(os.path.join(cellsCheckFolder, "background"))
        #     os.mkdir(os.path.join(cellsCheckFolder, "cell_shape"))
        #
        # Parameters.ImageProcessingParameter["sink"] = (os.path.join(cellsCheckFolder, 'cells-allpoints.npy'),
        #                                                os.path.join(cellsCheckFolder, 'intensities-allpoints.npy'))
        # #        Parameters.ImageProcessingParameter["detectSpotsParameter"]["correctIlluminationParameter"]["save"] = os.path.join(cellsCheckFolder, 'illumination_correction/Z\d{4}.tif')
        # #        Parameters.ImageProcessingParameter["detectSpotsParameter"]["removeBackgroundParameter"]["save"] = os.path.join(cellsCheckFolder, r'background/Z\d{4}.tif')
        # #        Parameters.ImageProcessingParameter["detectSpotsParameter"]["detectCellShapeParameter"]["save"] = os.path.join(cellsCheckFolder, r'cell_shape/Z\d{4}.tif')
        # Parameters.ImageProcessingParameter["detectSpotsParameter"]["removeBackgroundParameter"]["save"] = None
        # Parameters.ImageProcessingParameter["detectSpotsParameter"]["detectCellShapeParameter"]["save"] = None
        #
        # #        Parameters.ImageProcessingParameter["processMethod"] = "sequential";
        # detectCells(**Parameters.ImageProcessingParameter)
        #
        # points, intensities = io.readPoints(Parameters.ImageProcessingParameter["sink"])
        # points, intensities = thresholdPoints(points, intensities,
        #                                       threshold=Parameters.thresholdPointParameter["threshold"], row=(3, 3))
        # io.writePoints(Parameters.FilteredCellsFileCheck, (points, intensities))
        #
        # import ClearMap.Visualization.Plot as plot
        # pointSource = os.path.join(Parameters.workdir, Parameters.FilteredCellsFileCheck[0])
        # data = plot.overlayPoints(Parameters.ResamplingParametercFos["source"], pointSource, pointColor=None,
        #                           **Parameters.cFosDetectionRange)
        # io.writeData(os.path.join(cellsCheckFolder, 'cells_check.tif'), data)

        ###############################################################################
        # heatmap
        ###############################################################################

        if Parameters.analysisDir != None:
            if not os.path.exists(Parameters.analysisDir):
                os.mkdir(Parameters.analysisDir)

        points, intensities = io.readPoints(Parameters.ImageProcessingParameter["sink"])
        points, intensities = thresholdPoints(points, intensities,
                                              threshold=Parameters.thresholdPointParameter["threshold"],
                                              row=(3, 3))  # 200,9000#100,900
        io.writePoints(Parameters.FilteredCellsFile, (points, intensities))

        points = io.readPoints(Parameters.FilteredCellsFile[0])
        points = resamplePoints(**Parameters.ResamplePointsParameters)
        points = transformPoints(points, transformDirectory=Parameters.ChannelsAlignmentParameter["resultDirectory"],
                                 indices=False, resultDirectory=None, binary=False)
        points = transformPoints(points, transformDirectory=Parameters.TemplateAlignmentParameter["resultDirectory"],
                                 indices=False, resultDirectory=None, binary=False)
        io.writePoints(Parameters.TransformedCellsFile, points)

        points = io.readPoints(Parameters.TransformedCellsFile)
        intensities = io.readPoints(Parameters.FilteredCellsFile[1])

        vox = voxelize(points, Parameters.TemplateFile,
                       **Parameters.voxelizeParameter)  # TemplateFile TemplateFileNrrd

        if Parameters.analysisDir == None:
            io.writeData(os.path.join(Parameters.workdir, 'cells_heatmap.tif'), vox.astype('int32'))
        else:
            io.writeData(os.path.join(Parameters.analysisDir, 'cells_heatmap_{0}.tif'.format(Parameters.mouseNumber)),
                         vox.astype('int32'))

        if kwargs["playSound"]:
            ut.playSound(1)

    if folderOperation[7] == 'T':

        ############################################################################################################################################
        ############################################################################################################################################
        # Test

        ##        annotationFile = Parameters.AnnotationFile;
        ##
        ##        Path2Alignement = os.path.join(Parameters.workdir,'elastix_template_to_auto');
        ##        copyfile(os.path.join(Path2Alignement,'TransformParameters.1.txt'), os.path.join(Path2Alignement,'TransformParameters_noInterpolation.1.txt'));
        ##        transformFile = os.path.join(Path2Alignement,'TransformParameters_noInterpolation.1.txt');
        ##
        ##        for line in fileinput.input([transformFile], inplace=True):
        ##
        ##            if 'InitialTransformParametersFileName' in line :
        ##
        ##                pathLine = line;
        ##                line = line.replace(pathLine, '(InitialTransformParametersFileName "{}")'.format(os.path.join(Path2Alignement,'TransformParameters.0.txt')));
        ##
        ##            line = line.replace('(FinalBSplineInterpolationOrder 3)', '(FinalBSplineInterpolationOrder 0)');
        ##            sys.stdout.write(line);
        ##
        #        sinkDir = os.path.join(Parameters.workdir,"elastix_annotation_transformed");
        ##        sourceFile = annotationFile;
        ##
        ##        elx.transformData(sourceFile,sink=os.path.join(sinkDir,"annotation_transformed_template.tif"),\
        ##                          transformParameterFile=transformFile,resultDirectory=sinkDir);
        ##
        #        transformationResult = os.path.join(sinkDir,"annotation_transformed_template.tif");
        ##
        ##
        ##
        ##
        ##
        #        Path2Alignement = os.path.join(Parameters.workdir,'elastix_cfos_to_auto');
        ##        copyfile(os.path.join(Path2Alignement,'TransformParameters.0.txt'), os.path.join(Path2Alignement,'TransformParameters_noInterpolation.0.txt'));
        #        transformFile = os.path.join(Path2Alignement,'TransformParameters_noInterpolation.0.txt');
        ##
        ##        for line in fileinput.input([transformFile], inplace=True):
        ##
        ##            line = line.replace('(FinalBSplineInterpolationOrder 3)', '(FinalBSplineInterpolationOrder 0)');
        ##            sys.stdout.write(line);
        #
        #        sinkDir2 = os.path.join(Parameters.workdir,"elastix_annotation_transformed_to_raw");
        #
        #        elx.transformData(transformationResult,sink=os.path.join(sinkDir2,"annotation_transformed_template_to_cfos.tif"),\
        #                          transformParameterFile=transformFile,resultDirectory=sinkDir2);
        #
        #        result = os.path.join(sinkDir2,"annotation_transformed_template_to_cfos.tif");
        #        result = io.readData(result);
        ##        result = result + 32768;
        ##        result = np.uint16(result);
        #
        #        io.writeData(os.path.join(sinkDir2,"annotation_transformed_template_to_cfos_2.tif"), result);
        #
        ##        print("starting upsampling")
        ##
        ##        Parameters.UpsamplingParameterAnnotation["source"] = os.path.join(sinkDir2,"edges.tif");
        ##        print("reading data")
        ##        cFos = io.readData(os.path.join(Parameters.outdir, "Z0000.tif"));
        ##        Parameters.UpsamplingParameterAnnotation["dataSizeSink"] = (cFos.shape[0],cFos.shape[1],len(os.listdir(Parameters.outdir)));
        ##        Parameters.UpsamplingParameterAnnotation["resolutionSink"] = Parameters.DataResolution
        ##
        ##        print("upsampling data")
        ##        rs.resampleData(**Parameters.UpsamplingParameterAnnotation);

        ############################################################################################################################################
        ############################################################################################################################################
        # Annotation to Raw Data

        ResamplingParametercFos10m = Parameters.ResamplingParametercFos.copy();
        ResamplingParametercFos10m["sink"] = os.path.join(Parameters.workdir, 'cfos_resampled_10m.tif');
        ResamplingParametercFos10m["resolutionSink"] = Parameters.TemplateResolution10m;

        print("\n");
        print(ut.coloredMessage(ut.titleMessage('Starting resampling cFos'), 'darkgreen'));

        rs.resampleData(**ResamplingParametercFos10m);

        if kwargs["playSound"]:
            ut.playSound(1);

#        ######################################################################
#        
#        annotationFile = Parameters.AnnotationFile;
#          
#        Path2Alignement = os.path.join(Parameters.workdir,'elastix_template_to_auto');
#        copyfile(os.path.join(Path2Alignement,'TransformParameters.1.txt'), os.path.join(Path2Alignement,'TransformParameters_noInterpolation.1.txt'));
#        transformFile = os.path.join(Path2Alignement,'TransformParameters_noInterpolation.1.txt');
#        
#        for line in fileinput.input([transformFile], inplace=True):
#            
#            if 'InitialTransformParametersFileName' in line :
#                
#                pathLine = line;
#                line = line.replace(pathLine, '(InitialTransformParametersFileName "{}")'.format(os.path.join(Path2Alignement,'TransformParameters.0.txt')));
#                
#            line = line.replace('(FinalBSplineInterpolationOrder 3)', '(FinalBSplineInterpolationOrder 0)');
#            sys.stdout.write(line);
#          
#        sinkDir = os.path.join(Parameters.workdir,"elastix_annotation_transformed");
#        sourceFile = annotationFile;
#
#        elx.transformData(sourceFile,sink=os.path.join(sinkDir,"annotation_transformed_template.tif"),\
#                          transformParameterFile=transformFile,resultDirectory=sinkDir);
#                     
#        #######################################################################
#        
#        Path2Alignement = os.path.join(Parameters.workdir,'elastix_cfos_to_auto');
#        copyfile(os.path.join(Path2Alignement,'TransformParameters.0.txt'), os.path.join(Path2Alignement,'TransformParameters_noInterpolation.0.txt'));
#        transformFile = os.path.join(Path2Alignement,'TransformParameters_noInterpolation.0.txt');
#        
#        for line in fileinput.input([transformFile], inplace=True):
#            
#          line = line.replace('(FinalBSplineInterpolationOrder 3)', '(FinalBSplineInterpolationOrder 0)');
#          sys.stdout.write(line);
#          
#        sinkDir = os.path.join(Parameters.workdir,"elastix_annotation_transformed");
#        sourceDir = sinkDir;  
#        sourceFile = os.path.join(sourceDir, "annotation_transformed_template.tif");
#
#        elx.transformData(sourceFile,\
#                          sink=os.path.join(sinkDir,"annotation_transformed_cfos.tif"),transformParameterFile=transformFile,resultDirectory=sinkDir);
#
#        ######################################################################
#        
#        Parameters.UpsamplingParameterAnnotation["source"] = os.path.join(os.path.join(Parameters.workdir,"elastix_annotation_transformed"),"annotation_transformed_cfos.tif");
#        cFosResampled = io.readData(ResamplingParametercFos10m["sink"]);
#        Parameters.UpsamplingParameterAnnotation["dataSizeSink"] = (cFosResampled.shape[0],cFosResampled.shape[1],cFosResampled.shape[2]); #
#        
#        rs.resampleData(**Parameters.UpsamplingParameterAnnotation);

#        ############################################################################################################################################
#        ############################################################################################################################################
#        #EYFP pipeline 10m
#        
#        ######################################################################
#        #Resampling raw data
#        
#        if not os.path.exists(os.path.join(Parameters.workdir,"10m")) :
#            
#            os.mkdir(os.path.join(Parameters.workdir,"10m"));
#            
#        ResamplingParametercFos10m = Parameters.ResamplingParametercFos.copy();
#        ResamplingParametercFos10m["sink"] = os.path.join(os.path.join(Parameters.workdir,"10m"), 'raw_resampled_10m.tif');
#        ResamplingParametercFos10m["resolutionSink"] = Parameters.TemplateResolution10m;
#        
#        rs.resampleData(**ResamplingParametercFos10m);
#        
#        ######################################################################
#        #Resampling autofluo
#        
#        ResamplingParameterAutoFluo10m = Parameters.ResamplingParameterAutoFluo.copy();
#        ResamplingParameterAutoFluo10m["sink"] = os.path.join(os.path.join(Parameters.workdir,"10m"), 'auto_resampled_10m.tif');
#        ResamplingParameterAutoFluo10m["resolutionSink"] = Parameters.TemplateResolution10m;
#        
#        rs.resampleData(**ResamplingParameterAutoFluo10m);
#        
#        ######################################################################
#        #Align auto to cfos
#        
#        ChannelsAlignmentParameterBis = Parameters.ChannelsAlignmentParameter.copy();
#        ChannelsAlignmentParameterBis["movingImage"], ChannelsAlignmentParameterBis["fixedImage"] = ResamplingParametercFos10m["sink"], ResamplingParameterAutoFluo10m["sink"]
#        ChannelsAlignmentParameterBis["resultDirectory"] = os.path.join(os.path.join(Parameters.workdir,"10m"), 'elastix_auto_to_cfos')
#        
#        resultDirectory = elx.alignData(**ChannelsAlignmentParameterBis);
#        
#        resultFile = os.path.join(ChannelsAlignmentParameterBis["resultDirectory"], "result.0.mhd");
#        resultData = io.readData(resultFile);
#        io.writeData( os.path.join( os.path.dirname(resultFile) , "result.0.tif"), resultData);
#        
#        ######################################################################
#        #Align auto to template
#        
#        TemplateAlignmentParameterBis = Parameters.TemplateAlignmentParameter.copy();
#        TemplateAlignmentParameterBis["movingImage"], TemplateAlignmentParameterBis["fixedImage"] = ResamplingParameterAutoFluo10m["sink"], Parameters.TemplateFile10m
#        TemplateAlignmentParameterBis["resultDirectory"] = os.path.join(os.path.join(Parameters.workdir,"10m"), 'elastix_auto_to_template')
#        
#        resultDirectory  = elx.alignData(**TemplateAlignmentParameterBis);
#        
#        ######################################################################
#        #Tranform raw data to template
#        
#        Path2Alignement = os.path.join(os.path.join(Parameters.workdir,"10m"),'elastix_auto_to_template');
#        copyfile(os.path.join(Path2Alignement,'TransformParameters.1.txt'), os.path.join(Path2Alignement,'TransformParameters_noInterpolation.1.txt'));
#        transformFile = os.path.join(Path2Alignement,'TransformParameters_noInterpolation.1.txt');
#        
#        for line in fileinput.input([transformFile], inplace=True):
#            
#            if 'InitialTransformParametersFileName' in line :
#                
#                pathLine = line;
#                line = line.replace(pathLine, '(InitialTransformParametersFileName "{}")'.format(os.path.join(Path2Alignement,'TransformParameters.0.txt')));
#                
#            line = line.replace('(FinalBSplineInterpolationOrder 3)', '(FinalBSplineInterpolationOrder 0)');
#            sys.stdout.write(line);
#          
#        sinkDir = os.path.join(os.path.join(Parameters.workdir,"10m"),"elastix_raw_transformed");
#        source = os.path.join( os.path.dirname(resultFile) , "result.0.tif");
#
#        elx.transformData(source,sink=os.path.join(sinkDir,"raw_transformed.tif"),\
#                          transformParameterFile=transformFile,resultDirectory=sinkDir);

############################################################################################################################################
############################################################################################################################################
# EYFP pipeline 25m

######################################################################
# Align auto to cfos

#        ChannelsAlignmentParameterBis = Parameters.ChannelsAlignmentParameter.copy();
#        ChannelsAlignmentParameterBis["movingImage"], ChannelsAlignmentParameterBis["fixedImage"] = ResamplingParametercFos["sink"], ResamplingParameterAutoFluo["sink"]
#        ChannelsAlignmentParameterBis["resultDirectory"] = os.path.join(Parameters.workdir, 'elastix_auto_to_cfos')
#        
#        resultDirectory = elx.alignData(**ChannelsAlignmentParameterBis);
#        
#        resultFile = os.path.join(ChannelsAlignmentParameterBis["resultDirectory"], "result.0.mhd");
#        resultData = io.readData(resultFile);
#        io.writeData( os.path.join( os.path.dirname(resultFile) , "result.0.tif"), resultData);
#        
#        ######################################################################
#        #Align auto to template
#        
#        TemplateAlignmentParameterBis = Parameters.TemplateAlignmentParameter.copy();
#        TemplateAlignmentParameterBis["movingImage"], TemplateAlignmentParameterBis["fixedImage"] = ResamplingParameterAutoFluo["sink"], Parameters.TemplateFile
#        TemplateAlignmentParameterBis["resultDirectory"] = os.path.join(Parameters.workdir, 'elastix_auto_to_template')
#        
#        resultDirectory  = elx.alignData(**TemplateAlignmentParameterBis);
#        
#        ######################################################################
#        #Tranform raw data to template
#        
#        Path2Alignement = os.path.join(Parameters.workdir,'elastix_auto_to_template');
#        copyfile(os.path.join(Path2Alignement,'TransformParameters.1.txt'), os.path.join(Path2Alignement,'TransformParameters_noInterpolation.1.txt'));
#        transformFile = os.path.join(Path2Alignement,'TransformParameters_noInterpolation.1.txt');
#        
#        for line in fileinput.input([transformFile], inplace=True):
#            
#            if 'InitialTransformParametersFileName' in line :
#                
#                pathLine = line;
#                line = line.replace(pathLine, '(InitialTransformParametersFileName "{}")'.format(os.path.join(Path2Alignement,'TransformParameters.0.txt')));
#                
#            line = line.replace('(FinalBSplineInterpolationOrder 3)', '(FinalBSplineInterpolationOrder 0)');
#            sys.stdout.write(line);
#          
#        sinkDir = os.path.join(Parameters.workdir,"elastix_raw_transformed");
#        source = os.path.join( os.path.dirname(resultFile) , "result.0.tif");
#
#        elx.transformData(source,sink=os.path.join(sinkDir,"raw_transformed.tif"),\
#                          transformParameterFile=transformFile,resultDirectory=sinkDir);
#    
#    print("\n");
#    print(ut.coloredMessage(ut.titleMessage('Operations on Brain %s finished'),'darkgreen'));
#    
#    if kwargs["playSound"] :
#          ut.playSound(1);
