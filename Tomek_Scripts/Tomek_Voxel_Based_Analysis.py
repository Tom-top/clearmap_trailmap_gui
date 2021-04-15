#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 15:16:26 2019

@author: thomas.topilko
"""

import os;
import numpy as np;
import pandas as pd;
import itertools;
from scipy import stats;


import Parameters;
import ClearMap.IO as io;
from ClearMap.Analysis.Voxelization import voxelize;
import ClearMap.Alignment.Resampling as rs;
import ClearMap.Analysis.Statistics as stat;

class Class() :
    
    def __init__(self) :
        
        self.groups = {
                        "Pups" : [1,3,5,8], 
                        "Rubbers" : [2,4,6,9,10],
                        };
                
        self.dataPath = "/home/thomas.topilko/Desktop/tSNE";
        self.dataClearMapPath = os.path.join(self.dataPath, "ClearMap");
        self.dataBehaviorPath = os.path.join(self.dataPath, "Behavior");
        
        self.savingPath = os.path.join(self.dataPath, "Figures");
        
        self.PupsVoxelMatrix = [];
        self.RubbersVoxelMatrix = [];
        
        #Looping over groups
        for k, v in self.groups.items() :
            
            Buffer = [];
            
            #Looping over animals in group
            for animal in v :
                
                #Looping over folders in ClearMap data
                for f in os.listdir(self.dataClearMapPath) :
                    
                    if int(f.split("_")[0]) == int(animal) :

                        dataFile = os.path.join(self.dataClearMapPath, f);
                        voxelizedData = self.Voxelize(dataFile);
                
                Buffer.append(voxelizedData);
                
            if k == "Pups" :
                
                self.PupsVoxelMatrix.append(Buffer);
            
            elif k == "Rubbers" :
                
                self.RubbersVoxelMatrix.append(Buffer);
                
        
        self.PupsBehaviorMatrix = [];
        self.RubbersBehaviorMatrix = [];
        
        #Looping over groups
        for k, v in self.groups.items() :
            
            Buffer = [];
            
            #Looping over animals in group
            for animal in v :
                
                #Looping over folders in ClearMap data
                for f in os.listdir(self.dataBehaviorPath) :
                    
                    if f.split(".")[-1] == "xlsx" :
                        
                        if int(f.split(".")[0]) == int(animal) :
    
                            dataFile = os.path.join(self.dataBehaviorPath, f);
                            excelData = pd.read_excel(dataFile, header=None);
                            
                            Start = [i for i in excelData.loc[4][1:] if i > 0];
                            End = [i for i in excelData.loc[5][1:] if i > 0];
                            
                            Delta = [e-s for s,e in zip(Start, End)];
                            Delta = np.sum(Delta)
                        
                Buffer.append(Delta);
                
            if k == "Pups" :
                
                self.PupsBehaviorMatrix.append(Buffer);
            
            elif k == "Rubbers" :
                
                self.RubbersBehaviorMatrix.append(Buffer);
                  
    def Voxelize(self, File) :
        
        points = io.readPoints(File);
        vox = voxelize(points, "/home/thomas.topilko/Documents/ClearMap_GUI/ClearMap_Ressources/25um_Autofluo_Reference/Template_1-279.tif", **Parameters.voxelizeParameter); 
        
        return vox;
      
    def GenerateVoxelMatrix(self, mat) :
        
        mat = np.array(mat[0], dtype=float);
        minMat = np.min(mat, axis=0, keepdims=1);
        maxMat = np.max(mat, axis=0, keepdims=1);
        
        sink = (mat - minMat)/ (maxMat - minMat);
        
        pi = np.isnan(sink);
        sink[pi] = 0.;
        
        return sink;
            
    def GeneratetSNEMatrixBehavior(self, Matrix) :
        
        Matrix = Matrix[0];
        
        minMatrix = float(np.array(Matrix).min());
        maxMatrix = float(np.array(Matrix).max());
            
        normalizedMatrix = np.array([ (n - minMatrix) / (maxMatrix - minMatrix) if (maxMatrix - minMatrix) != 0 else 0 for n in np.array(Matrix) ]);
        
        return normalizedMatrix;
        
         
#Initialize the vba class 
behavior = "pup_interaction";
vba = Class();

#Normalize voxel data
normalizedVoxelMatrix = vba.GenerateVoxelMatrix(vba.PupsVoxelMatrix);

#Generate the behavior matrix
rawBehaviorMatrix = vba.GeneratetSNEMatrixBehavior(vba.PupsBehaviorMatrix);

#Reshape behavior matrix to fit the voxel matrix
shapedBehaviorMatrix = np.array([ np.full_like(np.array(normalizedVoxelMatrix[0,:,:,:],dtype=float),\
                                               float(normalizedVoxelMatrix[i])) for i in np.arange(len(rawBehaviorMatrix)) ]);

pcoeffs = np.full_like(normalizedVoxelMatrix[0], 0.);

for i,j,k in itertools.product(np.arange(normalizedVoxelMatrix.shape[1]),\
                               np.arange(normalizedVoxelMatrix.shape[2]),\
                               np.arange(normalizedVoxelMatrix.shape[3])):
    
    if np.isnan(stats.pearsonr(normalizedVoxelMatrix[:,i,j,k], normalizedVoxelMatrix[:,i,j,k])[0]) :
                pcoeffs[i,j,k] = 0.;
    else :
        pcoeffs[i,j,k] = stats.pearsonr(normalizedVoxelMatrix[:,i,j,k], normalizedVoxelMatrix[:,i,j,k])[0];

pcoeffsc = stat.colorPearsonsCoeffs(pcoeffs, cutoff=0.1);

io.writeData(os.path.join("/home/thomas.topilko/Desktop", 'pvalues_correlation_{0}.tif'.format(behavior)), rs.sagittalToCoronalData(pcoeffsc.astype('float32')));
np.save(os.path.join("/home/thomas.topilko/Desktop", 'pcoeffs_correlation_{0}.npy'.format(behavior)), pcoeffs)
