#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 20 10:08:29 2018

@author: thomas.topilko
"""

import os;
import sys;
from natsort import natsorted;
import numpy as np;

import ClearMap.Analysis.Statistics as stat;
import ClearMap.Analysis.Label as lbl;
import ClearMap.IO.IO as io;
import ClearMap.Alignment.Resampling as rs;
import ClearMap.Analysis.Tools.MultipleComparisonCorrection as FDR;

##############################################################################
# Miscellaneous
##############################################################################

cfosPrefix = ['cfos','cFos','CFOS','CFos','fos','FOS','Fos','ofcs','cofs','igg','eyfp','gfp','camKII', "projections"];
autoPrefix = ['auto','AUTO',"testEdU-lowmag"];

# Parameters for message outputs

colors = {
    'white':    "\033[1;37m",
    'yellow':   "\033[1;33m",
    'green':    "\033[1;32m",
    'blue':     "\033[1;34m",
    'cyan':     "\033[1;36m",
    'red':      "\033[1;31m",
    'magenta':  "\033[1;35m",
    'black':      "\033[1;30m",
    'darkwhite':  "\033[0;37m",
    'darkyellow': "\033[0;33m",
    'darkgreen':  "\033[0;32m",
    'darkblue':   "\033[0;34m",
    'darkcyan':   "\033[0;36m",
    'darkred':    "\033[0;31m",
    'darkmagenta':"\033[0;35m",
    'darkblack':  "\033[0;30m",
    'bold' :      "\033[1m",
    'off':        "\033[0;0m"
};

##############################################################################
### Utilities
##############################################################################

def playSound(n) :

  for i in range(n) :
    os.system("/usr/bin/canberra-gtk-play --id='bell'")

def coloredMessage(source,color) :
    
    ''' Function that makes a text message colored
    
    Arguments:
        source (str) = the text to be colored
        color (str) = the bash color code
    '''
    
    return colors[color]+source+colors['off'];
  
def tryObject(var) :
  
  ''' Function that check if a variable exists, if not creates it as None
    
  Arguments:
      var (variable) = object/variable to be tested
  '''
  
  try :
    
    var;
    
    #return var;
    
  except :
    
    print(coloredMessage('[/!\ Warning] {} is not defined !'.format(var),'darkred'));
    var = None;
    
    #return var;
    
def checkPathExists (path) :
  
  path2File = os.path.split(path)[0];
  f = path.split("_")[-1];
  fName = f.split(".")[0];
  
  if os.path.exists(path) :
    print(coloredMessage('[INFO] {0} file has been successfully detected'.format(fName),'darkgreen'));
    
  else :
    raise RuntimeError ('[INFO] No {0} file detected in {1}!'.format(fName,path2File));

def loadingBar(percent) :
    
    ''' Function that creates a loading bar
    
    Arguments:
        percent (int) = percentage of progression (example : 50 = 50%)
    
    Returns:
        pattern (list) = loading bar as a list
    '''
    
    pattern = ['.','.','.','.','.','.','.','.','.','.','.','.','.','.','.','.',\
               '.','.','.','.','.','.','.','.','.','.','.','.','.','.','.','.',\
               '.','.','.','.','.','.','.','.','.','.','.','.','.','.','.','.',\
               '.','.'];
    
    for n,i in enumerate(pattern) : 
        
        if n < percent/2 :
            pattern[n] = '#';
                   
        else :
            pattern[n] = '.';
                   
    pattern = "".join(pattern);
    
    return pattern;

def hoursMinutesSeconds(time) :
    
    ''' Function that splits hours minutes and seconds from an large input of seconds
    
    Arguments:
        time (int) = seconds
    
    Returns:
        hours (int) = time in hours
        minutes (int) = time in minutes
        seconds (int) = time in seconds
    '''
    
    remainingMinutes = time%3600;
    remainingSeconds = remainingMinutes%60;
    hours = int(time/3600);
    minutes = int(remainingMinutes/60);
    seconds = int(remainingSeconds);
    
    return hours,minutes,seconds;

def askBeforeExecuting(toolTip,verbose=True) : 
    
    ''' Function that serves as a check before executing a command
    
    Arguments:
        toolTip (str) = Description of the function
        verbose (bol) = True : outputs text in console, False : outputs in console are disabled
        sound (bol) = True : plays sound, False : sounds are disabled
        
    Returns:
        True/False (bol) = True : executes the function, False : doesn't
    '''
                          
    if sys.version_info[0] == 2 :
        sys.stdin.flush();
        Input = raw_input(coloredMessage("[INFO] Would you like to "+toolTip+" ? : please type"\
                                     +" 'Yes', 'No' to proceed : ","darkred"));
    
    if sys.version_info[0] == 3 :
        Input = input(coloredMessage("[INFO] Would you like to "+toolTip+" ? : please type"\
                                     +" 'Yes', 'No' to proceed : ","darkred"));
                                   
    if Input == 'Yes' or Input == 'yes' :       
        return True; 
    
    elif Input == 'No' or Input == 'no' :  
        
        if verbose == True :
            print(coloredMessage("[WARNING] The function was not executed \n"\
                                   ,"darkred"));
        return False;
    
    else :
        
        if verbose == True :
            print(coloredMessage("[WARNING] Argument "+Input+" unknown, the function was not"\
                                 +" executed \n","darkred"));
        return False;
    
def askForOperations(Mouse,verbose=True,sound=True) : 
    
    ''' Function that serves as a check before executing a command
    
    Arguments:
        toolTip (str) = Description of the function
        verbose (bol) = True : outputs text in console, False : outputs in console are disabled
        sound (bol) = True : plays sound, False : sounds are disabled
        
    Returns:
        True/False (bol) = True : executes the function, False : doesn't
    '''
                          
    if sys.version_info[0] == 2 :
        sys.stdin.flush();
        Input = raw_input(coloredMessage("[INFO] Operations on brain %s : "%(Mouse),"darkgreen"));
    
    if sys.version_info[0] == 3 :
        Input = input(coloredMessage("[INFO] Operations on brain %s : "%(Mouse),"darkgreen"));
    
    if Input == 'all' or Input == 'All' :
        Input = 'TTTTTT';
    elif Input == 'none' or Input == 'None' :  
        Input = 'FFFFFF'
                   
    return Input
  
def titleMessage(msg) :

  message = '[INFO] ########## '+str(msg)\
            +' ##########'
  
  return message;

def printToolTip() :
  
  print(coloredMessage("[INFO] The steps that can be run for each brain are : "\
                      +"Stiching, Resampling auto, Resampling cFos, Alignement_FOS_Auto"\
                      +",Alignement_Temp_Auto, Cell Detection. \n"\
                      +"The format to run those steps is : (T..F) T=True F=False, all = All True "\
                      "none = All False \n","darkgreen"));
                       
                       
###########################################################################
### Function that lauches the Manual/GUI modes
##############################################################################

def launchManual(fd,toolTip) :
  
  print(coloredMessage('[INFO] Manual mode selected, GUI_Clearmap not used','darkred'));
  
  operations = {};

  for folder in natsorted(os.listdir(fd)) :
      if os.path.isdir(os.path.join(folder,fd)) :
          if toolTip :
              printToolTip();
              toolTip = False;
          sequence = askForOperations(folder);
          operations[folder] = sequence;
          
  return operations;

def setPaths(fd,experiment,operations,prefix) :
  
  all_false = 'FFFFFFFF';
  
  for folder in natsorted(os.listdir(fd)) :
  
    auto = False;
    cfos = False;
    
    path2Folder = os.path.join(fd,folder);
    
    if os.path.isdir(path2Folder) : 
      
        if folder.split('_')[0] != folder :
          exp = folder.split('_')[0];
        elif folder.split('-')[0] != folder :
          exp = folder.split('-')[0];
        else :
          raise RuntimeError('No folder(s) with the name {} were found in {} \n'.format(experiment,fd));
          
        if exp == experiment : 
          
          if operations[folder] != all_false :
          
            workdir = path2Folder;
            print(coloredMessage('[INFO] Current Workdir is : %s'%(workdir),'darkgreen'));
            
            for subfolder in os.listdir(os.path.join(fd,folder)) :
                
                path2Subfolder = os.path.join(path2Folder,subfolder);
                
                if os.path.isdir(path2Subfolder) :
                  
                    timeStamp = subfolder.split('_')[-1];
                    
                    if len(timeStamp.split('-')) == 3 :
                    
                      try :
                        dataType = subfolder.split('_')[1];
                        #print(dataType);
                      except :
                        dataType = None;
                      
                      if dataType in prefix :
                          datadir = path2Subfolder;
                          nCFosFiles = len(os.listdir(datadir));   
                          print(coloredMessage('[INFO] Current Datadir is : %s'%(datadir),'darkgreen'));
                          cFosfile = os.path.join(datadir, r''+timeStamp+'_%s_UltraII\[(?P<row>\d{2}) x (?P<col>\d{2})\]_C00.ome.tif'%(dataType))
                          print(coloredMessage('[INFO] %d %s files have been detected'%(nCFosFiles,dataType),'darkgreen'));
                          cfos = True;
                      elif dataType in autoPrefix :
                          autodir = path2Subfolder;
                          nAutoFiles = len(os.listdir(autodir));   
                          print(coloredMessage('[INFO] Current Autodir is : %s'%(autodir),'darkgreen'));
                          AutoFluoFile = os.path.join(autodir, timeStamp+'_%s_UltraII_C00_xyz-Table Z\d{4}.ome.tif'%(dataType));
                          print(coloredMessage('[INFO] %d %s files have been detected \n'%(nAutoFiles,dataType),'darkgreen'));
                          auto = True;
                          
    return workdir,datadir,cFosfile,cfos,autodir,AutoFluoFile,auto;

def getSampleNames (directory, experiment, sep) :

    names = [];
    
    for folder in natsorted(os.listdir(directory)) :
        
      path2folder = os.path.join(directory,folder);
      if os.path.isdir(path2folder) :
        try :
          exp = folder.split(sep)[0];
        except :
          exp = None;
        if exp == experiment :
          names.append(sep.join(folder.split(sep)[1:]));
          
    return names;

class Group() :
  
  def __init__ (self, **kwargs) :
      
    self.name = kwargs["group"];
    
    if kwargs["analysisDir"] == None : 
        self.paths = [kwargs["homedir"]+os.sep+kwargs["experiment"]+kwargs["sep"]\
                      +str(y)+os.sep+'cells_heatmap.tif' for x,y in\
                      zip(kwargs["operation"], kwargs["sampleNames"]) if x == 'T'];
                      
        self.mouseNames = [x.split("/")[-2].split(kwargs["sep"])[-1] for x\
                       in self.paths];
        
    else :
        self.paths = [os.path.join(kwargs["savingDir"], 'cells_heatmap_{0}_{1}.tif'.\
                                   format(kwargs["experiment"], y)) for x,y in\
                                    zip(kwargs["operation"],kwargs["sampleNames"]) if x == 'T'];
                                   
        self.mouseNames = [x.split("/")[-1].split(kwargs["sep"])[-1].split(".")[-2] for x\
                       in self.paths];
        
    self.data = stat.readDataGroup(self.paths);
    self.mean = np.mean(self.data,axis = 0);

def combinliste(seq, k):
    p = []
    i, imax = 0, 2**len(seq)-1
    while i<=imax:
        s = []
        j, jmax = 0, len(seq)-1
        while j<=jmax:
            if (i>>j)&1==1:
                s.append(seq[j])
            j += 1
        if len(s)==k:
            p.append(s)
        i += 1 
    return p

def combinlisterep(seq, k):
    """Renvoie la liste des combinaisons avec répétition des objets de seq pris k à k"""
    # ajoute chaque objet de seq pour quils apparaissent chacun k fois
    seq2 = []
    for elem in seq:
        if elem not in seq2:
            for i in range(0,k):
                seq2.append(elem)
    # calcule la liste "normale" des combinaisons
    p = combinliste(seq2, k)
    # élimine de cette liste les éléments identiques (comme [1,2] et [1,2])
    p2 = []
    for x in p:
        if x not in p2 and len(set(x)) > 1:
            p2.append(x)
    # et renvoie le résultat
    return p2


def launchPairAnalysis(groups, combinations, sinkDir, annotationFile, cutoff=0.05) :
    
    groupNames = list(groups.keys());
    
    for c in combinations :
        
        print("\n");
        print(coloredMessage(titleMessage('Generating p-value maps {0} vs {1}'.format(groupNames[c[0]],\
                                                groupNames[c[1]])),'darkgreen'));
        
        pvals, psign = stat.tTestVoxelization(groups[groupNames[c[0]]].data.astype('float'),\
                                              groups[groupNames[c[1]]].data.astype('float'),\
                                              signed = True, pcutoff = cutoff);
                                              
        pvalsc = stat.colorPValues(pvals, psign, positive = [0,1], negative = [1,0]);
        
        io.writeData(os.path.join(sinkDir, 'pvalues_'+groups[groupNames[c[0]]].name+'_vs_'+groups[groupNames[c[1]]].name+'.tif'),\
                     rs.sagittalToCoronalData(pvalsc.astype('float32')));
        
        print("\n");
        print(coloredMessage(titleMessage('Starting region based stats {0} vs {1}'.format(groups[groupNames[c[0]]].name,\
                                                groups[groupNames[c[1]]].name)),'darkgreen'));
                                          
        ids, pc1, pci1 = stat.countPointsGroupInRegions(groups[groupNames[c[0]]].paths2, intensityGroup = groups[groupNames[c[0]]].i,\
                                                      returnIds = True, labeledImage = annotationFile, returnCounts = True, collapse=True);
        _, pc2, pci2 = stat.countPointsGroupInRegions(groups[groupNames[c[1]]].paths2, intensityGroup = groups[groupNames[c[1]]].i,\
                                              returnIds = True, labeledImage = annotationFile, returnCounts = True, collapse=True);
                                                      
        pvals, psign = stat.tTestPointsInRegions(pc1, pc2, pcutoff = None, signed = True, equal_var = True);
        pvalsi, psigni = stat.tTestPointsInRegions(pci1, pci2, pcutoff = None, signed = True, equal_var = True);
        
        id1 = pvals < 1;
        Ids = ids[id1];
        
        Pc1 = pc1[id1];
        Pc2 = pc2[id1];
        
        Psign = psign[id1];
        Pvals = pvals[id1];
        Qvals = FDR.estimateQValues(Pvals);
        
        dtypes = [('id','int64'),('mean1','f8'),('std1','f8'),('mean2','f8'),('std2','f8'),('pvalue', 'f8'),('qvalue', 'f8'),('psign', 'int64')];

        for i in groups[groupNames[c[0]]].mouseNames:
            dtypes.append(("{0}".format(i), 'f8'));
        for i in groups[groupNames[c[1]]].mouseNames:
            dtypes.append(("{0}".format(i), 'f8'));   
        dtypes.append(('name', 'U100'));
        
        cellTable = np.zeros(Ids.shape, dtype = dtypes)
        cellTable["id"] = Ids;
        cellTable["mean1"] = Pc1.mean(axis = 1)
        cellTable["std1"] = Pc1.std(axis = 1)
        cellTable["mean2"] = Pc2.mean(axis = 1)
        cellTable["std2"] = Pc2.std(axis = 1)
        cellTable["pvalue"] = Pvals;
        cellTable["qvalue"] = Qvals;
        
        cellTable["psign"] = Psign;
        for n, i in enumerate(groups[groupNames[c[0]]].mouseNames):
            cellTable["{0}".format(i)] = Pc1[:,n];
        for n, i in enumerate(groups[groupNames[c[1]]].mouseNames):
            cellTable["{0}".format(i)] = Pc2[:,n];
        cellTable["name"] = [("").join(x.split(",")) for x in lbl.labelToName(Ids)];
        
        ii = np.argsort(Pvals);
        cellTableSorted = cellTable.copy();
        cellTableSorted = cellTableSorted[ii];
        
        with open((sinkDir+os.sep+'counts-cells_table_'+groups[groupNames[c[0]]].name+'_vs_'+groups[groupNames[c[1]]].name+'.csv'),'w') as f:
            f.write(', '.join([str(item) for item in cellTable.dtype.names]));
            f.write('\n');
            for sublist in cellTableSorted:
                f.write(', '.join([str(item) for item in sublist]));
                f.write('\n');
            f.close();
            
        Pci1 = pci1[id1];
        Pci2 = pci2[id1];
        
        Psigni = psigni[id1];
        Pvalsi = pvalsi[id1];
        Qvalsi = FDR.estimateQValues(Pvalsi);
        
        intensityTable = np.zeros(Ids.shape, dtype = dtypes)
        intensityTable["id"] = Ids;
        intensityTable["mean1"] = Pci1.mean(axis = 1)
        intensityTable["std1"] = Pci1.std(axis = 1)
        intensityTable["mean2"] = Pci2.mean(axis = 1)
        intensityTable["std2"] = Pci2.std(axis = 1)
        intensityTable["pvalue"] = Pvalsi;
        intensityTable["qvalue"] = Qvalsi;
        
        intensityTable["psign"] = Psigni;
        for n, i in enumerate(groups[groupNames[c[0]]].mouseNames) :
            intensityTable["{0}".format(i)] = Pci1[:,n];
        for n, i in enumerate(groups[groupNames[c[1]]].mouseNames) :
            intensityTable["{0}".format(i)] = Pci2[:,n];
        intensityTable["name"] = [("").join(x.split(",")) for x in lbl.labelToName(Ids)];
        
        ii = np.argsort(Pvalsi);
        intensityTableSorted = intensityTable.copy();
        intensityTableSorted = intensityTableSorted[ii];
        
        with open((sinkDir+os.sep+'counts-intensities_table_'+groups[groupNames[c[0]]].name+'_vs_'+groups[groupNames[c[1]]].name+'.csv'),'w') as f:
            f.write(', '.join([str(item) for item in intensityTable.dtype.names]));
            f.write('\n');
            for sublist in intensityTableSorted:
                f.write(', '.join([str(item) for item in sublist]));
                f.write('\n');
            f.close();
  
  