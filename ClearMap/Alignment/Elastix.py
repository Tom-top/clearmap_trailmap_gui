# -*- coding: utf-8 -*-
"""
Interface to Elastix for alignment of volumetric data

The elastix documentation can be found `here <http://elastix.isi.uu.nl/>`_.

In essence, a transformation :math:`T(x)` is sought so that for a fixed image 
:math:`F(x)` and a moving image :math:`M(x)`:

.. math::
    F(x) = M(T(x))

Once the map :math:`T` is estimated via elastix, transformix maps an image
:math:`I(x)` from the moving image frame to the fixed image frame, i.e.:

.. math::
    I(x) \\rightarrow I(T(x)) 

To register an image onto a reference image, the fixed image is typically 
choosed to be the image to be registered, while the moving image is the 
reference image. In this way an object identified in the data at position x
is mapped via transformix as:

.. math::
    x \\rightarrow T(x)
    



Summary
-------
    * elastix finds a transformation :math:`T: \\mathrm{fixed image} \\rightarrow \\mathrm{moving image}`
    * the fixed image is image to be registered
    * the moving image is typically the reference image
    * the result folder may contain an image (mhd file) that is :math:`T^{-1}(\\mathrm{moving})`,
      i.e. has the size of the fixed image
    * transformix applied to data gives :math:`T^{-1}(\\mathrm{data})` !
    * transformix applied to points gives :math:`T(\\mathrm{points})` !
    * point arrays are assumed to be in (x,y,z) coordinates consistent with (x,y,z) array represenation of images in ClearMap
    
Main routines are: :func:`alignData`, :func:`transformData` and :func:`transformPoints`.
    
See Also:
    `Elastix documentation <http://elastix.isi.uu.nl/>`_
    :mod:`~ClearMap.Alignment.Resampling`
"""
#:copyright: Copyright 2015 by Christoph Kirst, The Rockefeller University, New York City
#:license: GNU, see LICENSE.txt for details.


import os
import tempfile
import shutil
import re

import numpy as np

import multiprocessing as mp


import ClearMap.Settings as settings
import ClearMap.IO as io

##############################################################################
# Initialization and Enviroment Settings
##############################################################################

ElastixBinary = None;
"""str: the elastix executable

Notes:
    - setup in :func:`initializeElastix`
"""

ElastixLib = None;
"""str: path to the elastix library

Notes:
    - setup in :func:`initializeElastix`
"""

TransformixBinary = None;
"""str: the transformix executable

Notes:
    - setup in :func:`initializeElastix`
"""
    
Initialized = False;
"""bool: True if the elastixs binarys and paths are setup

Notes:
    - setup in :func:`initializeElastix`
"""

    
def printSettings():
    """Prints the current elastix configuration
    
    See also:
        :const:`ElastixBinary`, :const:`ElastixLib`, :const:`TransformixBinary`, :const:`Initialized`
    """
    
    global ElastixBinary, ElastixLib, TransformixBinary, Initialized
    
    if Initialized:
        print ("ElastixBinary     = %s" % ElastixBinary);
        print ("ElastixLib        = %s" % ElastixLib);
        print ("TransformixBinary = %s" % TransformixBinary);
    else:
        print ("Elastix not initialized");


def setElastixLibraryPath(path = None): 
    """Add elastix library path to the LD_LIBRARY_PATH variable in linux
    
    Arguments:
        path (str or None): path to elastix root directory if None :const:`ClearMap.Settings.ElastixPath` is used.
    """
     
    if path is None:
        path = os.path.join(settings.ElastixPath, "lib")

    if 'LD_LIBRARY_PATH' in os.environ :
    #if os.environ.has_key('LD_LIBRARY_PATH'):
        lp = os.environ['LD_LIBRARY_PATH'];
        if not path in lp.split(':'):
            os.environ['LD_LIBRARY_PATH'] = lp + ':' + path;
    else:
        os.environ['LD_LIBRARY_PATH'] = path


def initializeElastix(path = None):
    """Initialize all paths and binaries of elastix

    Arguments:
        path (str or None): path to elastix root directory, if None 
        :const:`ClearMap.Settings.ElastixPath` is used.
        
    See also:
        :const:`ElastixBinary`, :const:`ElastixLib`, :const:`TransformixBinary`,
        :const:`Initialized`, :func:`setElastixLibraryPath`
    """
    
    global ElastixBinary, ElastixLib, TransformixBinary, Initialized
    
    if path is None:
        path = settings.ElastixPath;
    
    #search for elastix binary
    elastixbin = os.path.join(path, 'bin/elastix');
    if os.path.exists(elastixbin):
        ElastixBinary = elastixbin;
    else:
        raise RuntimeError("Cannot find elastix binary %s, set path in Settings.py accordingly!" % elastixbin);
    
    #search for transformix binarx
    transformixbin = os.path.join(path, 'bin/transformix');
    if os.path.exists(transformixbin):
      TransformixBinary = transformixbin;
    else:
      raise RuntimeError("Cannot find transformix binary %s set path in Settings.py accordingly!" % transformixbin);
    
    #search for elastix libs
    elastixlib = os.path.join(path, 'lib');
    if os.path.exists(elastixlib):
      ElastixLib = elastixlib;
    else:
      elastixlib = os.path.join(path, 'bin');
      if os.path.exists(elastixlib):
        ElastixLib = elastixlib;
      else:
        raise RuntimeError("Cannot find elastix libs in %s  set path in Settings.py accordingly!" % elastixlib);
    
    #set path
    setElastixLibraryPath(ElastixLib);
        
    Initialized = True;
    
    print ("Elastix sucessfully initialized from path: %s" % path);
    
    return path;



initializeElastix();


def checkElastixInitialized():
    """Checks if elastix is initialized
    
    Returns:
        bool: True if elastix paths are set.
    """
    
    global Initialized;
    
    if not Initialized:
        raise RuntimeError("Elastix not initialized: run initializeElastix(path) with proper path to elastix first");
    #print ElastixSettings.ElastixBinary;

    return True;



##############################################################################
# Basic interface routines
####################################txtfile##########################################

def transformFile(resultdir):
    """Finds and returns the transformation parameter file generated by elastix
    
    Notes:
        In case of multiple transformation parameter files the top level file is returned     
    
    Arguments:
        resultdir (str): path to directory of elastix results
        
    Returns:
        str: file name of the first transformation parameter file 
    """    
    
    files = os.listdir(resultdir);
    files = [x for x in files if re.match('TransformParameters.\d.txt', x)];
    files.sort();
    
    if files == []:
        raise RuntimeError('Cannot find a valid transformation file in ' + resultdir);
    
    return os.path.join(resultdir, files[-1])


def transformDirectoryAndFile(transformParameterFile = None, transformDirectory = None):  
  """Determines transformation directory and file from either
     
  Arguments:
      transformParameterFile (str or None): file name of the transformation parameter file
      transformDirectory (str or None): directory to the transformation parameter
      
  Notes:
      Only one of the two arguments need to be specified.
  """
  if transformParameterFile == None:
    if transformDirectory == None:
      ValueError('Neither the alignment directory nor the transformation parameter file is specified!'); 
    transformparameterdir  = transformDirectory
    transformparameterfile = transformFile(transformparameterdir);
  else:
    transformparameterdir  = os.path.split(transformParameterFile);
    transformparameterdir  = transformparameterdir[0];
    transformparameterfile = transformParameterFile;
    
  return transformparameterdir, transformparameterfile;


def setPathTransformFiles(resultdir):
    """Replaces relative with absolute path in the parameter files in the result directory

    Notes:
        When elastix is not run in the directory of the transformation files
        the aboslute path needs to be given in each transformation file 
        to point to the subsequent transformation files. This is done via this 
        routine
       
    Arguments:
        resultdir (str): path to directory of elastix results
    """

    #print resultdir 
    files = os.listdir(resultdir);
    files = [x for x in files if re.match('TransformParameters.\d.txt', x)];
    files.sort();
    
    if files == []:
      raise RuntimeError('Cannot find a valid transformation file in ' + resultdir);
    
    rec = re.compile("\(InitialTransformParametersFileName \"(?P<parname>.*)\"\)");
    
    for f in files:
      fh, tmpfn = tempfile.mkstemp();
      ff = os.path.join(resultdir, f);
      #print ff        
        
      with open(tmpfn, 'w') as newfile:
        with open(ff) as parfile:
          for line in parfile:
            #print line
            m = rec.match(line);
            if m != None:
              pn = m.group('parname');
              if pn != 'NoInitialTransform':
                pathn, filen = os.path.split(pn);
                filen = os.path.join(resultdir, filen);
                newfile.write(line.replace(pn, filen));
              else:
                newfile.write(line);
            else:
              newfile.write(line);
                            
      os.close(fh);
      os.remove(ff);
      shutil.move(tmpfn, ff);


def setMetricParameterFile(parameterFile, metric):
    """Replaces the metric in the parameter file
       
    Arguments:
        parameterFile (str): the parameter file
        metric (str): the metric to use
        
    Notes:
        Used to replace the metric when inverse transform is estimated
    """
    #TODO: Ideally have a full parameter parser, or switch to SimpleElastix
    fh, tmpfn = tempfile.mkstemp();     
    rec = re.compile("\(Metric \"(?P<parname>.*)\"\)");
    mset = False;
    
    with open(tmpfn, 'w') as newfile:
      with open(parameterFile) as parfile:
        for line in parfile:
          #print line
          m = rec.match(line);
          if m != None:
            pn = m.group('parname');
            newfile.write(line.replace(pn, metric));
            mset = True;
          else:
            newfile.write(line);
           
    if not mset:
      newfile.write("(Metric \"" + metric + "\")\n");
             
    os.close(fh);
    os.remove(parameterFile);
    shutil.move(tmpfn, parameterFile);



def resultDataFile(resultdir):
    """Returns the mhd result file in a result directory
    
    Arguments:
        resultdir (str): Path to elastix result directory.
        
    Returns:
        str: The mhd file in the result directory.
    
    """
    
    files = os.listdir(resultdir);
    files = [x for x in files if re.match('.*.mhd', x)];
    files.sort();
    
    if files == []:
      raise RuntimeError('Cannot find a valid result data file in ' + resultdir);
    
    return os.path.join(resultdir, files[0])


def transformFileSizeAndSpacing(transformfile):
    """Parse the image size and spacing from a transformation parameter file

    Arguments:
        transformfile (str): File name of the transformix parameter file.
        
    Returns:
        (float, float): the image size and spacing
    """
    
    resi = re.compile("\(Size (?P<size>.*)\)");
    resp = re.compile("\(Spacing (?P<spacing>.*)\)");
    
    si = None;
    sp = None;
    with open(transformfile) as parfile:
        for line in parfile:
            #print line;
            m = resi.match(line)
            if m != None:
                pn = m.group('size');
                si = pn.split();
                #print si
                
            m = resp.match(line);
            if m != None:
                pn = m.group('spacing');
                sp = pn.split();
                #print sp 
    
        parfile.close();
    
    si = [float(x) for x in si];
    sp = [float(x) for x in sp];
    
    return si, sp

    
def setTransformFileSizeAndSpacing(transformfile, size, spacing):
    """Replaces size and scale in the transformation parameter file
    
    Arguments:
        transformfile (str): transformation parameter file
        size (tuple): the new image size
        spacing (tuplr): the new image spacing 
    """
    
    resi = re.compile("\(Size (?P<size>.*)\)");
    resp = re.compile("\(Spacing (?P<spacing>.*)\)");
    
    fh, tmpfn = tempfile.mkstemp();
    
    si = [int(x) for x in size];
    
    with open(transformfile) as parfile:        
        with open(tmpfn, 'w') as newfile:
            for line in parfile:
                #print line
                
                m = resi.match(line)
                if m != None:
                    newfile.write("(Size %d %d %d)" % si);
                else:
                    m = resp.match(line)
                    if m != None:
                        newfile.write("(Spacing %d %d %d)" % spacing);
                    else:
                        newfile.write(line);
            
            newfile.close();               
            parfile.close();
            
            os.remove(transformfile);
            shutil.move(tmpfn, transformfile);
        

def rescaleSizeAndSpacing(size, spacing, scale):
    """Rescales the size and spacing
    
    Arguments:
        size (tuple): image size
        spacing (tuple): image spacing
        scale (tuple): the scale factor 
    
    Returns:
        (tuple, tuple): new size and spacing
    """   

    si = [int(x * scale) for x in size];
    sp = spacing / scale;
    
    return si, sp



##############################################################################
# Elastix Runs
##############################################################################

def alignData(fixedImage, movingImage, affineParameterFile, bSplineParameterFile = None, resultDirectory = None, processes = mp.cpu_count()):
    """Align images using elastix, estimates a transformation :math:`T:` fixed image :math:`\\rightarrow` moving image.
    
    Arguments:
        fixedImage (str): image source of the fixed image (typically the reference image)
        movingImage (str): image source of the moving image (typically the image to be registered)
        affineParameterFile (str or None): elastix parameter file for the primary affine transformation
        bSplineParameterFile (str or None): elastix parameter file for the secondary non-linear transformation
        resultDirectory (str or None): elastix result directory
        processes (int): number of threads to use
        
    Returns:
        str: path to elastix result directory
    """
    
    checkElastixInitialized();
    global ElastixBinary;
    
    # result directory
    if resultDirectory == None:
        resultDirectory = tempfile.gettempdir();
    
    if not os.path.exists(resultDirectory):
        os.mkdir(resultDirectory);
    
    # run elastix
    if bSplineParameterFile is None:
      cmd = '%s -threads %d -m %s -f %s -p %s -out %s' % (ElastixBinary, processes, movingImage, fixedImage, affineParameterFile, resultDirectory);
    elif affineParameterFile is None:
      cmd = '%s -threads %d -m %s -f %s -p %s -out %s' % (ElastixBinary, processes, movingImage, fixedImage, bSplineParameterFile, resultDirectory);
    else:
      cmd = '%s -threads %d -m %s -f %s -p %s -p %s -out %s' % (ElastixBinary, processes, movingImage, fixedImage, affineParameterFile, bSplineParameterFile, resultDirectory);
    
    res = os.system(cmd);
    
    if res != 0:
        raise RuntimeError('alignData: failed executing: ' + cmd);
    
    return resultDirectory


def transformData(source, sink = [], transformParameterFile = None, transformDirectory = None, resultDirectory = None):
  """Transform a raw data set to reference using the elastix alignment results
  
  If the map determined by elastix is 
  :math:`T \\mathrm{fixed} \\rightarrow \\mathrm{moving}`, 
  transformix on data works as :math:`T^{-1}(\\mathrm{data})`.
      
  Arguments:
      source (str or array): image source to be transformed
      sink (str, [] or None): image sink to save transformed image to. if [] return the default name of the data file generated by transformix.
      transformParameterFile (str or None): parameter file for the primary transformation, if None, the file is determined from the transformDirectory.
      transformDirectory (str or None): result directory of elastix alignment, if None the transformParameterFile has to be given.
      resultDirectory (str or None): the directorty for the transformix results
      
  Returns:
      array or str: array or file name of the transformed data
  """
  
  checkElastixInitialized();
  global TransformixBinary;    
  
  # image
  delete_image = None;
  if isinstance(source, np.ndarray):
    imgname = os.path.join(tempfile.gettempdir(), 'elastix_input.tif');
    io.writeData(source, imgname);
    delete_image = imgname;
  elif isinstance(source, str):
    if io.dataFileNameToType(source) == "TIF":
      imgname = source;
    else:
      imgname = os.path.join(tempfile.gettempdir(), 'elastix_input.tif');
      io.transformData(source, imgname);
      delete_image =  imgname;
  else:
    raise RuntimeError('transformData: source not a string or array');

  # result directory
  delete_result_directory = None;
  if resultDirectory == None:
    resultdirname = os.path.join(tempfile.tempdir, 'elastix_output');
    delete_result_directory = resultdirname;
  else:
    resultdirname = resultDirectory;
     
  if not os.path.exists(resultdirname):
    os.makedirs(resultdirname);
  
  # tranformation parameter
  transformparameterdir, transformParameterFile = transformDirectoryAndFile(transformParameterFile = transformParameterFile, transformDirectory = transformDirectory);
  
  setPathTransformFiles(transformparameterdir);
 
  #transformix -in inputImage.ext -out outputDirectory -tp TransformParameters.txx
  cmd = '%s -in %s -out %s -tp %s' % (TransformixBinary, imgname, resultdirname, transformParameterFile);
  
  res = os.system(cmd);
  
  if res != 0:
      raise RuntimeError('transformData: failed executing: ' + cmd);
  
  # read data and clean up
  if delete_image is not None:
      os.remove(delete_image);
  
  if sink == []:
    return resultDataFile(resultdirname);
  elif sink is None:
    resultfile = resultDataFile(resultdirname);
    result = io.readData(resultfile);
  elif isinstance(sink, str):
    resultfile = resultDataFile(resultdirname);
    result = io.convertData(resultfile, sink);
  else:
    raise RuntimeError('transformData: sink not valid!');
    
  if delete_result_directory is not None:
    shutil.rmtree(delete_result_directory);
  
  return result;


def deformationField(sink = [], transformParameterFile = None, transformDirectory = None, resultDirectory = None):
    """Create the deformation field T(x) - x
    
    The map determined by elastix is 
    :math:`T \\mathrm{fixed} \\rightarrow \\mathrm{moving}`
        
    Arguments:
        sink (str, [] or None): image sink to save the transformation field; if [] return the default name of the data file generated by transformix.
        transformParameterFile (str or None): parameter file for the primary transformation, if None, the file is determined from the transformDirectory.
        transformDirectory (str or None): result directory of elastix alignment, if None the transformParameterFile has to be given.
        resultDirectory (str or None): the directorty for the transformix results
        
    Returns:
        array or str: array or file name of the transformed data
    """
    checkElastixInitialized();
    global TransformixBinary;    
    
    # result directory
    delete_result_directory = None;
    if resultDirectory == None:
        resultdirname = os.path.join(tempfile.tempdir, 'elastix_output');
        delete_result_directory = resultdirname;
    else:
        resultdirname = resultDirectory;
        
    if not os.path.exists(resultdirname):
        os.makedirs(resultdirname);
       
    # setup transformation 
    transformparameterdir, transformParameterFile = transformDirectoryAndFile(transformParameterFile = transformParameterFile, transformDirectory = transformDirectory); 
    setPathTransformFiles(transformparameterdir);
   
    #transformix -in inputImage.ext -out outputDirectory -tp TransformParameters.txt
    cmd = '%s -def all -out %s -tp  %s' % (TransformixBinary, resultdirname, transformParameterFile)
    
    res = os.system(cmd);
    
    if res != 0:
        raise RuntimeError('deformationField: failed executing: ' + cmd);
    
    
    # read result and clean up
    if sink == []:
        return resultDataFile(resultdirname);
    elif sink is None:
        resultfile = resultDataFile(resultdirname);
        result = io.readData(resultfile);
    elif isinstance(sink, str):
        resultfile = resultDataFile(resultdirname);
        result = io.convertData(resultfile, sink);
    else:
        raise RuntimeError('deformationField: sink not valid!');
        
    if delete_result_directory is not None:
      shutil.rmtree(delete_result_directory);
    
    return result;


def deformationDistance(deformationField, sink = None, scale = None):
    """Compute the distance field from a deformation vector field
    
    Arguments:
        deformationField (str or array): source of the deformation field determined by :func:`deformationField`
        sink (str or None): image sink to save the deformation field to
        scale (tuple or None): scale factor for each dimension, if None = (1,1,1)
        
    Returns:
        array or str: array or file name of the transformed data
    """
    
    deformationField = io.readData(deformationField);
    
    df = np.square(deformationField);
    if not scale is None:
        for i in range(3):
            df[:,:,:,i] = df[:,:,:,i] * (scale[i] * scale[i]);
    df = np.sqrt(np.sum(df, axis = 3));
    
    return io.writeData(sink, df);
    

def writePoints(filename, points, indices = True, binary = True):
    """Write points as elastix/transformix point file
    
    Arguments:
        filename (str): file name of the elastix point file.
        points (array or str): source of the points.
        indices (bool): write as pixel indices or physical coordiantes
    
    Returns:
        str : file name of the elastix point file
    """
    
    points = io.readPoints(points);
    #points = points[:,[1,0,2]]; # points in ClearMap (y,x,z) -> permute to (x,y,z)
    
    if binary:
      with open(filename, 'wb') as pointfile:
        if indices:
          np.array(1, dtype = np.int64).tofile(pointfile)
        else:
          np.array(0, dtype = np.int64).tofile(pointfile)
          
        num_points = np.array(len(points), dtype = np.int64);
        num_points.tofile(pointfile);

        points = np.asarray(points, dtype = np.double);
        points.tofile(pointfile);

        pointfile.close();        
        
    else:
      with open(filename, 'w') as pointfile:
        if indices:
          pointfile.write('index\n')
        else:
          pointfile.write('point\n')
      
        pointfile.write(str(points.shape[0]) + '\n');
        np.savetxt(pointfile, points, delimiter = ' ', newline = '\n', fmt = '%.5e')
        pointfile.close();
    
    return filename;


def readPoints(filename, indices = True, binary = True):
    """Parses the output points from the output file of transformix
    
    Arguments:
        filename (str): file name of the transformix output file
        indices (bool): if True return pixel indices otherwise float coordinates
        
    Returns:
        points (array): the transformed coordinates     
    """
    
    if binary:
      with open(filename) as f:
        index = np.fromfile(f, dtype=np.int64, count = 1)[0];
        print (index)
        if index == 0:
          indices = False;
        else:
          indices = True;
        
        num_points = np.fromfile(f, dtype=np.int64, count = 1)[0];
        print (num_points)
        if num_points == 0:
          return np.zeros((0,3));
        
        points = np.fromfile(f, dtype = np.double);
        #print points.shape
        points = np.reshape(points, (num_points,3));
        
        f.close();
        
      return points;
    
    else: # text file
    
      with open(filename) as f:
        lines = f.readlines()
        f.close();
      
      num_points = len(lines);
      
      if num_points == 0:
        return np.zeros((0,3));
      
      points = np.zeros((num_points, 3));
      k = 0;
      for line in lines:
        ls = line.split();
        if indices:
          for i in range(0,3):
            points[k,i] = float(ls[i+22]);
        else:
          for i in range(0,3):
            points[k,i] = float(ls[i+30]);
        
        k += 1;
      
      return points;



def transformPoints(source, sink = None, transformParameterFile = None, transformDirectory = None, indices = True,
                    resultDirectory = None, tmpFile = None, binary = True):
    """Transform coordinates math:`x` via elastix estimated transformation to :math:`T(x)`

    Note the transformation is from the fixed image coorindates to the moving image coordiantes.
    
    Arguments:
        source (str): source of the points
        sink (str or None): sink for transformed points
        transformParameterFile (str or None): parameter file for the primary transformation, if None, the file is determined from the transformDirectory.
        transformDirectory (str or None): result directory of elastix alignment, if None the transformParameterFile has to be given.
        indices (bool): if True use points as pixel coordinates otherwise spatial coordinates.
        resultDirectory (str or None): elastic result directory
        tmpFile (str or None): file name for the elastix point file.
        
    Returns:
        array or str: array or file name of transformed points
    """
        
    global TransformixBinary;    
    
    checkElastixInitialized();    
    global ElastixSettings;

    # input point file
    if tmpFile == None:
      if binary:
        tempfile.gettempdir();
        tmpFile = os.path.join(tempfile.tempdir, 'elastix_input.bin');
      else:
        tempfile.gettempdir();
        tmpFile = os.path.join(tempfile.tempdir, 'elastix_input.txt');
    
    delete_point_file = None;
    if isinstance(source, str):
      if len(source) > 3 and source[-3:] in ['txt', 'bin']:
        if source[-3:] == 'txt':
          binary = False; 
        if source[-3] == 'bin':
          binary = True;
        pointfile = source;
      else:
        points = io.readPoints(source);
        pointfile = tmpFile;
        delete_point_file = tmpFile;
        writePoints(pointfile, points, indices = indices, binary = binary);
    elif isinstance(source, np.ndarray):
        pointfile = tmpFile;
        delete_point_file = tmpFile;
        writePoints(pointfile, source, indices = indices, binary = binary);
    else:
        raise RuntimeError('transformPoints: source not string or array!');
    
    # result directory
    delete_result_directory = None;
    if resultDirectory == None:
      outdirname = os.path.join(tempfile.tempdir, 'elastix_output');
      delete_result_directory = outdirname;
    else:
      outdirname = resultDirectory;
        
    if not os.path.exists(outdirname):
      os.makedirs(outdirname);
    
    #transform
    transformparameterdir, transformparameterfile = transformDirectoryAndFile(transformParameterFile = transformParameterFile,
                                                                              transformDirectory = transformDirectory);
    setPathTransformFiles(transformparameterdir);
    
    #run transformix   
    cmd = '%s -def %s -out %s -tp %s' % (TransformixBinary, pointfile, outdirname, transformparameterfile);
    
    res = os.system(cmd);
    
    if res != 0:
        raise RuntimeError('failed executing ' + cmd);
    
    # read data and clean up
    if delete_point_file is not None:
      os.remove(delete_point_file);
    
    #read data / file 
    if sink == []: # return sink as file name
      if binary:
        return os.path.join(outdirname, 'outputpoints.bin')
      else:
        return os.path.join(outdirname, 'outputpoints.txt')
    
    else:
        if binary:
          transpoints = readPoints(os.path.join(outdirname, 'outputpoints.bin'), indices = indices, binary = True);
        else:
          transpoints = readPoints(os.path.join(outdirname, 'outputpoints.txt'), indices = indices, binary = False); 
        
        if delete_result_directory is not None:
          shutil.rmtree(delete_result_directory);
        
        return io.writePoints(sink, transpoints);

        
        
def inverseTransform(fixedImage, affineParameterFile, bSplineParameterFile = None, transformParameterFile = None, transformDirectory = None, resultDirectory = None, processes = mp.cpu_count()):
    """Estimate inverse tranformation :math:`T^{-1}:` moving image :math:`\\rightarrow` fixed image.
    
    Arguments:
        fixedImage (str): image source of the fixed image (typically the reference image)
        affineParameterFile (str): the paramter file for the original affine transformation
        bSplineParameterFile (str): the paramter file for the original b-spline transformation
        transformDirectory (str): elastic result directory of the original transform
        resultDirectory (str or None): elastic result directory of the inverse transform
        
    Returns:
        str: path to elastix result directory
    """
    
    checkElastixInitialized();
    global ElastixBinary;
    
    # result directory
    if resultDirectory == None:
        resultDirectory = tempfile.gettempdir();
    
    if not os.path.exists(resultDirectory):
        os.mkdir(resultDirectory);
    
    # transformation files
    transformparameterdir, transformparameterfile = transformDirectoryAndFile(transformParameterFile = transformParameterFile, transformDirectory = transformDirectory);    
    setPathTransformFiles(transformparameterdir);
    
    #set metric of the parameter files
    if bSplineParameterFile is not None:
      _, bsplinefile = os.path.split(bSplineParameterFile);
      bsplinefile    = os.path.join(resultDirectory, bsplinefile);
      shutil.copyfile(bSplineParameterFile, bsplinefile);
      setMetricParameterFile(bsplinefile, 'DisplacementMagnitudePenalty');
    else:
      bsplinefile = None;
      
    if affineParameterFile is not None:
      _, affinefile = os.path.split(affineParameterFile);
      affinefile    = os.path.join(resultDirectory, affinefile);
      shutil.copyfile(affineParameterFile, affinefile);
      setMetricParameterFile(affinefile, 'DisplacementMagnitudePenalty');
    else:
      affinefile = None;
    
    # run elastix
    if bsplinefile is None:
      cmd = '%s -threads %d -m %s -f %s -t0 %s -p %s -out %s' % (ElastixBinary, processes, fixedImage, fixedImage, transformparameterfile, affinefile,  resultDirectory);
    elif affinefile is None:
      cmd = '%s -threads %d -m %s -f %s -t0 %s -p %s -out %s' % (ElastixBinary, processes, fixedImage, fixedImage, transformparameterfile, bsplinefile, resultDirectory);
    else:
      cmd = '%s -threads %d -m %s -f %s -t0 %s -p %s -p %s -out %s' % (ElastixBinary, processes, fixedImage, fixedImage, transformparameterfile, affinefile, bsplinefile, resultDirectory);    
    
    res = os.system(cmd);
    
    if res != 0:
        raise RuntimeError('inverseTransform: failed executing: ' + cmd);
    
    return resultDirectory
  




  

    
    
    
