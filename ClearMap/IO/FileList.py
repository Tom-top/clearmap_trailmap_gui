# -*- coding: utf-8 -*-
"""
Interface to read/write image stacks saved as a list of files

The filename is given as regular expression as described 
`here <https://docs.python.org/2/library/re.html>`_.

It is assumd that there is a single digit like regular expression in the file
name, i.e. ``\\d{4}`` indicates a placeholder for an integer with four digits using traling 0s
and ``\\d{}`` would jus asume an integer with variable size.

For example: ``/test\d{4}.tif`` or  ``/test\d{}.tif``
    
Examples:
    >>> import os, numpy
    >>> import ClearMap.Settings as settings
    >>> import ClearMap.IO.FileList as fl
    >>> filename = os.path.join(settings.ClearMapPath, 'Test/Data/FileList/test\d{4}.tif')  
    >>> data = numpy.random.rand(20,50,10);
    >>> data = data.astype('int32');
    >>> fl.writeData(filename, data);
    >>> img = fl.readData(filename);  
    >>> print img.shape
    (20, 50, 10)
"""
#:copyright: Copyright 2015 by Christoph Kirst, The Rockefeller University, New York City
#:license: GNU, see LICENSE.txt for details.

import numpy as np
import os
import re
import natsort
import multiprocessing 


import ClearMap.IO as io

from ClearMap.Utils.ProcessWriter import ProcessWriter;
 

def firstFile(source, sort = True):
  """Returns the first file that matches the soruce specification
  
  Arguments:
      source (str): file name as regular expression
      sort (bool): naturally sort the files if True, if False and arbirary file that matched first is returned
      
  Returns:
      string or None: first file name that matches the regular expression or None if no match found
    
  Note:
    For large file numbers speed can be imporved by setting the option sort = False
  """
  
  #get path        
  (fpath, fname) = os.path.split(source)
  search = re.compile(fname).search 
  
  if sort:
    for fname in natsort.natsorted(os.listdir(fpath)):
      if search(fname):
        return os.path.join(fpath, fname);
    return None;
  else:
    for fname in os.listdir(fpath):
      if search(fname):
        return os.path.join(fpath, fname);
    return None;


def readFileList(filename):
    """Returns list of files that match the regular expression
    
    Arguments:
        filename (str): file name as regular expression
    
    Returns:
        str, list: path of files, file names that match the regular expression
    """
    
    #get path        
    (fpath, fname) = os.path.split(filename)
    fnames = os.listdir(fpath);
    #fnames = [os.path.join(fpath, x) for x in fnames];
    
    searchRegex = re.compile(fname).search    
    fl = [ l for l in fnames for m in (searchRegex(l),) if m]  
    
    if fl == []:
        raise RuntimeError('no files found in ' + fpath + ' that match ' + fname + ' !');
    
    fl.sort();
        
    return fpath, fl;
    

def splitFileExpression(filename):
    """Split the regular expression at the digit place holder
    
    Arguments:
        filename (str): file name as regular expression
    
    Returns:
        tuple: file header, file extension, digit format
    """
    
    searchRegex = re.compile('.*\\\\d\{(?P<digit>\d)\}.*').search
    m = searchRegex(filename); 
    
    if not m is None:
        digits = int(m.group('digit'));
        searchRegex = re.compile('.*(?P<replace>\\\\d\{\d\}).*').search
        m = searchRegex(filename);       
        fs = filename.split(m.group('replace'));
        if not len(fs) == 2:
            raise RuntimeError('FileList: more than a single replacement for z indexing!');
        fileheader = fs[0];
        fileext    = fs[1];
        
    else: #fileheader is given
        digits = 4;
        fileheader = filename;
        fileext    = '.tif';
        
    digitfrmt = "%." + str(digits) + "d";        
        
    return (fileheader, fileext, digitfrmt);


def fileExpressionToFileName(filename, z):
    """Insert a number into the regular expression
    
    Arguments:
        filename (str): file name as regular expression
        z (int): z slice index
    
    Returns:
        str: file name
    """    
    
    (fileheader, fileext, digitfrmt) = splitFileExpression(filename);
    return fileheader + (digitfrmt % z) + fileext;    


def dataSize(filename, **args):
    """Returns size of data stored as a file list
    
    Arguments:
        filename (str): file name as regular expression
        x,y,z (tuple): data range specifications
    
    Returns:
        tuple: data size
    """
    
    fp, fl = readFileList(filename);
    nz = len(fl);
    
    d2 = io.dataSize(os.path.join(fp, fl[0]));
    if not len(d2) == 2:
        raise RuntimeError("FileList: importing multiple files of dim %d not supported!" % len(d2));
    #d2 = (d2[1], d2[0]);
    
    dims = d2 + (nz,);
    return io.dataSizeFromDataRange(dims, **args);
   

def dataType(filename, **args):
    """Returns type of data stored as a file list
    
    Arguments:
        filename (str): file name as regular expression
        x,y,z (tuple): data range specifications
    
    Returns:
        tuple: data size
    """
    
    fp, fl = readFileList(filename);
    nz = len(fl);
    
    d2 = io.dataSize(os.path.join(fp, fl[0]));
    if not len(d2) == 2:
        raise RuntimeError("FileList: importing multiple files of dim %d not supported!" % len(d2));
    #d2 = (d2[1], d2[0]);
    
    dims = d2 + (nz,);
    return io.dataSizeFromDataRange(dims, **args);
   
   
   
def dataZSize(filename, z = all, **args):
    """Returns size of data stored as a file list
    
    Arguments:
        filename (str): file name as regular expression
        z (tuple): z data range specification
    
    Returns:
        int: z data size
    """
    
    fp, fl = readFileList(filename);
    nz = len(fl);
    return io.toDataSize(nz, r = z);



def readDataFiles(filename, x = all, y = all, z = all, **args):
    """Read data from individual images assuming they are the z slices

    Arguments:
        filename (str): file name as regular expression
        x,y,z (tuple): data range specifications
    
    Returns:
        array: image data
    """   
    
    fpath, fl = readFileList(filename);
    nz = len(fl);
    
    #read first image to get data size and type
    rz = io.toDataRange(nz, r = z);
    sz = io.toDataSize(nz, r = z);
    fn = os.path.join(fpath, fl[rz[0]]); 
    img = io.readData(fn, x = x, y = y);
    nxy = img.shape;
    data = np.zeros( (sz,) + nxy , dtype = img.dtype);
    data[0,:,:] = img;
    
    for i in range(rz[0]+1, rz[1]):
        fn = os.path.join(fpath, fl[i]);    
        data[i-rz[0],:,:] = io.readData(fn, x = x, y = y);
    data = data.transpose([1,2,0]);
    
    return data;

     
def readData(filename, **args):
    """Read image stack from single or multiple images
    
    Arguments:
        filename (str): file name as regular expression
        x,y,z (tuple): data range specifications
    
    Returns:
        array: image data
    """
    
    if os.path.exists(filename):
         return io.readData(filename, **args);
    else:
         return readDataFiles(filename, **args);



def writeData(filename, data, startIndex = 0):
    """Write image stack to single or multiple image files
    
    Arguments:
        filename (str): file name as regular expression
        data (array): image data
        startIndex (int): index of first z-slice
    
    Returns:
        str: file name as regular expression
    """
    
    #create directory if not exsits
    io.createDirectory(filename); 

    #check for the \d{xx} part of the regular expression -> if not assume file header    
    (fileheader, fileext, digitfrmt) = splitFileExpression(filename);

    d = len(data.shape);
    if d == 2:
        fname = fileheader + (digitfrmt % startIndex) + fileext;
        io.writeData(fname, data);
        return fname;
    else:
        nz = data.shape[2];
        for i in range(nz):
            fname = fileheader + (digitfrmt % (i + startIndex)) + fileext;
            io.writeData(fname, data[:,:,i]);
        return filename;


def copyData(source, sink):
    """Copy a data file from source to sink when for entire list of files
    
    Arguments:
        source (str): file name pattern of source
        sink (str): file name pattern of sink
    
    Returns:
        str: file name patttern of the copy
    """ 
    
    (fileheader, fileext, digitfrmt) = splitFileExpression(sink);
    
    fp, fl = readFileList(source);
    
    for i in range(len(fl)):
        io.copyFile(os.path.join(fp, fl[i]), fileheader + (digitfrmt % i) + fileext);
    
    return sink



def _convertToNPYParallel(arg):
    """Convert helper function to use for parallel conversion of image slices into a numpy memmap"""
    
    source = arg[0];
    sink = arg[1];
    z = arg[2];
    ii = arg[3];
    nn = arg[4];
    
    if ii is not None:
      pw = ProcessWriter(ii);
      pw.write("convertingData: converting image %d / %d" % (ii, nn))    
      #pw.write('%s -> %s' % (fileSource, fileSink));

    data = io.readData(source);
    sink = np.lib.format.open_memmap(sink, 'r+');
    sink[:,:,z] = data[:,:];
    

def convertToNPY(source, sink, verbose = True, processes = all):
  """Converts list of image files to a numpy array on disk
  
  Arguments:
    source (str): file expression of source
    sink (str): filename
  
  Return:
    str: filename of numpy file
  """
  
  # read files
  fpath, flist = readFileList(source);
  nz = len(flist);
  if nz == 0:
    raise RuntimeError('No files matching %s found  to convert !' % source);
  
  firstfile = os.path.join(fpath, flist[0]);
  shape = io.dataSize(firstfile) + (nz,);
  
  data = io.readData(firstfile);
  dtype = data.dtype;
  #TODO: dtype = io.dataType(firstfile);
  
  # create memmap
  memmap = np.lib.format.open_memmap(sink, 'w+', shape = shape, dtype = dtype, fortran_order = True);
    
  if processes is None:
    processes = 1;
  if processes is all:
    processes = multiprocessing.cpu_count();
    
  if processes > 1: #parallel processing
    pool = multiprocessing.Pool(processes=processes);
    argdata = [];
   
    for i,z in enumerate(range(nz)):
      if verbose:
        argdata.append( (os.path.join(fpath, flist[z]), sink, z, (i+1), (nz+1)) );    
      else:
        argdata.append( (os.path.join(fpath, flist[z]), sink, z, None, None) );
      
    pool.map(_convertToNPYParallel, argdata);
    
  else: # sequential processing
    for i,z in enumerate(range(nz)):
      if verbose:
        print ("convertData: converting image %d / %d" % (i+1, nz+1));
      
      fileSource = os.path.join(fpath, flist[z]);
      data = io.readData(fileSource);
      memmap[:,:,z] = data[:,:];
        
    memmap.flush();
    
  return sink;


def _cropParallel(arg):
    """Crop helper function to use for parallel cropping of image slices"""
    
    fileSource = arg[0];
    fileSink = arg[1];
    x = arg[2];
    y = arg[3];
    ii = arg[4];
    nn = arg[5];
    
    if ii is not None:
        pw = ProcessWriter(ii);
        pw.write("cropData: corpping image %d / %d" % (ii, nn))    
        #pw.write('%s -> %s' % (fileSource, fileSink));
    
    data = io.readData(fileSource, x = x, y = y);
    io.writeData(fileSink, data);

  
def cropData(source, sink = None, x = all, y = all, z = all, adjustOverlap = False, verbose = True, processes = all):
  """Crop source from start to stop point
  
  Arguments:
    source (str or array): filename or data array of source
    sink (str or None): filename or sink
    x,y,z (tuple or all): the range to crop the data to
    adjustOverlap (bool): correct overlap meta data if exists
  
  Return:
    str or array: array or filename with cropped data
  """
  
  if sink is None:
    return readDataFiles(source, x = x, y = y, z = z);
  else: # sink assumed to be file expression
  
    if not io.isFileExpression(sink):
      raise RuntimeError('cropping data to different format not supported!')
      
    fileheader, fileext, digitfrmt = splitFileExpression(sink);

    #read first image to get data size and type
    fp, fl = readFileList(source);    
    nz = len(fl);
    rz = io.toDataRange(nz, r = z);
  
    if adjustOverlap: #change overlap in first file 
      try: 
        fn = os.path.join(fp, fl[0]);
        info = io.readMetaData(fn, info = ['description', 'overlap', 'resolution']);
        description = str(info['description']);
        overlap = numpy.array(info['overlap'], dtype = float);
        resolution = numpy.array(info['resolution'], dtype = float);
        
      except:
        raise RuntimeWarning('could not modify overlap!')
      
      fullsize = io.dataSize(fn);
      data = io.readData(fn, x = x, y = y);
      
      #overlap in pixels
      poverlap = overlap[:2] / resolution[:2];
      print (poverlap)
      
      #cropped pixel
      xr = io.toDataRange(fullsize[0], r = x);
      yr = io.toDataRange(fullsize[1], r = y);

      print (xr)
      print (yr)
      print (fullsize)

      poverlap[0] = poverlap[0] - xr[0] - (fullsize[0] - xr[1]);
      poverlap[1] = poverlap[1] - yr[0] - (fullsize[1] - yr[1]);
      print (poverlap)
      
      #new overlap in microns
      overlap = poverlap * resolution[:2];
      
      #check for consistency      
      if numpy.abs(fullsize[0]-xr[1] - xr[0]) > 1 or numpy.abs(fullsize[1]-yr[1] - yr[0]) > 1:
        raise RuntimeWarning('cropping is inconsistent with overlap )modification!');

      #change image description
      import ClearMap.IO.TIF as CMTIF
      description = CMTIF.changeOMEMetaDataString(description, {'overlap': overlap});
      print (len(description))
      
      
      #write first file
      fnout = fileheader + (digitfrmt % 0) + fileext;
      io.writeData(fnout, data, info  = description);
      
      zr = range(rz[0]+1, rz[1]);
    else:
      zr = range(rz[0], rz[1]);
    
    print (zr)
    nZ = len(zr); 
    
    if processes is None:
      processes = 1;
    if processes is all:
      processes = multiprocessing.cpu_count();
    
    if processes > 1: #parallel processing
      pool = multiprocessing.Pool(processes=processes);
      argdata = [];
   
      for i,z in enumerate(zr):
        if verbose:
          argdata.append( (os.path.join(fp, fl[z]), fileheader + (digitfrmt % (i+1)) + fileext, x, y, (i+1), (nZ+1)) );    
        else:
          argdata.append( (os.path.join(fp, fl[z]), fileheader + (digitfrmt % (i+1)) + fileext, x, y, None, None) );
      
      pool.map(_cropParallel, argdata);
    
    else: # sequential processing
      for i,z in enumerate(zr):
        if verbose:
          print ("cropData: corpping image %d / %d" % (i+1, nZ+1));
        
        fileSource = os.path.join(fp, fl[z]);
        data = io.readData(fileSource, x = x, y = y);
        
        fileSink = fileheader + (digitfrmt % (i+1)) + fileext
        io.writeData(fileSink, data);
    
    return sink;



def readMetaData(source, info = all, sort = True):
  """Reads the meta data from the image files
  
  Arguments:
    source: the data source
    info (list or all): optional list of keywords
    sort (bool): if True use first file to infer meta data, otherwise arbitrary file
  
  Returns:
    object: an object with the meta data
  """
    
  firstfile = firstFile(source, sort = sort);
  
  mdata = io.readMetaData(firstfile, info = info);
  
  if 'size' in mdata.keys():
    mdata['size'] = dataSize(source);

  return mdata;


def test():    
    """Test FileList module"""  
    import ClearMap.IO.FileList as self
    reload(self)
    
    from iDISCO.Parameter import iDISCOPath
    import os
    import numpy

    basedir = iDISCOPath();
    fn = os.path.join(basedir,'Test/Data/FileList/test\d{4}.tif')  

    data = numpy.random.rand(20,50,10);
    data[5:15, 20:45, 2:9] = 0;
    data = 20 * data;
    data = data.astype('int32');
    
    print ("writing raw image to: " + fn);    
    self.writeData(fn, data);

    print ("Loading raw image from: " + fn);
    img = self.readData(fn);  
    print ("Image size: " + str(img.shape))
    
    diff = img - data;
    print (diff.max(), diff.min())

    
    fn = os.path.join(basedir,'Test/Data/OME/16-17-27_0_8X-s3-20HF_UltraII_C00_xyz-Table Z\d{4}.ome.tif')        
    
    fp, fl = self.readFileList(fn);
    print ("Found " + str(len(fl)) + " images!")    
    
    #dataSize
    print ("dataSize  is %s" % str(self.dataSize(fn)))
    print ("dataZSize is %s" % str(self.dataZSize(fn)))
    
    print ("dataSize  is %s" % str(self.dataSize(fn, x = (10,20))))
    print ("dataZSize is %s" % str(self.dataZSize(fn)))
    
    
    img = self.readData(fn, z = (17,all));  
    print ("Image size: " + str(img.shape)   ) 
    
    import iDISCO.Visualization.Plot as plt
    plt.plotTiling(img)
        

if __name__ == "__main__":
    
    test();
    
