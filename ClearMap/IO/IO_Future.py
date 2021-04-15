# -*- coding: utf-8 -*-
"""
IO interface to read image and array data

This is the main module to distribute the reading and writing of individual data formats to the specialized sub-modules.

A data source is either a specification to image or other array data.
  
See :mod:`ClearMap.IO` for details.
"""
__author__    = 'Christoph Kirst <ckirst@rockefeller.edu>'
__license__   = 'MIT License <http://www.opensource.org/licenses/mit-license.php>'


# notes: it might be a good idea to have wrapper classes for each data type which deal with the underlying data format etc.
#        and then have an io module that generates the data class from the input
#        e.g. ArrayNYP, ImageFileList, etc... base class would be Data -> Image , Array etc, then Image -> ImageFileList, ImageRaw etc..
#        classes would have info printing routines, metadata etc..

import re
import numpy
import importlib

from .FileUtils import *
from .Region import *

import ClearMap.Utils.ParameterTools as par

##############################################################################
# File associations 
##############################################################################

arrayFileExtensions = ["csv", "txt", "npy", "vtk", "ims"];
"""list of extensions supported as an array file"""

arrayFileTypes = ["CSV", "NPY", "VTK", "Imaris"];
"""list of array file types"""

arrayFileExtensionToType = {"csv" : "CSV", "txt" : "CSV", "npy" : "NPY", "vtk" : "VTK", "ims" : "Imaris"};
"""map from point file extensions to point file types"""


imageFileExtensions = ["tif", "tiff", "mhd", "raw", "ims", "nrrd"];
"""list of extensions supported as a image data file"""

imageFileTypes = ["FileList", "TIF", "RAW", "NRRD", "Imaris"]
"""list of image data file types"""

imageFileExtensionToType = { "tif" : "TIF", "tiff" : "TIF", "raw" : "RAW", "mhd" : "RAW", "nrrd": "NRRD", "ims" : "Imaris"}
"""map from image file extensions to image file types"""


dataFileExtensions = arrayFileExtensions + imageFileExtensions
"""list of extensions supported data file"""

dataFileTypes = arrayFileTypes + imageFileTypes
"""list of data file types"""

dataFileExtensionToType = par.joinParameter(arrayFileExtensionToType, imageFileExtensionToType);
"""map from data file extensions to data file types"""

##############################################################################
# Basic Querries 
##############################################################################

class iostr(str):
    """Wrapper to identify a string pointing to a data file
    
    Notes:
      used in parameter dictionaries to indicate an io parameter
    """
    pass


def isFileExpression(source):
    """Checks if filename is a regular expression denoting a file list
    
    Arguments:
        source (str): source file name
        
    Returns:
        bool: True if source is regular expression with a digit placeholder
    """    

    if not isinstance(source, str):
      return False;    
      
    ext = fileExtension(source);
    if not ext in dataFileExtensions:
      return False
    
    #sepcified number of digits
    searchRegex = re.compile('.*\\\\d\{(?P<digit>\d)\}.*').search
    m = searchRegex(source); 
    if not m is None:
      return True;
    
    #digits without trailing zeros \d* or 
    searchRegex = re.compile('.*\\\\d\*.*').search
    m = searchRegex(source); 
    if not m is None:
      return True;
    
    #digits without trailing zeros \d{} or 
    searchRegex = re.compile('.*\\\\d\{\}.*').search
    m = searchRegex(source); 
    if not m is None:
      return True;
      
    return False;

     
def isDataFile(source, exists = False):
    """Checks if a file has a valid data file extension useable in *ClearMap*
     
    Arguments:
        source (str): source file name
        exists (bool): if true also checks if source exists 
        
    Returns:
        bool: true if source is an data file usable in *ClearMap*
    """   
    
    if not isinstance(source, str):
        return False;    
    
    fext = fileExtension(source);
    if fext in dataFileExtensions:
      if not exists:
        return True;
      else:
        return exitsFile(source) or existsFileExpression(source);
    else:
        return False;
        

def isImageFile(source, exists = False):
    """Checks if a file has a valid image file extension useable in *ClearMap*
     
    Arguments:
        source (str): source file name
        exists (bool): if true also checks if source exists 
        
    Returns:
        bool: true if source is an image file usable in *ClearMap*
    """   
    
    if not isinstance(source, str):
        return False;    
    
    fext = fileExtension(source);
    if fext in imageFileExtensions:
      if not exists:
        return True;
      else:
        return exitsFile(source) or existsFileExpression(source);
    else:
        return False;


def isArrayFile(source, exists =False):
    """Checks if a file is a valid array data file
     
    Arguments:
        source (str): source file name
        exists (bool): if true also checks if source exists 
        
    Returns:
        bool: true if source is a array data file
    """     
    
    if not isinstance(source, str):
        return False;
    
    fext = fileExtension(source);
    if fext in pointFileExtensions:
      if not exists:
        return True;
      else:
        return exitsFile(source) or existsFileExpression(source);
    else:
        return False;


def isDataSource(source, exists = False):
  """Checks if source is a valid data source for use in *ClearMap*
   
  Arguments:
      source (str): source file name or array
      exists (bool): if true also checks if source exists 
      
  Returns:
      bool: true if source is an data source usable in *ClearMap*
  """  
  
  return (not exists and source is None) or isinstance(source, numpy.ndarray) or isDataFile(source, exists = exists);


def isSource(source, exists = False):
  """Checks if source is a valid source for use in *ClearMap*
   
  Arguments:
      source (str): source file name or array
      exists (bool): if true also checks if source exists 
      
  Returns:
      bool: true if source is an data source usable in *ClearMap*
  """  
  
  return isDataSource(source = source, exists = exists)


##############################################################################
# Type to IO Module conversions
##############################################################################

def arrayFileNameToType(filename):
    """Returns type of an array file
    
    Arguments:
        filename (str): file name
        
    Returns:
        str: array data type in :const:`arrayFileTypes`
    """       
    
    #if not self.isFile(filename):
    #    raise RuntimeError("Cannot find point file %s" % filename);
    #else:
    fext = fileExtension(filename);
    if fext in arrayFileExtensions:
        return arrayFileExtensionToType[fext];
    else:
       raise RuntimeError("Cannot determine type of array file %s with extension %s" % (filename, fext));     


def imageFileNameToType(filename):
    """Returns type of a image file
    
    Arguments:
        filename (str): file name
        
    Returns:
        str: image data type in :const:`imageFileTypes`
    """       
    
    #if not self.isFile(filename):
    #    raise RuntimeError("Cannot find point file %s" % filename);
    #else:
    fext = fileExtension(filename);
    if fext in imageFileExtensions:
        return imageFileExtensionToType[fext];
    else:
       raise RuntimeError("Cannot determine type of image file %s with extension %s" % (filename, fext));     


def dataFileNameToType(filename):
    """Returns type of a data file
    
    Arguments:
        filename (str): file name
        
    Returns:
        str: image or array data type in :const:`dataFileTypes`
        
    Notes:
        image are tried to be used first
    """      
    
    try:
      t = imageFileNameToType(filename);
    except:
      try:
        t = arrayFileNameToType(filename);
      except:
        fext = fileExtension(filename);
        raise RuntimeError("Cannot determine type of data file %s with extension %s" % (filename, fext));
  
    return t;


def arrayFileNameToModule(filename):
    """Return the module that handles io for a array file
        
    Arguments:
        filename (str): file name
        
    Returns:
        object: sub-module that handles a specific point file type
    """ 

    ft = arrayFileNameToType(filename);
    return importlib.import_module("ClearMap.IO." + ft);


def imageFileNameToModule(filename):
    """Return the module that handles io for a image file
        
    Arguments:
        filename (str): file name
        
    Returns:
        object: sub-module that handles a specific point file type
    """ 

    ft = imageFileNameToType(filename);
    return importlib.import_module("ClearMap.IO." + ft);
    
    
def dataFileNameToModule(filename):
    """Returns the module that handles io for a data file
        
    Arguments:
        filename (str): file name
        
    Returns:
        object: sub-module that handles a specific data type
    """          
    
    ft = dataFileNameToType(filename);
    return importlib.import_module("ClearMap.IO." + ft);


##############################################################################
# Data Sizes
##############################################################################

def arraySize(source, **args):
    """Returns size of the array source
       
    Arguments:
        source (array, tuple or str): source data or tuple with size specification
        
    Returns:
        tuple: size of the image data taking into account region specifications
    """ 
    
    if isinstance(source, str):
        mod = arrayFileNameToModule(source);
        return mod.size(source, **args);
    elif isinstance(source, numpy.ndarray):
        return sizeFromRegion(source.shape, **args);
    elif isinstance(source, tuple):
        return sizeFromRegion(source, **args);
    else:
        raise RuntimeError("arraySize: argument not a string, tuple or array!");  
    
def imageSize(source, **args):
    """Returns size of the image
       
    Arguments:
        source (array, tuple or str): source data or tuple size specification
        
    Returns:
        tuple: size of the image data taking into account region specifications
    """ 

    if isinstance(source, str):
        mod = imageFileNameToModule(source);
        return mod.size(source, **args);
    elif isinstance(source, numpy.ndarray):
        return sizeFromRegion(source.shape, **args);
    elif isinstance(source, tuple):
        return sizeFromRegion(source, **args);
    else:
        raise RuntimeError("imageSize: argument not a string, tuple or array!");  

    
def dataSize(source, **args):
    """Returns size of the data
       
    Arguments:
        source (array or str): source data
        x,y,z (tuple or all): range specifications, ``all`` is full range
        
    Returns:
        tuple: size of the image data taking into account region specifications
    """ 
    
    try:
      size = imageSize(source, **args);
    except:
      try:
        size = arraySize(source, **args);
      except:
         raise RuntimeError("dataSize: cannot determine size of the source!");
    
    return size;


def size(source, **args):
    """Returns size of the source taking into account any region specifications
    
    Arguments:
        source (array or str): source data
        x,y,z (tuple or all): range specifications, ``all`` is full range
        
    Returns:
        tuple: size of the image data taking into account region specifications
    """ 
      
    return dataSize(source, **args);


def arraySizeZ(source, **args):
    """Returns size of the array in the third dimension, None if array dimension < 3 
           
    Arguments:
        source (array or str): source data
        z (tuple or all): z-range specification, ``all`` is full range
        
    Returns:
        int: size of the image data in z after reading and range reduction
    """ 
      
    if isinstance(source, str):
        mod = arrayFileNameToModule(source);
        return mod.sizeZ(source, **args);
    elif isinstance(source, numpy.ndarray):
        size = dataSize(source.shape, **args);
    elif isinstance(source, tuple):
        size = dataSize(source, **args);
    else:
        raise RuntimeError("arraySizeZ: argument not a string, tuple or array!");        

    if len(size) > 2:
      return size[2];
    else:
      return None;
    

def imageSizeZ(source, **args):
    """Returns size of the image in the third dimension, None if image dimension < 3 
           
    Arguments:
        source (array or str): source data
        z (tuple or all): z-range specification, ``all`` is full range
        
    Returns:
        int: size of the image data in z after reading and range reduction
    """ 
      
    if isinstance(source, str):
        mod = imageFileNameToModule(source);
        return mod.sizeZ(source, **args);
    elif isinstance(source, numpy.ndarray):
        size = dataSize(source.shape, **args);
    elif isinstance(source, tuple):
        size = dataSize(source, **args);
    else:
        raise RuntimeError("imageSizeZ: argument not a string, tuple or array!");        

    if len(size) > 2:
      return size[2];
    else:
      return None;
        
        
def dataSizeZ(source, **args):
    """Returns size of the data in the third dimension, None if data dimension < 3 
           
    Arguments:
        source (array or str): source data
        z (tuple or all): z-range specification, ``all`` is full range
        
    Returns:
        int: size of the image data in z after reading and range reduction
    """ 
   
    try:
      size = imageSizeZ(source, **args);
    except:
      try:
        size = arraySizeZ(source, **args);
      except:
         raise RuntimeError("dataSizeZ: cannot determine size of the source!");
    
    return size;


def sizeZ(source, **args):
    """Returns size of the source in the third dimension, None if data dimension < 3 
           
    Arguments:
        source (array or str): source data
        z (tuple or all): z-range specification, ``all`` is full range
        
    Returns:
        int: size of the image data in z after reading and range reduction
    """ 
      
    return dataSizeZ(source, **args);



##############################################################################
# Read / Write Images
##############################################################################


def readImage(source, **args):
    """Read image from one of the supported formats
    
    Arguments:
        source (str, array or None): full data array, if numpy array simply reduce its range
        region (tuple or Region): region specifications
        x,z,y (tuple or all): specific region specifications, ``all`` is full range
        **args: further arguments specific to image data reader
    
    Returns:
        array: image data as numpy array, None if source is None
    
    See Also:
        :func:`writeImage`, :func:`readArray`, :func:`readData`, :func:`read`
    """  
    
    if source is None:
        return None;   
    elif isinstance(source, str):
        mod = imageFileNameToModule(source);
        return mod.read(source, **args);
    elif isinstance(source, numpy.ndarray ):
        return dataFromRegion(source, **args);
    else:
        raise RuntimeError('readImage: cannot infer format of the requested data/file.');


def writeImage(sink, data, **args):
    """Write data to one of the supported image formats
    
    Arguments:
        sink (str, array or None): the destination for the data, if None the data is returned directly
        data (array or None): data to be written
        **args: further arguments specific to image data format writer
    
    Returns:
        array, str or None: data or file name of the written data
    
    See Also:
        :func:`readImage`, :func:`writeArray`, :func:`writeData`, :func:`write`
    """ 
    
    if sink is None: #dont write to disk but return the data
        return data; #shallow copy
    elif isinstance(sink, numpy.ndarray):
        sink[:] = data[:]; #deep copy !
        return sink;
    elif isinstance(sink, str):
        mod = imageFileNameToModule(sink);    
        return mod.write(sink, data, **args);
    else:
        raise RuntimeError('writeImage: cannot infer format of the requested image sink to write to.');


def copyImage(source, sink):
    """Copy an image file from source to sink, which can consist of multiple files
    
    Arguments:
        source (str): file name of source
        sink (str): file name of sink
    
    Returns:
        str: name of the copied file
    
    See Also:
        :func:`copyArray`, :func:`copyData`, :func:`copy`, :func:`convertData`
    """     
    
    mod = imageFileNameToModule(source);
    return mod.copy(source, sink);

 
def combineImage(*args):
    """Concatenate single channel arrays to one multi channel array with channels in last dim
    
    Arguments:
        *args (arrays): arrays to be concatenated
    
    Returns:
        array: concatenated multi-channel array
    """
    
    data = numpy.array(args);
    return data.rollaxis(data, 0, data.ndim);



##############################################################################
# Read / Write Arrays
##############################################################################


def readArray(source, **args):
    """Read array from one of the supported array formats
    
    Arguments:
        source (str, array or None): full data array, if numpy array simply reduce its range
        region (tuple or Region): region specifications
        x,z,y (tuple or all): specific region specifications, ``all`` is full range
        **args: further arguments specific to image data reader
    
    Returns:
        array: array data as numpy array, None if source is None
    
    See Also:
        :func:`writeArray`, :func:`readImage`, :func:`readData`, :func:`read`
    """  
    
    if source is None:
        return None;   
    elif isinstance(source, str):
        mod = arrayFileNameToModule(source);
        return mod.read(source, **args);
    elif isinstance(source, numpy.ndarray ):
        return dataFromRegion(source, **args);
    else:
        raise RuntimeError('readArray: cannot infer format of the requested data/file.');


def writeArray(sink, data, **args):
    """Write data to one of the supported array formats
    
    Arguments:
        sink (str, array or None): the destination for the data, if None the data is returned directly
        data (array or None): data to be written
        **args: further arguments specific to image data format writer
    
    Returns:
        array, str or None: data or file name of the written data
    
    See Also:
        :func:`readArray`, :func:`writeImage`, :func:`writeData`, :func:`write`
    """ 
    
    if sink is None: #dont write to disk but return the data
        return data; #shallow copy
    elif isinstance(sink, numpy.ndarray):
        sink[:] = data[:]; #deep copy !
        return sink;
    elif isinstance(sink, str):
        mod = arrayFileNameToModule(sink);    
        return mod.write(sink, data, **args);
    else:
        raise RuntimeError('writeArray: cannot infer format of the requested array sink to write to.');


def copyArray(source, sink):
    """Copy an array file from source to sink, which can consist of multiple files
    
    Arguments:
        source (str): file name of source
        sink (str): file name of sink
    
    Returns:
        str: name of the copied file
    
    See Also:
        :func:`copyImage`, :func:`copyData`, :func:`copy`, :func:`convertData`
    """     
    
    mod = arrayFileNameToModule(source);
    return mod.copy(source, sink);
 
def combineArray(*args):
    """Concatenate single channel arrays to one multi channel array with channels in last dim
    
    Arguments:
        *args (arrays): arrays to be concatenated
    
    Returns:
        array: concatenated multi-channel array
    """
    
    data = numpy.array(args);
    return data.rollaxis(data, 0, data.ndim);


##############################################################################
# Read / Write Data
##############################################################################


def readData(source, **args):
    """Read data from one of the supported data formats
    
    Arguments:
        source (str, array or None): full data array, if numpy array simply reduce its range
        region (tuple or Region): region specifications
        x,z,y (tuple or all): specific region specifications, ``all`` is full range
        **args: further arguments specific to image data reader
    
    Returns:
        array: data as numpy array, None if source is None
    
    See Also:
        :func:`writeData`, :func:`readImage`, :func:`readArray`, :func:`read`
    """  
    
    if source is None:
        return None;   
    elif isinstance(source, str):
        mod = dataFileNameToModule(source);
        return mod.read(source, **args);
    elif isinstance(source, numpy.ndarray ):
        return dataFromRegion(source, **args);
    else:
        raise RuntimeError('readData: cannot infer format of the requested data/file.');


def writeData(sink, data, **args):
    """Write data to one of the supported data formats
    
    Arguments:
        sink (str, array or None): the destination for the data, if None the data is returned directly
        data (array or None): data to be written
        **args: further arguments specific to image data format writer
    
    Returns:
        array, str or None: data or file name of the written data
    
    See Also:
        :func:`readData`, :func:`writeImage`, :func:`writeArray`, :func:`write`
    """ 
    
    if sink is None: #dont write to disk but return the data
        return data; #shallow copy
    elif isinstance(sink, numpy.ndarray):
        sink[:] = data[:]; #deep copy !
        return sink;
    elif isinstance(sink, str):
        mod = dataFileNameToModule(sink);    
        return mod.write(sink, data, **args);
    else:
        raise RuntimeError('writeData: cannot infer format of the requested data sink to write to.');


def copyData(source, sink):
    """Copy a data file from source to sink, which can consist of multiple files
    
    Arguments:
        source (str): file name of source
        sink (str): file name of sink
    
    Returns:
        str: name of the copied file
    
    See Also:
        :func:`copyImage`, :func:`copyArray`, :func:`copy`, :func:`convertData`
    """     
    
    mod = dataFileNameToModule(source);
    return mod.copy(source, sink);
 
def combineData(*args):
    """Concatenate single channel data to one multi channel array with channels in last dim
    
    Arguments:
        *args (arrays): arrays to be concatenated
    
    Returns:
        array: concatenated multi-channel array
    """
    
    data = numpy.array(args);
    return data.rollaxis(data, 0, data.ndim);


##############################################################################
# Read / Write Generic
##############################################################################


def read(source, **args):
    """Read data from one of the supported data formats
    
    Arguments:
        source (str, array or None): full data array, if numpy array simply reduce its range
        region (tuple or Region): region specifications
        x,z,y (tuple or all): specific region specifications, ``all`` is full range
        **args: further arguments specific to image data reader
    
    Returns:
        array: data as numpy array, None if source is None
    
    See Also:
        :func:`write`, :func:`readImage`, :func:`readArray`, :func:`readData`
    """  
    
    return readData(source, **args);


def write(sink, data, **args):
    """Write data to one of the supported data formats
    
    Arguments:
        sink (str, array or None): the destination for the data, if None the data is returned directly
        data (array or None): data to be written
        **args: further arguments specific to image data format writer
    
    Returns:
        array, str or None: data or file name of the written data
    
    See Also:
        :func:`readData`, :func:`writeArray`, :func:`writeImage`, :func:`writeData`
    """ 
    
    return writeData(sink, data, **args);


def copy(source, sink):
    """Copy a data file from source to sink, which can consist of multiple files
    
    Arguments:
        source (str): file name of source
        sink (str): file name of sink
    
    Returns:
        str: name of the copied file
    
    See Also:
        :func:`copyImage`, :func:`copyArray`, :func:`copyData`, :func:`convert`
    """     
    
    return copyData(source, sink);


def convert(source, sink, **args):
    """Transforms data from source format to sink format
    
    Arguments:
        source (str): file name of source
        sink (str): file name of sink
    
    Returns:
        str: name of the copied file
        
    Warning:
        Not optimized for large image data sets yet
    
    See Also:
        :func:`copyImage`, :func:`combineImage`
    """      

    if source is None:
        return None;
    
    elif isinstance(source, str):
        if sink is None:        
            return read(source, **args);
        elif isinstance(sink, str):
            #if args == {} and dataFileNameToType(source) == dataFileNameToType(sink):
            #    return copy(source, sink);
            #else:
            data = read(source, **args); #TODO: improve for large data sets
            return write(sink, data);
        else:
            raise RuntimeError('convert: unknown sink format!');
            
    elif isinstance(source, numpy.ndarray):
        if sink is None:
            return dataFromRegion(source, **args);
        elif isinstance(sink,  str):
            data = dataFromRange(source, **args);
            return writeData(sink, data);
        else:
            raise RuntimeError('convert: unknown sink format!');
    
    else:
      raise RuntimeError('convert: unknown srouce format!');
    

##############################################################################
# Other
##############################################################################

def writeTable(filename, table):
    """Writes a numpy array with column names to a csv file.
    
    Arguments:
        filename (str): filename to save table to
        table (annotated array): table to write to file
        
    Returns:
        str: file name
    """
    with open(filename,'w') as f:
        for sublist in table:
            f.write(', '.join([str(item) for item in sublist]));
            f.write('\n');
        f.close();

    return filename;
