# -*- coding: utf-8 -*-
"""
File Utils

This module provides utilities for file management used by various IO modules

See :mod:`ClearMap.IO`.
"""
__author__    = 'Christoph Kirst <ckirst@rockefeller.edu>'
__license__   = 'MIT License <http://www.opensource.org/licenses/mit-license.php>'

import os
import re
import shutil
import natsort

import ClearMap.Utils.RegularExpression as cre

##############################################################################
### Basic file queries
##############################################################################

def fileExtension(filename):
    """Returns file extension if exists
    
    Arguments:
        filename (str): file name
        
    Returns:
        str: file extension or None
    """
    
    if not isinstance(filename, str):
        return None;
    
    fext = filename.split('.');
    if len(fext) < 2:
        return None;
    else:
        return fext[-1];


def readFileList(filename, sort = True):
    """Returns the list of files that match the regular expression
    
    Arguments:
        filename (str): file name as regular expression
        sort (bool): sort files naturally
    
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
        raise RuntimeError('No files found in ' + fpath + ' matching ' + fname + ' !');
    
    #fl.sort();
    if sort:
      fl = natsort.natsorted(fl);
        
    return fpath, fl;


##############################################################################
### Checks on files
##############################################################################

    
def isFile(source):
    """Checks if source is an existing file, returns false if it is a directory
    
    Arguments:
        source (str): source file name
        
    Returns:
        bool: true if source is a real file   
    """
    
    if not isinstance(source, str):
        return False;
  
    if os.path.exists(source):
        if os.path.isdir(source):
            return False;
        else:
            return True;
    else:
        return False;

     
def isDirectory(source):
    """Checks if source is an exsiting directory
    
    Arguments:
        source (str): source file name
        
    Returns:
        bool: true if source is a real file   
    """
    
    if not isinstance(source, str):
        return False;
  
    if os.path.exists(source):
        if os.path.isdir(source):
            return True;
        else:
            return False;
    else:
        return False;
      

def isFileExpression(source, check = True):
    """Checks if source is an exsiting file expression with at least one file matching
    
    Arguments:
        source (str): source file name
        
    Returns:
        bool: true if source is a real file matching the file expression
    """
    
    if not isinstance(source, str):
        return False;
      
    if not cre.isExpression(source, nPatterns = 0, exclude = ['any']):
      return False
    
    if isFile(source):
      return False;

    if check:
      try:
        fp, fl = readFileList(source, sort = False);
      except:
        return False;
      return len(fl) > 0;
    else:
      return True;

    
def isFileName(source):
    """Checks if source is a file name and not a directory or file expression
    
    Arguments:
        source (str): source file name
        
    Returns:
        bool: true if source is a real file   
    """
    
    if not isinstance(source, str):
        return False;

    if os.path.exists(source):
        if os.path.isdir(source):
            return False;
        else:
            return True;
    else:
        if not isFileExpression(source):
          return True;
        else:
          return False;
    

##############################################################################
### File manipulation
##############################################################################

def createDirectory(filename, create = True, split = True):
    """Creates the directory of the file if it does not exists
     
    Arguments:
        filename (str): file name
        
    Returns:
        str: directory name
    """       
    if split:
      dirname, fname = os.path.split(filename);
    else:
      dirname = filename;    
    
    if create and not os.path.exists(dirname):
        os.makedirs(dirname);
    
    return dirname;

    
def copyFile(source, sink):
    """Copy a file from source to sink
    
    Arguments:
        source (str): file name of source
        sink (str): file name of sink
    
    Returns:
        str: name of the copied file
    
    See Also:
        :func:`ClearMap.IO.IO.copyData`, :func:`ClearMap.IO.IO.convertData`
    """ 
    
    shutil.copy(source, sink);
    return sink;



if __name__ == '__main__':
  import ClearMap.IO.FileUtils as fu
  import importlib
  importlib.reload(fu)