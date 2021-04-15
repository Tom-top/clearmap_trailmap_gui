# -*- coding: utf-8 -*-
"""
Module providing routines to check and convert between regular and file expressions 

Utility module for io modules :mod:`~ClearMap.IO`, :mod:`~ClearMap.IO.FileList`, and :mod:`~ClearMap.IO.FileArray`
"""
__author__    = 'Christoph Kirst <ckirst@rockefeller.edu>'
__license__   = 'MIT License <http://www.opensource.org/licenses/mit-license.php>'
__docformat__ = 'rest'


import fnmatch
import sre_parse as sre

import ClearMap.Utils.InverseRegularExpression as ire;


def isExpression(expression, groups = None, nPatterns = None, exclude = None, verbose = False):
  """Checks if the regular expression fullfill certain criteria
  
  Arguments:
    expression (str): regular expression to check
    groups (list or None): list of group names that should be present
    nPatterns (int or None): number of patterns to expect
    exclude (list or None): exculde these tokens
    verbose (bool): if True, print reason for expression to not fullfil desired criteria
    
  Returns
    bool: True if the expression fullfills the desired criteria
  """
  
  #parse regular expression 
  p = sre.parse(expression);

  #group patterns
  gd = p.pattern.groupdict
    
  if groups is None:
    groups = [];
  elif groups is all:
    groups = gd.keys();
  
  for gn in gd.keys():
    if gn not in groups:
      if verbose:
        print('Expression does contain a non required group %s!' % gn);
      return False;
  
  for gn in groups:
    if gn not in gd.keys():
      if verbose:
        print('Expression does not contain required group %s!' % gn);
      return False;

  if exclude is None:
    exclude = [];
  
  n = 0;
  for i,l in enumerate(p):
    lit = l[0];
    if lit != sre.LITERAL and lit not in exclude:
      n = n + 1;
  
  if nPatterns is not None:
    if nPatterns > 0:
      if n != nPatterns:
        if verbose:
          print('Expression has more than %d regular expression patterns!' % nPatterns);
        return False;
    
    elif nPatterns < 0:
      if n <= -nPatterns:
        if verbose:
           print('Expression has less than %d regular expression patterns!' % -nPatterns);
        return False;
    
    else:
      if n == 0:
        if verbose:
           print('Expression has no regular expression patterns!');
        return False;
  
  return True;


def expressionToPattern(expression):
  """Convert regular expression to a parsed pattern for manipulation
  
  Arguments:
    expression (str): regular expression
    
  Returns:
    object: parsed pattern
  """
  return sre.parse(expression);


def patternToExpression(pattern):
  """Convert pattern to regular expression
  
  Arguments:
    pattern (object): regular expression pattern
    
  Returns:
    str: regular expression
  """  
  ire.patternToExpression(pattern);


def replace(expression, replace = None):
  """Replaces patterns in a regular expression with given strings
  
  Arguments:
    expression (str): regular expression
    replace (dict): replacements of patterns given as entries pos : str or groupname : str
  
  Returns:
    str: converted regular expression
  """

  if not isinstance(replace, dict): 
    return expression;
  
  #parse regular expression 
  pattern = sre.parse(expression);

  #group patterns
  gd = pattern.pattern.groupdict
  gdi = {y:x for x,y in gd.items()};
  
  gkeys = gdi.keys();
  rkeys = replace.keys();
  
  newpattern = [];
  for p in pattern:
    if p[0] == sre.SUBPATTERN:
      sid = p[1][0];
      if sid is not None:
        if sid in rkeys:
          newpattern += [(ire.TEXT, replace[sid])];
        elif sid in gkeys:
          gn = gdi[sid];
          if gn in rkeys:
            newpattern += [(ire.TEXT, replace[gn])];
          else:
            newpattern.append(p);
        else:
          newpattern.append(p);
      else:
        newpattern.append(p);
    else:
      newpattern.append(p);
  
  pattern.data = newpattern;
  return ire.patternToExpression(pattern);
      

def insertGroups(expression, groups = None):
  """Inserts group names into a regular expression at specified sub patterns
  
  Arguments:
    expression (str): regular expression
    groups (dict or None): dictionary specifining the group names as groupid : groupname
    
  Returns:
    str: regular expression with named groups
  """
  
  if not isinstance(groups, dict):
    return expression;

  pattern = sre.parse(expression);
  
  #group patterns
  gd = pattern.pattern.groupdict
  for i,n in groups.items():
    gd[n] = i;
  
  pattern.pattern.groupdict = gd;
  return ire.patternToExpression(pattern);


def expressionToGlob(expression, replace = all, default = '*'):
  """Converts regular expression to a glob expression to search for files"""
  
  if not isinstance(replace, dict):
    replace = {};
  
  #parse regular expression 
  pattern = sre.parse(expression);

  #group patterns
  gd = pattern.pattern.groupdict
  gdi = {y:x for x,y in gd.items()};
  
  gkeys = gdi.keys();
  rkeys = replace.keys();
  
  defp = [(ire.TEXT, default)];
  
  newpattern = [];
  for p in pattern:
    if p[0] == sre.SUBPATTERN:
      sid = p[1][0];
      if sid is not None:
        if sid in rkeys:
          newpattern += [(ire.TEXT, replace[sid])];
        elif sid in gkeys:
          gn = gdi[sid];
          if gn in rkeys:
            newpattern += [(ire.TEXT, replace[gn])];
          else:
            newpattern += defp;
        else:
          newpattern += defp;
      else:
        newpattern += defp;
    else:
      newpattern.append(p);
  
  pattern.data = newpattern;
  return ire.patternToExpression(pattern);


def globToExpression(expression, groups = None):
  """Converts a glob expression to a regular expression
  
  Arguments:
    expression (str): glob expression
    groups (dict or None): a dicitonary specifying named groups in the form id: name
  
  Returns:
    str: regular expression
  """
  
  regex =  fnmatch.translate(expression);
  regex =  insertGroups(regex, groups = groups);
  if groups is None:
    regex = ire.patternToExpression(sre.parse(regex));
  return regex;



def test():  
  import ClearMap.Utils.RegularExpression as cre
  reload(cre);
  source = r'/path/to/files/image_file_(?P<row>\d{2})_(?P<col>\d{2})_(?P<z>\d{4}).tif'
  cre.isExpression(source, groups = ['row', 'col', 'z'], exclude = ['any'])

  cre.isExpression(source, groups = all, nPatterns=3, verbose = True)
  
  cre.isExpression(source, groups = all, nPatterns=3, verbose = True,  exclude = ['any'])
  
  cre.insertGroups(r'test/([7-8]).tif', groups = {1 : 'z'})
  
  cre.expressionToGlob(r'test/([7-8]).tif')
  
  cre.globToExpression('test/*.tif')

if __name__ == '__main__':
  test();