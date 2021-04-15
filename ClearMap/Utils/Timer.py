# -*- coding: utf-8 -*-
"""
Module provides tools for timing information

Example:
  >>> import ClearMap.Utils.Timer as timer
  >>> t = timer.Timer();
  >>> for i in range(100000000):
  >>>   x = 10 + i;
  >>> t.printElapsedTime('test')
"""
__author__    = 'Christoph Kirst <ckirst@rockefeller.edu>'
__license__   = 'MIT License <http://www.opensource.org/licenses/mit-license.php>'
__copyright__ = 'Copyright (c) 2017 by Christoph Kirst, The Rockefeller University, New York City'


import time

import ClearMap.Utils.Sound as snd;


class Timer(object):
  """Class to stop time and print results in formatted way
  
  Attributes:
      time (float): the time since the timer was started
  """
  def __init__(self, verbose=False):
    self.verbose = verbose;
    self.start();

  def start(self):
    """Start the timer"""
    self.time = time.time();

  def reset(self):
    """Reset the timer"""
    self.time = time.time();

  def elapsedTime(self, head = None, asstring = True):
    """Calculate elapsed time and return as formated string
    
    Arguments:
        head (str or None): prefix to the string
        asstring (bool): return as string or float
    
    Returns:
        str or float: elapsed time
    """
    
    t = time.time();
    
    if asstring:
      t = self.formatElapsedTime(t - self.time);
      if head != None:
        return head + ": elapsed time: " + t;
      else:
        return "Elapsed time: " + t;
    else:
      return t - self.time;
  
  def printElapsedTime(self, head = None, beep = False):
    """Print elapsed time as formated string
    
    Arguments:
        head (str or None): prefix to the string
    """
    print(self.elapsedTime(head = head));
    if beep:
      snd.beep();
  
  def formatElapsedTime(self, t):
    """Format time to string
    
    Arguments:
        t (float): time in seconds prefix
    
    Returns:
        str: time as hours:minutes:seconds
    """
    m, s = divmod(t, 60);
    h, m = divmod(m, 60);
    ms = 1000 * (s % 1);
    
    return "%d:%02d:%02d.%03d" % (h, m, s, ms);


def timeit(method):
  def timed(*args, **kw):
      ts = time.time()
      result = method(*args, **kw)
      te = time.time()
      
      m, s = divmod(te-ts, 60);
      h, m = divmod(m, 60);
      ms = 1000 * (s % 1);
      
      print ("%r took %d:%02d:%02d.%03d" % (method.__name__, h, m, s, ms))
      return result

  return timed


if __name__ == "__main__":
  
  import ClearMap.Utils.Timer as timer
  t = timer.Timer();
  for i in range(100000000):
    x = 10 + i;
  t.printElapsedTime('test')