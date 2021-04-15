# -*- coding: utf-8 -*-
"""
Module provides tools for color manipulation
"""
__author__    = 'Christoph Kirst <ckirst@rockefeller.edu>'
__license__   = 'MIT License <http://www.opensource.org/licenses/mit-license.php>'
__copyright__ = 'Copyright (c) 2017 by Christoph Kirst, The Rockefeller University, New York City'


def hex_to_rgb(color):
  """Return rgb color from hex"""
  c = color.lstrip('#');  
  return [int(c[i:i+2], 16) for i in (0, 2 ,4)];


def rgb_to_hex(red, green, blue):
  """Return hex color from rgb."""
  return '#%02x%02x%02x' % (red, green, blue)
  
  
def writePAL(filename, cols):
  """Write a pal pallette file from a list of colors for use with e.g. imaris"""
  with open(filename, 'w') as f:
    for c in cols:
      f.write(('%3.3f' % c[0]).rjust(10) +  ('%3.3f' % c[1]).rjust(11) + ('%3.3f' % c[2]).rjust(11));
      f.write('\n');
    f.close();

def writeLUT(filename, cols):
  """Write a lut lookup table file from a list of colors for use with e.g. imagej"""
  with open(filename, 'w') as f:
    for c in cols:
      f.write(('%d' % c[0]).rjust(4) +  ('%d' % c[1]).rjust(5) + ('%d' % c[2]).rjust(5));
      f.write('\n');
    f.close();
 

if __name__ == "__main__":
  
  import ClearMap.Visualization.Color as color
  reload(color);
  
  color.hex_to_rgb('FFFFFF')