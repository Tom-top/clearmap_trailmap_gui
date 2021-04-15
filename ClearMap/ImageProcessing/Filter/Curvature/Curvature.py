# -*- coding: utf-8 -*-
"""
Module to calculate various curvature and tube measures in 3D
"""
__author__    = 'Christoph Kirst <ckirst@rockefeller.edu>'
__license__   = 'MIT License <http://www.opensource.org/licenses/mit-license.php>'
__copyright__ = 'Copyright (c) 2017 by Christoph Kirst, The Rockefeller University, New York City'

import numpy as np

def hessian(data):
  """Returns the hessian matrix at each location calculatd via finaite differences
  
  Arguments
  ---------
    data : ndarray
      the 3d image data
    

  Returns
  -------
    ndarray:
      5d array with the hessian matrix in the first two dimensions
  """    
  
  if data.ndim != 3:
    raise RuntimeError('hessian expexts data to be 3d image');
    
  mm = np.pad(data, 1, 'edge');
  h = np.zeros((3,3) + data.shape);
  
  c = 2*data;
  h[0,0] = mm[0:-2,1:-1,1:-1] - c + mm[2:,1:-1,1:-1];
  h[1,1] = mm[1:-1,0:-2,1:-1] - c + mm[1:-1,2:,1:-1];
  h[2,2] = mm[1:-1,1:-1,0:-2] - c + mm[1:-1,1:-1,2:];
  
  h[0,1] = h[1,0] = (mm[2:,2:,1:-1] - mm[:-2,2:,1:-1] - mm[2:,:-2,1:-1] + mm[:-2,:-2,1:-1]) / 4;
  h[0,2] = h[2,0] = (mm[2:,1:-1,2:] - mm[:-2,1:-1,2:] - mm[2:,1:-1,:-2] + mm[:-2,1:-1,:-2]) / 4;
  h[1,2] = h[2,1] = (mm[1:-1,2:,2:] - mm[1:-1,:-2,2:] - mm[1:-1,2:,:-2] + mm[1:-1,:-2,:-2]) / 4;

  return h;


def eigenvalues3D(m):
    """Find the coefficients of the characteristic polynomial:
       http://en.wikipedia.org/wiki/Eigenvalue_algorithm

		// The matrix looks like:
		/*
			A  B  C
			B  D  E
			C  E  F
		*/
    """
    
    A = m[0,0];
    B = m[0,1];
    C = m[0,2];
    D = m[1,1];
    E = m[1,2];
    F = m[2,2];

    a = -1;
    b = A + D + F;
    
    c = B * B + C * C + E * E - A * D - A * F - D * F;
    
    d = A * D * F - A * E * E	- B * B * F + 2 * B * C * E - C * C * D;

    third = 0.333333333333333333333333333333333333;

    #Now use the root-finding formula described here:
    #http://en.wikipedia.org/wiki/Cubic_equation#oot-finding_formula

    q = (3*a*c - b*b) / (9*a*a);
    r = (9*a*b*c - 27*a*a*d - 2*b*b*b) / (54*a*a*a);

    discriminant = q*q*q + r*r;
    #print discriminant.shape
    
    eigenValues = np.zeros((3,) + m.shape[2:]);
    #print eigenValues.shape

    #ids = discriminant > 0;
    #if np.any(ids):
        #should never happen for symmetric matrix
    #    return None;

    ids = discriminant < 0;
    
    rr = r[ids];
      
    rootThree = 1.7320508075688772935;
    innerSize = np.sqrt( rr*rr - discriminant[ids] );
      
    ids2 = rr > 0;
    innerAngle = np.zeros(rr.shape);
    
    innerAngle[ids2] = np.arctan(np.sqrt(-discriminant[ids][ids2]) / rr[ids2] );
    
    ids2 = np.logical_not(ids2);
    innerAngle[ids2] = ( np.pi - np.arctan( np.sqrt(-discriminant[ids][ids2]) / -rr[ids2] ) );
       
    # So now s is the cube root of innerSize * e ^ (   innerAngle * i )
    # and t is the cube root of innerSize * e ^ ( - innerAngle * i )       
       
    stSize = np.power(innerSize,third);
      
    sAngle = innerAngle / 3;
    #tAngle = - innerAngle / 3;
     
    sPlusT = 2 * stSize * np.cos(sAngle);
    
    eigenValues[0][ids] = ( sPlusT - (b[ids] / (3*a)) );
    firstPart = - (sPlusT / 2) - (b[ids] / 3*a);
    lastPart = - rootThree * stSize * np.sin(sAngle);
      
    eigenValues[1][ids] = ( firstPart + lastPart );
    eigenValues[2][ids] = ( firstPart - lastPart );

    
    #The discriminant is zero (or very small),
    #so the second two solutions are the same:
    ids = discriminant == 0;
    rr = r[ids];    
    
    ids2 = rr >= 0;
        
    
    sPlusT = np.zeros(rr.shape);
    sPlusT[ids2] = 2 * np.power(rr[ids2],third);
    ids2 = np.logical_not(ids2);
    sPlusT[ids2] = -2 * np.power(-rr[ids2],third);
        
    bOver3A = b[ids] / (3 * a);
      
    eigenValues[0][ids] = ( sPlusT   - bOver3A );
    eigenValues[1][ids] = (-sPlusT/2 - bOver3A );
    eigenValues[2][ids] = eigenValues[1][ids];
      
    return eigenValues;


from scipy import ndimage as ndi    


def eigenvalues(m, sigma = 1.0, sort = True):
  if sigma is not None:
    smoothed = ndi.gaussian_filter(np.asarray(m, dtype = float), sigma=sigma);
  else:
    smoothed = m;

  h3 = hessian(smoothed);
  e3 = eigenvalues3D(h3);
  e3a = np.abs(e3);
  
  if sort:
    index = list(np.ix_(*[np.arange(i) for i in e3a.shape]))
    index[0] = e3a.argsort(axis = 0)
    return e3[index];
  else:
    return e3;


def tubeness(m, sigma = 1.0):   
  if sigma is not None:
    smoothed = ndi.gaussian_filter(np.asarray(m, dtype = float), sigma=sigma);
  else:
    smoothed = m;

  h3 = hessian(smoothed);
  e3 = eigenvalues3D(h3);
  e3a = np.abs(e3);
  
  index = list(np.ix_(*[np.arange(i) for i in e3a.shape]))
  index[0] = e3a.argsort(axis = 0)
  e3s = e3[index];
  
  tb = np.zeros(m.shape);
  ids = np.logical_and(e3s[1] < 0, e3s[2] < 0);
  tb[ids] = np.sqrt(e3s[1][ids] * e3s[2][ids]);

  return tb;


def frangi(m, alpha = 0.5, beta  = 0.5, gamma  = 100, sigma = 1.0):
  if sigma is not None:
    smoothed = ndi.gaussian_filter(np.asarray(m, dtype = float), sigma=sigma);
  else:
    smoothed = m;

  h3 = hessian(smoothed);
  e3 = eigenvalues3D(h3);
  e3a = np.abs(e3);
  
  index = list(np.ix_(*[np.arange(i) for i in e3a.shape]))
  index[0] = e3a.argsort(axis = 0)
  e3s = e3[index];
  
  ra = np.abs(e3s[1]) / np.abs(e3s[2]);
  rb = np.abs(e3s[0]) / np.sqrt(np.abs(e3s[1] * e3s[2]));
  s  = np.sqrt(np.square(e3s[0]) + np.square(e3s[1]) + np.square(e3s[2]))

  plate = 1 - np.exp(- np.square(ra) / (2 * np.square(alpha)));
  blob  = np.exp(-np.square(rb) / (2 * np.square(beta)));
  background = 1 - np.exp(-np.square(s) / (2 * np.square(gamma)));

  f = plate * blob * background;
  
  ids = np.logical_and(e3s[1] >= 0, e3s[2] >= 0);
  f[ids] = 0;
  
  return f;

