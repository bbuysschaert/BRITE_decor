# -*- coding: utf-8 -*-
"""
Routines to perform a 1D and nD autocorrelation of a given dataset
    
Last update 18 March 2016

@author: Bram Buysschaert
"""

#===============================================================================
# 				Packages
#===============================================================================
import numpy as np
#===============================================================================
# 				Code
#===============================================================================
def acf1D(x, **kwargs):
  """
  Determines the autocorrelation of a 1D array.
  
  Returns: the autocorrelation of x
  
  @param x: measurements of x [???]
  @type x: numpy array of length N
  
  @return acor: the autocorrelation of x
  @rtype acor: numpy array of length N
  """
  #Calculating the autocorrelation function
  n = len(x)
  variance = x.var()
  x = x-x.mean()
  r = np.correlate(x, x, mode = 'full')[-n:]
  acor = r/(variance*(np.arange(n, 0, -1)))
  
  return acor

def acfnD(x, **kwargs):
    """
    Determines the autocorrelation of a nD array.
    
    NOTE routine was found at http://stackoverflow.com/questions/4503325/autocorrelation-of-a-multidimensional-array-in-numpy
    
    Returns: the autocorrelation of x
    
    @param x: measurements of x [???]
    @type x: numpy nDarray of length NxM
    
    @return acor: the autocorrelation of x
    @rtype acor: numpy NDarray of length NxM    
    """

    # used for transposes
    t = np.roll(range(x.ndim), 1)

    # pairs of indexes
    # the first is for the autocorrelation array
    # the second is the shift
    ii = [list(enumerate(range(1, s - 1))) for s in x.shape]

    # initialize the resulting autocorrelation array
    acor = np.empty(shape=[len(s0) for s0 in ii])

    # iterate over all combinations of directional shifts
    for i in product(*ii):
        # extract the indexes for
        # the autocorrelation array 
        # and original array respectively
        i1, i2 = np.asarray(i).T

        x1 = x.copy()
        x2 = x.copy()

        for i0 in i2:
            # clip the unshifted array at the end
            x1 = x1[:-i0]
            # and the shifted array at the beginning
            x2 = x2[i0:]

            # prepare to do the same for 
            # the next axis
            x1 = x1.transpose(t)
            x2 = x2.transpose(t)

        # normalize shifted and unshifted arrays
        x1 -= x1.mean()
        x1 /= x1.std()
        x2 -= x2.mean()
        x2 /= x2.std()

        # compute the autocorrelation directly
        # from the definition
        acor[tuple(i1)] = (x1 * x2).mean()

    return acor