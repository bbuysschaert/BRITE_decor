# -*- coding: utf-8 -*-
"""
Routines to perform different median-sigma based outlier rejection.
    
Last update 17 March 2016

TODO update the routines, so you can do only upper / lower clipping.

NOTE: percentageFILTERonRELATION is much faster than medianFILTERonRELATION, since it only calculates one lowess filter

@author: Bram Buysschaert
"""

#===============================================================================
# 				Packages
#===============================================================================
import numpy as np

import statsmodels.api as sm
#===============================================================================#
# 			Class for colored console printing			#
#===============================================================================#
"""
Copied from http://stackoverflow.com/questions/22886353/printing-colors-in-python-terminal
"""
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
#===============================================================================
# 				Code
#===============================================================================
def medianFILTER(param, sigma, **kwargs):
    """
    Routine to perform a median filtering of param, with a threshold of sigma*median.
    
    NOTE: you could also just use scipy.stats.sigmaclip(param, sigma, sigma), but does this for CLIPiteration = np.inf and it uses the mean, not median
    
    NOTE: python sets are unstructured and unordered lists
    
    Returns: The indexes of the outliers, compared to the original array.   
    
    @param param: param measurements [???]
    @type param: numpy array of length N
    @param sigma: threshold for the rejection []
    @type sigma: numpy.float
    
    @return IDXoutliers: array of the indexes of the outliers
    @rtype: numpy array of length K (dtype='int32')
    
    @kwargs: CLIPiteration: maximum number of iterations you want to do the median clipping - Default is 100 ~np.inf [integer]
    """
    # Reading in the kwargs
    maxITER = np.int(kwargs.get('CLIPiteration', 100)) #[integer]
    
    # Doing a deep copy, so you have two arrays to work with. One to compare and on to refer too.
    paramCOPY = np.copy(param)
    
    # Performing the clipping itself
    iteration = 0
    while (iteration < maxITER) and (any(np.abs(paramCOPY - np.median(paramCOPY)) >= sigma) == True):
      IDXoutliersITERATION = np.where(np.abs(paramCOPY - np.median(paramCOPY)) >= sigma)[0] # [0] needed because of the set-up of np.where.
      paramCOPY = np.delete(paramCOPY, IDXoutliersITERATION)
      iteration +=1
      
    # Using sets to compare paramCOPY to param, find the removed elements, and lastly determine their original indexes. WARNING might not optimal for very long arrays.    
    paramREMOVED = list(set(param) - set(paramCOPY))
    
    IDXoutliers = set() # Empty set
    for rr in range(len(paramREMOVED)):
      IDXoutliersELEMENT = np.where(param == paramREMOVED[rr])[0] # Does not have to be ONE value, in case you have a non-unique entries in your array == headaches...
      IDXoutliers.update(set(list(IDXoutliersELEMENT))) # Append a new set of a list of a numpy array (easiest way to get around all lengths of arrays)
    
    # Converting the set back to a numpy array
    IDXoutliers = np.unique(np.array(list(IDXoutliers), dtype='int32'))
    return sorted(IDXoutliers)

def medianFILTERonRELATION(param1, param2, sigma, **kwargs):
    """
    Routine to perform a median filtering of param2, with a threshold of sigma*median, which also shows a clear trend with param1. We approximate this trend with a lowess filter (i.e. smoothing).
    
    NOTE: python sets are unstructured and unordered lists
    
    Returns: The indexes of the outliers, compared to the original array.   
    
    @param param1: param1 measurements [???]
    @type param1: numpy array of length N
    @param param2: param2 measurements [???]
    @type param2: numpy array of length N
    @param sigma: threshold for the rejection []
    @type sigma: numpy.float
    
    @return IDXoutliers: array of the indexes of the outliers
    @rtype: numpy array of length K (dtype='int32')
    
    @kwargs: CLIPiteration: maximum number of iterations you want to do the median clipping - Default is 100 ~np.inf [integer]
    @kwargs: LOWESSfrac: fraction of the data you wish to use for the lowess filter - default is 0.1 (should be in range ~0.15 and ~0.35)
    """
    # Reading in the kwargs
    maxITER = np.int(kwargs.get('CLIPiteration', 100)) #[integer]
    lowessFRAC = kwargs.get('LOWESSfrac', 0.1) # []
    # Doing a deep copy, so you have two arrays to work with. One to compare and on to refer too.
    param1COPY, param2COPY = np.copy(param1), np.copy(param2)
    
    # Performing the clipping itself. We use a lowess filter to smooth the relation with param1.
    iteration = 0
    
    lowessRELATION = sm.nonparametric.lowess(param2COPY, param1COPY, frac=lowessFRAC)
    param1RELATION, param2RELATION = lowessRELATION[:,0], lowessRELATION[:,1]
    
    while (iteration < maxITER) and (any(np.abs((param2COPY - param2RELATION) - np.median(param2COPY - param2RELATION)) >= sigma) == True):
      IDXoutliersITERATION = np.where(np.abs((param2COPY - param2RELATION) - np.median(param2COPY - param2RELATION)) >= sigma)[0] # [0] needed because of the set-up of np.where.
      
      param1COPY, param2COPY = np.delete(param1COPY, IDXoutliersITERATION), np.delete(param2COPY, IDXoutliersITERATION)
      
      lowessRELATION = sm.nonparametric.lowess(param2COPY, param1COPY, frac=lowessFRAC)
      param1RELATION, param2RELATION = lowessRELATION[:,0], lowessRELATION[:,1]
      iteration +=1
    
    print iteration
    # Using sets to compare paramCOPY to param, find the removed elements, and lastly determine their original indexes. WARNING might not optimal for very long arrays.    
    paramREMOVED = list(set(param2) - set(param2COPY))
    
    IDXoutliers = set() # Empty set
    for rr in range(len(paramREMOVED)):
      IDXoutliersELEMENT = np.where(param2 == paramREMOVED[rr])[0] # Does not have to be ONE value, in case you have a non-unique entries in your array == headaches...
      IDXoutliers.update(set(list(IDXoutliersELEMENT))) # Append a new set of a list of a numpy array (easiest way to get around all lengths of arrays)
    
    # Converting the set back to a numpy array
    IDXoutliers = np.unique(np.array(list(IDXoutliers), dtype='int32'))
    
    return sorted(IDXoutliers)
  
if __name__ == '__main__':
    # TESTING needed
    
    testARRAY = np.arange(0,10,0.1)
    print 'outliers are at indexes ', medianFILTER(testARRAY, 2.) 
    
    arrayIN2 = np.arange(10000)
    testARRAY2 = np.random.rand(10000)*5 + 5 + 0.1* arrayIN2
    
    print 'outliers are at indexes ', medianFILTERonRELATION(arrayIN2, testARRAY2, 2.) 
    
