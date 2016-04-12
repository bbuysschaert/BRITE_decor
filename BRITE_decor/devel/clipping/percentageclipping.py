# -*- coding: utf-8 -*-
"""
Routines to perform different percentage based outlier rejection.

NOTE: percentageFILTERonRELATION is much faster than medianFILTERonRELATION, since it only calculates one lowess filter
    
Last update 17 March 2016

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
def percentageFILTER(param, percentageLOW, percentageUP, **kwargs):
    """
    Routine to perform a filtering of param, with a upper and lower threshold percentage.
    
    NOTE: python sets are unstructured and unordered lists
    
    Returns: The indexes of the outliers, compared to the original array.   
    
    @param param: param measurements [???]
    @type param: numpy array of length N
    @param percentageLOW: lower threshold percentage for the rejection [0-100]
    @type percentageLOW: numpy.float
    @param percentageUP: upper threshold percentage for the rejection [0-100]
    @type percentageUP: numpy.float
    
    @return IDXoutliers: array of the indexes of the outliers
    @rtype: numpy array of length K (dtype='int32')
    """
    
    # Determine the values of the upper and lower percentiles.
    lowerLIMIT, upperLIMIT = np.percentile(param2 - param2RELATION, [percentageLOW, 100-percentageUP])
    
    # numpy.where doesn't like nested statements. So, do it in two steps...
    IDXoutliers = np.where(param2 >= upperLIMIT)[0]
    IDXoutliers = np.append(IDXoutliers, np.where(param2 <= lowerLIMIT)[0])
    
    return sorted(IDXoutliers)

def percentageFILTERonRELATION(param1, param2, percentageLOW, percentageUP, **kwargs):
    """
    Routine to perform a filtering of param2, with a upper and lower threshold percentage, which also shows a clear trend with param1. We approximate this trend with a lowess filter (i.e. smoothing).
    
    NOTE: it is assumed that you choose LOWESSfrac, to be not influenced by the outliers. In case of a high LOWESSfrac and high percentages, this will pose a problem!
    
    
    Returns: The indexes of the outliers, compared to the original array.   
    
    @param param1: param1 measurements [???]
    @type param1: numpy array of length N
    @param param2: param2 measurements [???]
    @type param2: numpy array of length N
    @param percentageLOW: lower threshold percentage for the rejection [0-100]
    @type percentageLOW: numpy.float
    @param percentageUP: upper threshold percentage for the rejection [0-100]
    @type percentageUP: numpy.float
    
    @return IDXoutliers: array of the indexes of the outliers
    @rtype: numpy array of length K (dtype='int32')
    
    @kwargs: LOWESSfrac: fraction of the data you wish to use for the lowess filter - default is 0.1 (should be in range ~0.15 and ~0.35)
    """
    # Reading in the kwargs
    lowessFRAC = kwargs.get('LOWESSfrac', 0.1) # []
    
    # Doing the smoothing with the lowess filter.
    lowessRELATION = sm.nonparametric.lowess(param2, param1, frac=lowessFRAC)
    param1RELATION, param2RELATION = lowessRELATION[:,0], lowessRELATION[:,1]
    
    # Determine the values of the upper and lower percentiles.
    lowerLIMIT, upperLIMIT = np.percentile(param2 - param2RELATION, [percentageLOW, 100-percentageUP])
    # numpy.where doesn't like nested statements. So, do it in two steps...
    IDXoutliers = np.where((param2 - param2RELATION) >= upperLIMIT)[0]
    IDXoutliers = np.append(IDXoutliers, np.where((param2 - param2RELATION) <= lowerLIMIT)[0])
    
    return sorted(IDXoutliers)