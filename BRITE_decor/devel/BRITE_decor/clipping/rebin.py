# -*- coding: utf-8 -*-
"""
Routines to perform rebinning of data
    
Last update 18 March 2016

@author: Bram Buysschaert
"""

#===============================================================================
# 				Packages
#===============================================================================
import numpy as np

import scipy.interpolate as scInterp

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
def rebinUNIQUE(param1, param2, **kwargs):
    """
    Rebinning of a given param2, according to unique values in the param1 array.
    This should be a general script, which does not assume anything on both arrays, except that they have an equal length
    
    Returns: an array with the unique param1 values, a matrix with bin values (mean, median, std, percentiles) for param2
    
    @param param1: parameter 1 values (the values of which you take unique values)
    @type param1: numpy array of length N
    @param param2: parameter 2 values (the values of which you take bin values)
    @type param2: numpy array of length N
    
    @return: param1UNIQUE
    @rtype: numpy array of length M
    @return: param2MATRIX
    @rtype: numpy matrix of size 5xM
    """
    param1UNIQUE = np.unique(param1)
    param2MATRIX = np.zeros((5,len(param1UNIQUE))) #mean, median, std, 15.85%, 84.15%
    
    for uu in range(len(param1UNIQUE)):
      param2BIN = param2[param1==param1UNIQUE[uu]]
      param2MATRIX[:,uu] = [np.mean(param2BIN), np.median(param2BIN), np.std(param2BIN), np.percentile(param2BIN, 15.85), np.percentile(param2BIN, 84.15)]
    
    return param1UNIQUE, param2MATRIX