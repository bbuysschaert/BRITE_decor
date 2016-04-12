# -*- coding: utf-8 -*-
"""
Functions you might want to use, while using the bounded least-squares minimisation with the lmfit package
    
Last update 16 March 2016

@author: Bram Buysschaert
"""

#===============================================================================
# 				Packages
#===============================================================================
import numpy as np
#===============================================================================
# 				Code
#===============================================================================
def lmfit_sinslope(params, x):
    """
    Function to be used for the lmfit LS minimisation routine to fit a given function to the data. The function here is a sine function with a *linear* background.
    
    @param param: lmfit parameter containing the information needed to fit
    @type time: lmfit parameter class
    @param x: parameter for which we want to determine a function; y=param(x)
    @type x: numpy array of length N
    
    Returns: The values of the function lmfit_sinslope.
    """
    frequency = params['frequency'].value
    constant = params['constant'].value
    phase = params['phase'].value
    amplitude = params['amplitude'].value
    slope = params['slope'].value

    return constant + slope * x +  amplitude * np.sin(2. * np.pi * (x * (frequency) + phase))

def lmfit_sinslope_vs_data(params, x, signal):
    """
    Function comparing the data with the returned values of the lmfit_sin function, which is a sine function with a *linear* background.
    
    @param param: lmfit parameter containing the information needed to fit
    @type time: lmfit parameter class
    @param x: parameter for which we want to determine a function; y=param(x)
    @type x: numpy array of length N
    @param signal: parameter which we want to represent with the function y=param(x)
    @type signal: numpy array of length N
    
    Returns: Residuals between the data and the values of the function lmfit_sinslope
    """
    sin_vs_data = signal - lmfit_sinslope(params, x)

    return sin_vs_data
      
def lmfit_sin(params, x):
    """
    Function to be used for the lmfit LS minimisation routine to fit a given function to the data. The function here is a sine function with a *constant* background.
    
    @param param: lmfit parameter containing the information needed to fit
    @type time: lmfit parameter class
    @param x: parameter for which we want to determine a function; y=param(x)
    @type x: numpy array of length N
    
    Returns: The values of the function lmfit_sin.
    """
    frequency = params['frequency'].value
    constant = params['constant'].value
    phase = params['phase'].value
    amplitude = params['amplitude'].value

    return constant + amplitude * np.sin(2. * np.pi * (x * (frequency) + phase))

def lmfit_sin_vs_data(params, x, signal):
    """
    Function comparing the data with the returned values of the lmfit_sin function, which is a sine function with a *constant* background.
    
    @param param: lmfit parameter containing the information needed to fit
    @type time: lmfit parameter class
    @param x: parameter for which we want to determine a function; y=param(x)
    @type x: numpy array of length N
    @param signal: parameter which we want to represent with the function y=param(x)
    @type signal: numpy array of length N
    
    Returns: Residuals between the data and the values of the function lmfit_sin
    """
    sin_vs_data = signal - lmfit_sin(params, x)

    return sin_vs_data   