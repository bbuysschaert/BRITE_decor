# -*- coding: utf-8 -*-
"""
Some general fitting routines for lmfit, which could provide to be useful. For example, when trying to determine interpixel variations.
    
Last update 17 March 2016

@author: Bram Buysschaert
"""

#===============================================================================
# 				Packages
#===============================================================================
import numpy as np

from lmfit import minimize, Parameters, Parameter, report_fit, conf_interval2d, conf_interval, report_ci

import BRITE.fitfunctions.ff_lmfit
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

def sineFIT(param1, param2, **kwargs):
    """
    Bounded LS minimisation for the fitting of a sinefunction to a given dataset. It is assumed that you can represent param2 as a function of param1. The lmfit package is used to perform the bounded LS minimisation.
    
    DANGER this routine assumes you have a fixed frequency of 1 [unit**-1]. For example:
    - if param1 is time and param2 is flux, you will have a sine with a frequency of 1 c/d. 
    - if param1 is position and param2 is flux, you will have a sine with a frequency of 1 c/pix.
    
    TODO provide options through kwargs to set the boundaries
    
    Returns: The optimal lmfit parameter class.
    
    @param param1: param1 measurements [???]
    @type param1: numpy array of length N
    @param param2: param2 measurements [???]
    @type param2: numpy array of length N
    
    @return: paramSINE
    @rtype: lmfit parameter class
    
    @kwargs: show_ME: Boolean to indicate if you want to see the report_fit - Default is False [Boolean]
    """
    show_ME = kwargs.get('show_ME', False)
    
    # Determination of the guesses to start your bounded LS fit with.
    constantGUESS = np.median(param2) 				# param2
    amplitudeGUESS = np.max(np.abs(constantGUESS-param2))/2. 	# param2
    frequencyGUESS = 1. 					# DANGER Here is the frequency assumption assumption.
    phaseGUESS = 0.1 						# Using the param1[np.where(param2==np.max(param2))]-param1[0]%1-0.5 is best, when there is no scatter on param2
    
    paramSINE = Parameters()
    #Make a params class for lmfit. 
    #		  	(Name,		Value,		Vary,	Min,				Max,				Expr)
    paramSINE.add_many(('amplitude',	amplitudeGUESS,	True,	amplitudeGUESS*0.1,		amplitudeGUESS*1.2,		None),
		      ('frequency',	frequencyGUESS,	False,	frequencyGUESS-0.05,		frequencyGUESS+0.05,		None), # DANGER Here is the frequency assumption assumption. (It is set to non-vary.)
		      ('constant',	constantGUESS,	True,	-abs(constantGUESS)*1.5,	abs(constantGUESS)*1.5,		None),
		      ('phase',		phaseGUESS,	True,	0.,				1.,				None))
    resultSIN = minimize(BRITE.fitfunctions.ff_lmfit.lmfit_sin_vs_data, paramSINE, args=(param1, param2))
    
    if show_ME:
      report_fit(paramsSIN, show_correl=False)
    
    return paramSINE

def sineFITwLINEARbackground(param1, param2, **kwargs):
    """
    Bounded LS minimisation for the fitting of a sinefunction to a given dataset *with a linear background instead of a constant background*. It is assumed that you can represent param2 as a function of param1. The lmfit package is used to perform the bounded LS minimisation.
    
    DANGER this routine assumes you have a fixed frequency of 1 [unit**-1]. For example:
    - if param1 is time and param2 is flux, you will have a sine with a frequency of 1 c/d. 
    - if param1 is position and param2 is flux, you will have a sine with a frequency of 1 c/pix.
    
    Returns: The optimal lmfit parameter class.
    
    TODO provide options through kwargs to set the boundaries
    
    Returns: The optimal lmfit parameter class.
    
    @param param1: param1 measurements [???]
    @type param1: numpy array of length N
    @param param2: param2 measurements [???]
    @type param2: numpy array of length N
    
    @return: paramSINE
    @rtype: lmfit parameter class
    
    @kwargs: show_ME: Boolean to indicate if you want to see the report_fit - Default is False [Boolean]
    """
    show_ME = kwargs.get('show_ME', False)
    # Determination of the guesses to start your bounded LS fit with.
    constantGUESS = np.median(param2) 				# param2
    amplitudeGUESS = np.max(np.abs(constantGUESS-param2))/2. 	# param2
    frequencyGUESS = 1. 					# DANGER Here is the frequency assumption assumption.
    phaseGUESS = 0.1 						# Using the param1[np.where(param2==np.max(param2))]-param1[0]%1-0.5 is best, when there is no scatter on param2
    slopeGUESS, constantGUESS = np.polyfit(param1, param2, 1)
    
    paramSINE = Parameters()
    #Make a params class for lmfit. 
    #		  	(Name,		Value,		Vary,	Min,				Max,				Expr)
    paramSINE.add_many(('amplitude',	amplitudeGUESS,	True,	amplitudeGUESS*0.1,		amplitudeGUESS*1.2,		None),
		      ('frequency',	frequencyGUESS,	False,	frequencyGUESS-0.05,		frequencyGUESS+0.05,		None), # DANGER Here is the frequency assumption assumption. (It is set to non-vary.)
		      ('constant',	constantGUESS,	True,	-abs(constantGUESS)*1.5,	abs(constantGUESS)*1.5,		None),
		      ('phase',		phaseGUESS,	True,	0.,				1.,				None),
		      ('slope',		slopeGUESS,	True,	-2*slopeGUESS,			+2*slopeGUESS,			None))
    resultSIN = minimize(BRITE.fitfunctions.ff_lmfitlmfit_sinslope_vs_data, paramSINE, args=(param1, param2))
    if show_ME:
      report_fit(paramsSIN, show_correl=False)
    
    return paramSINE  