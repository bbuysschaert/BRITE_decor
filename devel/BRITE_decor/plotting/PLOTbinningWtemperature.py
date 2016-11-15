# -*- coding: utf-8 -*-
"""
Routines to visualise the binning of the timeseries for the temperature dependent PSF changes, which produce instrumental effects between flux and position.
    
Last update 15 November 2016

@author: Bram Buysschaert
"""

#===============================================================================
# 				Packages
#===============================================================================
import pylab as pl
import numpy as np

import scipy.interpolate as scInterp

from matplotlib.ticker import NullFormatter, MultipleLocator, FormatStrFormatter, NullLocator
from matplotlib import patches, rc, rcParams, pyplot
import matplotlib.gridspec as gridspec

import statsmodels.api as sm
#===============================================================================
# 				Settings					
#===============================================================================
from string import whitespace as ws
import re
rc('xtick', labelsize=20)
rc('ytick', labelsize=20) 
rc("lines", markeredgewidth=2.0)
rc("axes", linewidth=2.0)
rcParams["font.size"] = 25

# Enabling a greek 'nu' for frequency, if necessary.
p = re.compile('(%s)' % ('|'.join([c for c in ws]))) 
s = " \ nu "
nu = p.sub('',s)
#===============================================================================
# 			Class for colored console printing			
#===============================================================================
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

def PLOTopenTIMEwithTEMPERATURE(time, flux, temperature, temperatureLONG, xPOS, yPOS, timeINDEXES, timeINDEXES_reason, **kwargs):
    """
    Routine going with openTIMEwithTEMPERATURE to visualize the binning process needed to detrend the temperature dependent PSF changes.
    
    Returns: Figure with -hopefully- enough diagnostics to determine the strength / weaknessess of the detrending
    
    @param time: time measurements [d]
    @type time: numpy array of length N
    @param flux: flux measurements [adu]
    @type flux: numpy array of length N
    @param temperature: temperature measurements [deg]
    @type temperature: numpy array of length N
    @param temperatureLONG: the longterm trend of the CCD temperature [deg]
    @type temperatureLONG: numpy array of length N
    @param xPOS: CCD x position measurements [pixel]
    @type xPOS: numpy array of length N
    @param yPOS: CCD y position measurements [pixel]
    @type yPOS: numpy array of length N
    @param timeINDEXES: end index (of your time array for that bin)
    @type timeINDEXES: numpy array (dtype='int32') of length K
    @param timeINDEXES_reason: reason to end your temperature bin
    @type timeINDEXES_reason: list of length K (containing strings)
    
    """
    figBINNING = pl.figure(figsize=(16,16))
    gsBINNING = gridspec.GridSpec(4, 1,height_ratios=[3,1,1,1])
    axTEMP = figBINNING.add_subplot(gsBINNING[0,0]) 
    axFLUX = figBINNING.add_subplot(gsBINNING[1,0], sharex=axTEMP) 
    axPOSX = figBINNING.add_subplot(gsBINNING[2,0], sharex=axTEMP) 
    axPOSY = figBINNING.add_subplot(gsBINNING[3,0], sharex=axTEMP) 
    
    axTEMP.plot(time, temperature, 'k.', ms=6, alpha=.4)
    try:
      axTEMP.plot(time, temperatureLONG, 'r-', lw=3, alpha=.6)
    except:
      lowessFRAC = kwargs.get('LOWESSfrac', 0.2) # []
      lowessLONG = sm.nonparametric.lowess(temperature, time, frac=lowessFRAC, delta=0.1)
      timeLONG = lowessLONG[:,0]
      axTEMP.plot(timeLONG, temperatureLONG, 'r-', lw=3, alpha=.6)      
    pl.tick_params('both',length=10,width=2,which='major'); pl.tick_params('both',length=10,width=1,which='minor')
    pyplot.locator_params(axis = 'x', nbins = 5); pyplot.locator_params(axis = 'y', nbins = 5)  
    axFLUX.plot(time, flux, 'k.', ms=6, alpha=.4)
    pl.tick_params('both',length=10,width=2,which='major'); pl.tick_params('both',length=10,width=1,which='minor')
    pyplot.locator_params(axis = 'x', nbins = 5); pyplot.locator_params(axis = 'y', nbins = 5)  
    axPOSX.plot(time, xPOS, 'k.', ms=6, alpha=.4)
    pl.tick_params('both',length=10,width=2,which='major'); pl.tick_params('both',length=10,width=1,which='minor')
    pyplot.locator_params(axis = 'x', nbins = 5); pyplot.locator_params(axis = 'y', nbins = 5)  
    axPOSY.plot(time, yPOS, 'k.', ms=6, alpha=.4)
    pl.tick_params('both',length=10,width=2,which='major'); pl.tick_params('both',length=10,width=1,which='minor')
    pyplot.locator_params(axis = 'x', nbins = 5); pyplot.locator_params(axis = 'y', nbins = 5)  
    for TT in range(len(timeINDEXES)):
      if timeINDEXES_reason[TT] == 'temperature':
	axTEMP.axvline(x=time[timeINDEXES[TT]], color='g', lw=5, alpha=.4); axFLUX.axvline(x=time[timeINDEXES[TT]], color='g', lw=5, alpha=.4)
	axPOSX.axvline(x=time[timeINDEXES[TT]], color='g', lw=5, alpha=.4); axPOSY.axvline(x=time[timeINDEXES[TT]], color='g', lw=5, alpha=.4)
      elif timeINDEXES_reason[TT] == 'gap':
	axTEMP.axvline(x=time[timeINDEXES[TT]], color='b', lw=5, alpha=.4); axFLUX.axvline(x=time[timeINDEXES[TT]], color='b', lw=5, alpha=.4)
	axPOSX.axvline(x=time[timeINDEXES[TT]], color='b', lw=5, alpha=.4); axPOSY.axvline(x=time[timeINDEXES[TT]], color='b', lw=5, alpha=.4)
      elif timeINDEXES_reason[TT] == 'end':
	axTEMP.axvline(x=time[timeINDEXES[TT]], color='c', lw=5, alpha=.4); axFLUX.axvline(x=time[timeINDEXES[TT]], color='c', lw=5, alpha=.4)
	axPOSX.axvline(x=time[timeINDEXES[TT]], color='c', lw=5, alpha=.4); axPOSY.axvline(x=time[timeINDEXES[TT]], color='c', lw=5, alpha=.4)
    axTEMP.set_ylabel('Temp [C]'); axFLUX.set_ylabel('Flux [adu]'); axPOSX.set_ylabel('xPOS [pixel]'); axPOSY.set_ylabel('yPOS')
    axPOSY.set_xlabel('Time [d]'); axPOSY.set_xlim([np.min(time)-0.2, np.max(time)+0.2])
