# -*- coding: utf-8 -*-
"""
Routines to visualise the detrending for the position - flux instrumental effects.
    
Last update 21 March 2016

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
def PLOTdetrendPOSITIONfluxDIAGinformCRIT(flux, position, AICtck, BICtck, **kwargs):
    """
    Routine going with detrendPOSITIONflux to visualize diagnostics in case the AIC and BIC do not favor the same fit.
    
    Returns: Figure with -hopefully- enough diagnostics to determine the strength / weaknessess of the detrending
    
    @param flux: flux measurements [adu]
    @type flux: numpy array of length N
    @param position: CCD position measurements along a position axis [pixel]
    @type position: numpy array of length N
    @param AICtck: spline tck for the favored fit by the AIC
    @type AICtck: tuple
    @param BICtck: spline tck for the favored fit by the BIC
    @type BICtck: tuple
    """
    
    # Calculating the corrections
    fluxCORRECTIONaic = scInterp.splev(position, AICtck)
    fluxCORRECTIONbic = scInterp.splev(position, BICtck)
    
    # Setting up the figure window
    # ----------------------------
    figDIAGN = pl.figure(figsize=(16,16))
    axAIC = figDIAGN.add_subplot(221)
    axBIC = figDIAGN.add_subplot(222, sharey=axAIC, sharex=axAIC)
    axAICres = figDIAGN.add_subplot(223, sharey=axAIC, sharex=axAIC)
    axBICres = figDIAGN.add_subplot(224, sharey=axAIC, sharex=axAIC)
    # -- 
    axAIC.plot(position, flux, 'k.', ms=6, alpha=.4)
    axAIC.plot(position, fluxCORRECTIONaic, 'r.', ms=8, alpha=.6)
    pl.tick_params('both',length=10,width=2,which='major'); pl.tick_params('both',length=10,width=1,which='minor')
    pyplot.locator_params(axis = 'x', nbins = 5); pyplot.locator_params(axis = 'y', nbins = 5); axAIC.set_title('best AIC fit'); axAIC.set_ylabel('Flux [adu]')
    # -- 
    axAICres.plot(position, flux - fluxCORRECTIONaic, 'k.', ms=6, alpha=.4)
    pl.tick_params('both',length=10,width=2,which='major'); pl.tick_params('both',length=10,width=1,which='minor')
    pyplot.locator_params(axis = 'x', nbins = 5); pyplot.locator_params(axis = 'y', nbins = 5); axAICres.set_ylabel('Res. Flux [adu]')
    # -- 
    axBIC.plot(position, flux, 'k.', ms=6, alpha=.4)
    axBIC.plot(position, fluxCORRECTIONbic, 'r.', ms=8, alpha=.6)
    pl.tick_params('both',length=10,width=2,which='major'); pl.tick_params('both',length=10,width=1,which='minor')
    pyplot.locator_params(axis = 'x', nbins = 5); pyplot.locator_params(axis = 'y', nbins = 5); axBIC.set_title('best BIC fit')
    # -- 
    axBICres.plot(position, flux - fluxCORRECTIONbic, 'k.', ms=6, alpha=.4)
    pl.tick_params('both',length=10,width=2,which='major'); pl.tick_params('both',length=10,width=1,which='minor')
    pyplot.locator_params(axis = 'x', nbins = 5); pyplot.locator_params(axis = 'y', nbins = 5)
    
    # Settings
    # --------
    axAIC.set_xlim([np.min(position)-0.1, np.max(position)+0.1])
    axAICres.set_xlabel('CCD position [pixel]'); axBICres.set_xlabel('CCD position [pixel]')  
      
    return

def PLOTdetrendPOSITIONfluxFULL(time, flux, position, TCK, **kwargs):
    """
    Routine going with detrendPOSITIONflux to visualize the detrending of the photometry for the instrumental correction between position and flux of BRITE photometry. 
    
    Returns: Figure with -hopefully- enough diagnostics to determine the strength / weaknessess of the detrending
    
    @param flux: flux measurements [adu]
    @type flux: numpy array of length N
    @param position: CCD position measurements along a position axis [pixel]
    @type position: numpy array of length N
    @param TCK: spline tck for the favored fit
    @type TCK: tuple
    """
    
    # Calculating the corrections
    # ---------------------------
    fluxCORRECTION = scInterp.splev(position, TCK)

    
    # Setting up the figure window
    # ----------------------------
    figPOScorr = pl.figure(figsize=(16,16))
    gsPOScorr = gridspec.GridSpec(2, 2,height_ratios=[1,1], width_ratios=[3,2])
    axTIMEorig = figPOScorr.add_subplot(gsPOScorr[0,0])						# Initial flux with time
    axTIMEcorr = figPOScorr.add_subplot(gsPOScorr[1,0], sharex=axTIMEorig, sharey=axTIMEorig)						# Corrected flux with time
    axPOSorig = figPOScorr.add_subplot(gsPOScorr[0,1], sharey=axTIMEorig)						# Initial flux with position
    axPOScorr = figPOScorr.add_subplot(gsPOScorr[1,1], sharex=axPOSorig, sharey=axTIMEorig)						# Corrected flux with position
    
    # Panels related to time
    # ----------------------
    axTIMEorig.plot(time, flux, 'k.', ms=6, alpha=.4)
    pl.tick_params('both',length=10,width=2,which='major'); pl.tick_params('both',length=10,width=1,which='minor')
    pyplot.locator_params(axis = 'x', nbins = 5); pyplot.locator_params(axis = 'y', nbins = 5) 
    axTIMEorig.set_title('Original'); axTIMEorig.set_ylabel('Flux [adu]')
    
    axTIMEcorr.plot(time, flux - fluxCORRECTION, 'k.', ms=6, alpha=.4)
    pl.tick_params('both',length=10,width=2,which='major'); pl.tick_params('both',length=10,width=1,which='minor')
    pyplot.locator_params(axis = 'x', nbins = 5); pyplot.locator_params(axis = 'y', nbins = 5) 
    axTIMEcorr.set_title('After correction'); axTIMEcorr.set_ylabel('Flux [adu]'); axTIMEcorr.set_ylabel('Time [d]')
    # Panels related to position
    # --------------------------------------
    axPOSorig.plot(position, flux, 'k.', ms=6, alpha=.4)
    axPOSorig.plot(position, fluxCORRECTION, 'r.', ms=6, alpha=.8)
    pl.tick_params('both',length=10,width=2,which='major'); pl.tick_params('both',length=10,width=1,which='minor')
    pyplot.locator_params(axis = 'x', nbins = 5); pyplot.locator_params(axis = 'y', nbins = 5); axPOSorig.set_title('Correction') 
    
    axPOScorr.plot(position, flux - fluxCORRECTION, 'k.', ms=6, alpha=.4)
    pl.tick_params('both',length=10,width=2,which='major'); pl.tick_params('both',length=10,width=1,which='minor')
    pyplot.locator_params(axis = 'x', nbins = 5); pyplot.locator_params(axis = 'y', nbins = 5); axPOScorr.set_title('Residuals correction'); axPOScorr.set_ylabel('CCD position [pixel]')
    
    # Settings
    # --------
    axTIMEorig.set_xlim([np.min(time), np.max(time)]); axTIMEorig.set_ylim([np.min(flux)*1.2, np.max(flux)*1.2])
    axPOSorig.set_xlim([np.min(position) + 0.1, np.max(position) + 0.1])
    axTIMEcorr.set_xlabel('Time [d]'); axTIMEorig.set_ylabel('Flux [adu]'); axTIMEcorr.set_ylabel('Flux [adu]'); axPOScorr.set_xlabel('Position [pixel]')
      
    return