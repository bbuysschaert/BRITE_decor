# -*- coding: utf-8 -*-
"""
Routines to visualise the detrending for the temperature dependent PSF changes, which produce instrumental effects between flux and position.
    
Last update 21 March 2016

@author: Bram Buysschaert
"""

#===============================================================================
# 				Packages
#===============================================================================
# import CheckMatplotlib # LOCAL ROUTINE, remove
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
def PLOTdetrendTEMPpsfDIAGinformCRIT(flux, AICpos, AICtck, BICpos, BICtck, posORDER, **kwargs):
    """
    Routine going with detrendTEMPpsfFULL to visualize diagnostics in case the AIC and BIC do not favor the same fit.
    
    Returns: Figure with -hopefully- enough diagnostics to determine the strength / weaknessess of the detrending
    
    @param flux: flux measurements [adu]
    @type flux: numpy array of length N
    @param AICpos: CCD position measurements for the favored fit by the AIC [pixel]
    @type AICpos: numpy array of length N
    @param AICtck: spline tck for the favored fit by the AIC
    @type AICtck: tuple
    @param BICpos: CCD position measurements for the favored fit by the BIC [pixel]
    @type BICpos: numpy array of length N
    @param BICtck: spline tck for the favored fit by the BIC
    @type BICtck: tuple
    @param posORDER: indicating which axis is AICpos / BICpos
    @type posORDER: string
    """
    if not posORDER in ['xx', 'xy', 'yx', 'yy']:
      print bcolors.FAIL +  '\t\ERROR wrong posORDER parameter.\n\t\tOnly accepting "xy" or "yx"...\n\t\tDid not plot anything and now exiting...' + bcolors.ENDC
      return
    
    # Calculating the corrections
    fluxCORRECTIONaic = scInterp.splev(AICpos, AICtck)
    fluxCORRECTIONbic = scInterp.splev(BICpos, BICtck)
    
    # Setting up the figure window
    # ----------------------------
    figDIAGN = pl.figure(figsize=(16,16))
    axAIC = figDIAGN.add_subplot(221)
    axBIC = figDIAGN.add_subplot(222, sharey=axAIC)
    axAICres = figDIAGN.add_subplot(223, sharey=axAIC, sharex=axAIC)
    axBICres = figDIAGN.add_subplot(224, sharey=axAIC, sharex=axBIC)
    # -- 
    axAIC.plot(AICpos, flux, 'k.', ms=6, alpha=.4)
    axAIC.plot(AICpos, fluxCORRECTIONaic, 'r.', ms=8, alpha=.6)
    pl.tick_params('both',length=10,width=2,which='major'); pl.tick_params('both',length=10,width=1,which='minor')
    pyplot.locator_params(axis = 'x', nbins = 5); pyplot.locator_params(axis = 'y', nbins = 5); axAIC.set_title('best AIC fit'); axAIC.set_ylabel('Flux [adu]')
    # -- 
    axAICres.plot(AICpos, flux - fluxCORRECTIONaic, 'k.', ms=6, alpha=.4)
    pl.tick_params('both',length=10,width=2,which='major'); pl.tick_params('both',length=10,width=1,which='minor')
    pyplot.locator_params(axis = 'x', nbins = 5); pyplot.locator_params(axis = 'y', nbins = 5); axAICres.set_ylabel('Res. Flux [adu]')
    # -- 
    axBIC.plot(BICpos, flux, 'k.', ms=6, alpha=.4)
    axBIC.plot(BICpos, fluxCORRECTIONbic, 'r.', ms=8, alpha=.6)
    pl.tick_params('both',length=10,width=2,which='major'); pl.tick_params('both',length=10,width=1,which='minor')
    pyplot.locator_params(axis = 'x', nbins = 5); pyplot.locator_params(axis = 'y', nbins = 5); axBIC.set_title('best BIC fit')
    # -- 
    axBICres.plot(BICpos, flux - fluxCORRECTIONbic, 'k.', ms=6, alpha=.4)
    pl.tick_params('both',length=10,width=2,which='major'); pl.tick_params('both',length=10,width=1,which='minor')
    pyplot.locator_params(axis = 'x', nbins = 5); pyplot.locator_params(axis = 'y', nbins = 5)
    
    # Settings
    # --------
    axAIC.set_xlim([np.min(AICpos)-0.1, np.max(AICpos)+0.1]); axBIC.set_xlim([np.min(BICpos)-0.1, np.max(BICpos)+0.1])
    
    if posORDER == 'xx':
      axAICres.set_xlabel('xPOS [pixel]'); axBICres.set_xlabel('xPOS [pixel]')
    elif posORDER == 'xy':
      axAICres.set_xlabel('xPOS [pixel]'); axBICres.set_xlabel('yPOS [pixel]')  
    elif posORDER == 'yx':
      axAICres.set_xlabel('yPOS [pixel]'); axBICres.set_xlabel('xPOS [pixel]') 
    elif posORDER == 'yy':
      axAICres.set_xlabel('yPOS [pixel]'); axBICres.set_xlabel('yPOS [pixel]')   
      
    return

def PLOTdetrendTEMPpsfFULL(time, flux, POS1, tck1, POS2, tck2, posORDER, **kwargs):
    """
    Routine going with detrendTEMPpsfFULL to visualize the detrending of the photometry for the temperature dependent PSF of BRITE photometry. 
    
    Returns: Figure with -hopefully- enough diagnostics to determine the strength / weaknessess of the detrending
    
    @param time: time measurements [d]
    @type time: numpy array of length N
    @param flux: flux measurements [adu]
    @type flux: numpy array of length N
    @param POS1: CCD position measurements along the first position axis [pixel]
    @type POS1: numpy array of length N
    @param tck1: spline tck for the first detrending
    @type tck1: tuple
    @param POS2: CCD position measurements along the second position axis [pixel]
    @type POS2: numpy array of length N
    @param tck1: spline tck for the second detrending
    @type tck2: tuple
    @param posORDER: indicating which axis is POS1 / POS2
    @type posORDER: string
    """
    if not posORDER in ['xy', 'yx']:
      print bcolors.FAIL +  '\t\ERROR wrong posORDER parameter.\n\t\tOnly accepting "xy" or "yx"...\n\t\tDid not plot anything and now exiting...' + bcolors.ENDC
      return
    
    # Calculating the corrections
    # ---------------------------
    fluxCORRECTION1 = scInterp.splev(POS1, tck1)
    fluxCORRECTION2 = scInterp.splev(POS2, tck2)

    
    # Setting up the figure window
    # ----------------------------
    figPSFcorr = pl.figure(figsize=(20,16))
    gsPSFcorr = gridspec.GridSpec(3, 3,height_ratios=[1,1,1], width_ratios=[3,3,2])
    axTIME0 = figPSFcorr.add_subplot(gsPSFcorr[0,0])						# Flux with time initially
    axTIME1 = figPSFcorr.add_subplot(gsPSFcorr[1,0], sharex=axTIME0, sharey=axTIME0)		# Flux with time after first correction
    axTIME2 = figPSFcorr.add_subplot(gsPSFcorr[2,0], sharex=axTIME0, sharey=axTIME0)		# Flux with time after second correction
    # --
    axCORR1 = figPSFcorr.add_subplot(gsPSFcorr[0,1], sharey=axTIME0)				# First correction
    axCORR2 = figPSFcorr.add_subplot(gsPSFcorr[1,1], sharey=axTIME0)				# Second correction
    # --
    axCORR1res = figPSFcorr.add_subplot(gsPSFcorr[0,2], sharey=axTIME0, sharex=axCORR1)		# Residuals after first correction
    axCORR2res = figPSFcorr.add_subplot(gsPSFcorr[1,2], sharey=axTIME0, sharex=axCORR2)		# Residuals after second correction
    # --
    axFINALres1 = figPSFcorr.add_subplot(3,4,11, sharey=axTIME0, sharex=axCORR1)	# Final residuals with POS1
    axFINALres2 = figPSFcorr.add_subplot(3,4,12, sharey=axTIME0, sharex=axCORR2)	# Final residuals with POS2
    
    # Panels related to time
    # ----------------------
    axTIME0.plot(time, flux, 'k.', ms=6, alpha=.4)
    pl.tick_params('both',length=10,width=2,which='major'); pl.tick_params('both',length=10,width=1,which='minor')
    pyplot.locator_params(axis = 'x', nbins = 5); pyplot.locator_params(axis = 'y', nbins = 5) 
    axTIME0.set_title('Original')
    
    axTIME1.plot(time, flux - fluxCORRECTION1, 'k.', ms=6, alpha=.4)
    pl.tick_params('both',length=10,width=2,which='major'); pl.tick_params('both',length=10,width=1,which='minor')
    pyplot.locator_params(axis = 'x', nbins = 5); pyplot.locator_params(axis = 'y', nbins = 5) 
    axTIME1.set_title('After 1st correction'); axTIME1.set_ylabel('Flux [adu]')
    
    axTIME2.plot(time, flux - fluxCORRECTION1 - fluxCORRECTION2, 'k.', ms=6, alpha=.4)
    pl.tick_params('both',length=10,width=2,which='major'); pl.tick_params('both',length=10,width=1,which='minor')
    pyplot.locator_params(axis = 'x', nbins = 5); pyplot.locator_params(axis = 'y', nbins = 5)  
    axTIME2.set_title('After 2nd correction'); axTIME2.set_xlabel('Time [d]')
    
    # Panels related to the first correction
    # --------------------------------------
    axCORR1.plot(POS1, flux, 'k.', ms=6, alpha=.4)
    axCORR1.plot(POS1, fluxCORRECTION1, 'r.', ms=6, alpha=.8)
    pl.tick_params('both',length=10,width=2,which='major'); pl.tick_params('both',length=10,width=1,which='minor')
    pyplot.locator_params(axis = 'x', nbins = 5); pyplot.locator_params(axis = 'y', nbins = 5); axCORR1.set_title('1st correction') 
    
    axCORR1res.plot(POS1, flux - fluxCORRECTION1, 'k.', ms=6, alpha=.4)
    pl.tick_params('both',length=10,width=2,which='major'); pl.tick_params('both',length=10,width=1,which='minor')
    pyplot.locator_params(axis = 'x', nbins = 5); pyplot.locator_params(axis = 'y', nbins = 5); axCORR1res.set_title('Residuals 1st correction') 
    
    # Panels related to the first correction
    # --------------------------------------    
    axCORR2.plot(POS2, flux - fluxCORRECTION1, 'k.', ms=6, alpha=.4)
    axCORR2.plot(POS2, fluxCORRECTION2, 'r.', ms=6, alpha=.8)
    pl.tick_params('both',length=10,width=2,which='major'); pl.tick_params('both',length=10,width=1,which='minor')
    pyplot.locator_params(axis = 'x', nbins = 5); pyplot.locator_params(axis = 'y', nbins = 5); axCORR2.set_title('2nd correction') 
    
    axCORR2res.plot(POS2, flux - fluxCORRECTION1 - fluxCORRECTION2, 'k.', ms=6, alpha=.4)
    pl.tick_params('both',length=10,width=2,which='major'); pl.tick_params('both',length=10,width=1,which='minor')
    pyplot.locator_params(axis = 'x', nbins = 5); pyplot.locator_params(axis = 'y', nbins = 5); axCORR2res.set_title('Residuals 2nd correction') 
    
    # Panels related to the final residuals
    # -------------------------------------    
    axFINALres1.plot(POS1, flux - fluxCORRECTION1 - fluxCORRECTION2, 'k.', ms=6, alpha=.4)
    pl.tick_params('both',length=10,width=2,which='major'); pl.tick_params('both',length=10,width=1,which='minor')
    pyplot.locator_params(axis = 'x', nbins = 5); pyplot.locator_params(axis = 'y', nbins = 5); axFINALres1.set_title('Final residuals') 
    
    axFINALres2.plot(POS2, flux - fluxCORRECTION1 - fluxCORRECTION2, 'k.', ms=6, alpha=.4)
    pl.tick_params('both',length=10,width=2,which='major'); pl.tick_params('both',length=10,width=1,which='minor')
    pyplot.locator_params(axis = 'x', nbins = 5); pyplot.locator_params(axis = 'y', nbins = 5); axFINALres2.set_title('Final residuals')  
    
    
    
    # Settings
    # --------
    axTIME0.set_xlim([np.min(time), np.max(time)]); axTIME0.set_ylim([np.min(flux)*1.2, np.max(flux)*1.2])
    axCORR1.set_xlim([np.min(POS1) + 0.1, np.max(POS1) + 0.1]); axCORR2.set_xlim([np.min(POS2) + 0.1, np.max(POS2) + 0.1]); 
    if posORDER == 'xy':
      axCORR1.set_xlabel('xPOS [pixel]'); axCORR1res.set_xlabel('xPOS [pixel]'); axFINALres1.set_xlabel('xPOS [pixel]')
      axCORR2.set_xlabel('yPOS [pixel]'); axCORR2res.set_xlabel('yPOS [pixel]'); axFINALres2.set_xlabel('yPOS [pixel]')
    elif posORDER == 'yx':
      axCORR1.set_xlabel('yPOS [pixel]'); axCORR1res.set_xlabel('yPOS [pixel]'); axFINALres1.set_xlabel('yPOS [pixel]')
      axCORR2.set_xlabel('xPOS [pixel]'); axCORR2res.set_xlabel('xPOS [pixel]'); axFINALres2.set_xlabel('xPOS [pixel]')      
      
    return
