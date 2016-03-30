# -*- coding: utf-8 -*-
"""
This is an example routine to :
- show you how some of the detrending could be done;
- get some actual work done.

The testing here is as following:
1. Read the clipped data;
2. Read the PSF corrected data
3. Perform the detrending with CCD temperature;
4. Perform the detrending with CCD x-position;
5. Perform the detrending with CCD y-position;
6. Save the corrected data.

TODO add more commenting
TODO include the correlation matrix to decide the order of the corrections and see if a correction is actually necessary

Last update 21 March 2016


@author: Bram Buysschaert
"""

#===============================================================================
# 				Packages
#===============================================================================
#import CheckMatplotlib # LOCAL ROUTINE, remove
import numpy as np

import scipy.interpolate as scInterp

from matplotlib.ticker import NullFormatter, MultipleLocator, FormatStrFormatter, NullLocator
from matplotlib import patches, rc, rcParams, pyplot
import matplotlib.gridspec as gridspec

import pylab as pl

import statsmodels.api as sm

# BRITE routines
import BRITE.inout.load as loadBRITE
import BRITE.inout.save as saveBRITE

import BRITE.detrending.detrendPositionFlux as POSdetrendBRITE
import BRITE.detrending.detrendTempFlux as TEMPdetrendBRITE
#===============================================================================
# 				Settings
#===============================================================================
from string import whitespace as ws
import re
rc('xtick', labelsize=20) 
rc('ytick', labelsize=20) 
rc("lines", markeredgewidth=2.0)
rc("axes", linewidth=2.0)
#rc('text', usetex=True)
rcParams["font.size"] = 25

p = re.compile('(%s)' % ('|'.join([c for c in ws])))
s = " \ nu "
nu = p.sub('',s)
#===============================================================================
# 				Code
#===============================================================================

if __name__ == '__main__':
    # This executes a routine. Normally, you will always change parameters here, unless someone defines an amazing routine to provide input from the command line. However, it is not recommended as you have way too many different options to run the different routines.
    
    # WARNING you should do a proper position clipping before you try to run these routines, for an optimal result
    
    # Some parameters you can change to specify your data
    #fileCLIP = 'HD37742_OrionII_BHr_setup7_clipEXAMPLE.dat'
    fileCLIP = 'HD37742_BHr_CLIPPED_OrionII.dat'
    filePSFcorr_flux = 'HD37742_BHr_OrionII_decorr4PSF_FLUX.dat'
    fileCORR = 'HD37742_BHr_OrionII_corrected.dat'
    
    pathEXAMPLE = '/home/bram/python/BRITE/exampleBRITEroutines'
    
    # Read in the clipped file
    HJD, fluxRAW, xCCD, yCCD, tempCCD, expTIME, nSTACKS = loadBRITE.loadCLIPPED(fileCLIP, pathEXAMPLE)
    print 'I have read in a file'
    
    # Read in the PSF corrected flux
    HJD, fluxRAW, fluxCORRECTED, PSFcorrection = loadBRITE.loadPSFdetrend(filePSFcorr_flux, pathEXAMPLE)
    print 'I have read in a file'
    
    # Performing the detrending with temperature
    tckTEMPcorrection = TEMPdetrendBRITE.detrendTEMPflux(HJD, tempCCD, fluxCORRECTED, show_ME=True)
    TEMPcorrection = scInterp.splev(tempCCD, tckTEMPcorrection)
    fluxCORRECTED = fluxCORRECTED - TEMPcorrection
    print 'I have corrected the flux for temperature changes'
    
    # Performing the detrending with xCCD
    tckxCCDcorrection = POSdetrendBRITE.detrendPOSITIONflux(HJD, xCCD, fluxCORRECTED, show_ME=True)
    xCCDcorrection = scInterp.splev(xCCD, tckxCCDcorrection)
    fluxCORRECTED = fluxCORRECTED - xCCDcorrection
    print 'I have corrected the flux for xCCD changes'
    
    # Performing the detrending with yCCD
    tckyCCDcorrection = POSdetrendBRITE.detrendPOSITIONflux(HJD, yCCD, fluxCORRECTED, show_ME=True)
    yCCDcorrection = scInterp.splev(yCCD, tckyCCDcorrection)
    fluxCORRECTED = fluxCORRECTED - yCCDcorrection
    print 'I have corrected the flux for yCCD changes'
    
    
    
    figCORRECTION = pl.figure(figsize=(16,16))
    axTEMPcorr = figCORRECTION.add_subplot(311)
    axxCCDcorr = figCORRECTION.add_subplot(312, sharex=axTEMPcorr, sharey=axTEMPcorr)
    axyCCDcorr = figCORRECTION.add_subplot(313, sharex=axTEMPcorr, sharey=axTEMPcorr)
    
    axTEMPcorr.plot(HJD, TEMPcorrection, 'k.', ms=6, alpha=.4)
    axxCCDcorr.plot(HJD, xCCDcorrection, 'k.', ms=6, alpha=.4)
    axyCCDcorr.plot(HJD, yCCDcorrection, 'k.', ms=6, alpha=.4)
    
    
    figCORRECTION2 = pl.figure(figsize=(16,16))
    axRAW = figCORRECTION2.add_subplot(221)
    axFINAL = figCORRECTION2.add_subplot(223, sharex=axRAW, sharey=axRAW)
    axRAWfreq = figCORRECTION2.add_subplot(222)
    axFINALfreq = figCORRECTION2.add_subplot(224, sharex=axRAWfreq, sharey=axRAWfreq)
    
    axRAW.plot(HJD, fluxRAW, 'k.', ms=6, alpha=.4)
    axFINAL.plot(HJD, fluxCORRECTED, 'k.', ms=6, alpha=.4)
    
    # OWN routines, which are private
    import ivs.timeseries.pergrams as pg
    print 'I am starting the frequency calculations'
    freqORIG, amplORIG = pg.scargle(HJD, fluxRAW - np.mean(fluxRAW), df=0.001, fn=25)
    print '1/2'
    freqCORR, amplCORR = pg.scargle(HJD, fluxCORRECTED - np.mean(fluxCORRECTED), df=0.001, fn=25)
    
    axRAWfreq.plot(freqORIG, amplORIG, 'k-', lw=3, alpha=.4)
    axFINALfreq.plot(freqCORR, amplCORR, 'k-', lw=3, alpha=.4)
    
    pl.show()
    
    
    
