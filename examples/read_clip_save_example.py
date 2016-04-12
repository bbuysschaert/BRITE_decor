# -*- coding: utf-8 -*-
"""
This is an example routine to :
- test if all the different subroutines are working;
- provide an example on how to use them;
- get some actual work done.

The testing here is as following:
1. Read *one* setup file;
2. Correct the timing of the observations;
3. Perform outlier clipping on the different parameter spaces (no 2D position outliers yet);
4. Plot the outlier rejection during the  process;
5. Save the clipped data.

TODO add more commenting

Last update 21 March 2016


@author: Bram Buysschaert
"""

#===============================================================================
# 				Packages
#===============================================================================
import numpy as np

import scipy.interpolate as scInterp

from matplotlib.ticker import NullFormatter, MultipleLocator, FormatStrFormatter, NullLocator
from matplotlib import patches, rc, rcParams, pyplot
import matplotlib.gridspec as gridspec

import pylab as pl

import statsmodels.api as sm

# BRITE routines
import BRITE_decor.inout.load as loadBRITE
import BRITE_decor.inout.save as saveBRITE

import BRITE_decor.timing.hjdcorrection as hjdcorrectionBRITE
import BRITE_decor.timing.orbit as orbitBRITE

import BRITE_decor.clipping.medianclipping as medianclipBRITE
import BRITE_decor.clipping.percentageclipping as percentageclipBRITE

import BRITE_decor.gui.outlierrejection as outlierguiBRITE
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

def adjust_timing(HJD, JD, exposureTIME, numberSTACKS, **kwargs):
    """
    A simple routine to convert the HJD to the proper mid-exposure HJDs.
    
    @kwargs: STACKtime: time it takes to stack two observations - Default is 14 [s] (appropriate for almost all observations).
    """
    stackTIME = kwargs.get('STACKtime', 14.0) #[s] 
    
    # Step 1: Get a representation for the current heliocentric correction.
    heliocentricCORRECTION = (HJD - JD) * 24. * 3600. #s      
    tckHELLCORR = hjdcorrectionBRITE.determineHELLCORR(JD, heliocentricCORRECTION)
    
    # Step 2: Determine the mid-observation JD times.
    JD__midOBSERVATION = hjdcorrectionBRITE.determineMIDobsTIME(JD, exposureTIME, stackTIME, np.unique(numberSTACKS))
    
    # Step 3: Determine the heliocentric correction for the corrected JD times.
    heliocentricCORRECTION__midOBSERVATION = scInterp.splev(JD__midOBSERVATION, tckHELLCORR)
    
    # Step 4: Apply the heliocentric correction to the corrected JD times.
    HJD__midOBSERVATION = JD__midOBSERVATION + heliocentricCORRECTION__midOBSERVATION/(3600.*24.)
    
    return HJD__midOBSERVATION   
    
def clip_qualityFLAG(qFLAG, **kwargs):
    """
    Clip on the quality flag provided in DR2 data and return the indexes of the 'outliers'.
    
    NOTE: You can also use the formula from Andrzej's cookbook, and calculate the qFLAG yourself. Both calculations actually differ slightly.
    """
    
    IDXoutliers_qFLAG = np.where(qFLAG == 0)[0]
    
    return IDXoutliers_qFLAG
    
def perform_clipping(time, flux, xPOS, yPOS, temperature, qFLAG, exposureTIME, numberSTACKS, **kwargs):
    """
    The actual clipping routine.
    
    WARNING: The order in which you perform the clipping is actually important!
    1. Aperature flag removal for DR2 data (in case DR1 data you will need to do a stronger correcting in later steps).
    2. Position outliers.
    3. Temperature outliers.
    4. Remove full satellite orbits not having a significant amount of measurements.
    5. Flux outliers (caused by smearing and CTI; so they are most often at the lower end of the data).
    6. Remove full satellite orbits for those who have a very different mean value in that passage.
    
    @kwargs: SATELLITEname: name of the satellite that took the data you are trying to detrend - Default is BHr [string] NOTE: should actually be more of an arg instead of kwarg
    """
    # Reading in the kwars
    satelliteNAME = kwargs.get('SATELLITEname', 'BHr')
    if satelliteNAME == 'BHr':
      satelliteORBITperiod = 97.0972 #[min]
    elif satelliteNAME == 'BLb':
      satelliteORBITperiod = 99.6651 #[min]
    elif satelliteNAME == 'BTr':
      satelliteORBITperiod = 98.2428 #[min]
    elif satelliteNAME == 'BAb':
      satelliteORBITperiod = 100.3617 #[min]
    elif satelliteNAME == 'UBr':
      satelliteORBITperiod = 100.3708 #[min]
    
    # Deep copy of the original arrays, in case you want to plot them at the end for comparison reasons.
    timeORIG, fluxORIG, xPOSORIG, yPOSORIG, temperatureORIG, qFLAGORIG, exposureTIMEORIG, numberSTACKSORIG = np.copy(time), np.copy(flux), np.copy(xPOS), np.copy(yPOS), np.copy(temperature), np.copy(qFLAG), np.copy(exposureTIME), np.copy(numberSTACKS)
    
    # Step 1: Remove qFLAG outliers.
    IDXoutliers = clip_qualityFLAG(qFLAG)
    time, flux, xPOS, yPOS, temperature, qFLAG, exposureTIME, numberSTACKS = np.delete(time,IDXoutliers), np.delete(flux,IDXoutliers), np.delete(xPOS,IDXoutliers), np.delete(yPOS,IDXoutliers), np.delete(temperature,IDXoutliers), np.delete(qFLAG,IDXoutliers), np.delete(exposureTIME,IDXoutliers), np.delete(numberSTACKS,IDXoutliers)
    
    # Step 2: Remove position outliers.
    # This is done in a non unique way, we do it for x and for y, and merge the outlier arrays.
    # You should take care of the long-term trend seen for both positions with time.
    xPOSlowess = sm.nonparametric.lowess(xPOS, time, frac=0.3)
    time_l, xPOS_l = xPOSlowess[:,0], xPOSlowess[:,1]
    
    yPOSlowess = sm.nonparametric.lowess(yPOS, time, frac=0.3)
    time_l, yPOS_l = xPOSlowess[:,0], yPOSlowess[:,1]
    
    # NOTE-1 this is only a crude example. Here, I do not give the user the posibility to adjust the percentage of the array considered as outliers (i.e. 10%).
    # NOTE-2 another posibility is to use a scipy.interpolate.splrep (which might give you a tigher fit to the seen relation). However, you will have to do this yourself then.
    
    IDXoutliers_x = percentageclipBRITE.percentageFILTERonRELATION(time, xPOS, 5., 5., LOWESSfrac=0.25)
    IDXoutliers_y = percentageclipBRITE.percentageFILTERonRELATION(time, yPOS, 5., 5., LOWESSfrac=0.25)
    
    figPOSoutl = pl.figure(figsize=(16,16))
    axX = figPOSoutl.add_subplot(211)
    axX.plot(time, xPOS, 'kx', alpha=.4)
    axX.plot(time_l, xPOS_l, 'r.', alpha=.4)
    axX.plot(time[IDXoutliers_y], xPOS[IDXoutliers_y], 'ys', ms=14)
    axX.plot(time[IDXoutliers_x], xPOS[IDXoutliers_x], 'bs', ms=8)
    axX.set_ylabel('xPOS [pixel]')
    
    axY = figPOSoutl.add_subplot(212, sharex=axX)
    axY.plot(time, yPOS, 'kx', alpha=.4)
    axY.plot(time_l, yPOS_l, 'r.', alpha=.4)
    axY.plot(time[IDXoutliers_y], yPOS[IDXoutliers_y], 'ys', ms=14)
    axY.plot(time[IDXoutliers_x], yPOS[IDXoutliers_x], 'bs', ms=8)
    axY.set_xlabel('Time [d]'); axY.set_ylabel('yPOS [pixel]')
    
    pl.show()
    
    # Determine the unique set of outliers, by merging both outlier arrays
    time, flux, xPOS, yPOS, temperature, qFLAG, exposureTIME, numberSTACKS = np.delete(time,IDXoutliers), np.delete(flux,IDXoutliers), np.delete(xPOS,IDXoutliers), np.delete(yPOS,IDXoutliers), np.delete(temperature,IDXoutliers), np.delete(qFLAG,IDXoutliers), np.delete(exposureTIME,IDXoutliers), np.delete(numberSTACKS,IDXoutliers)
    
    
    # Step 3: Remove temperature outliers.
    # This is done using a GUI.
    IDXoutliers = outlierguiBRITE.interactiveOUTLIER(time, temperature)
    
    time, flux, xPOS, yPOS, temperature, qFLAG, exposureTIME, numberSTACKS = np.delete(time,IDXoutliers), np.delete(flux,IDXoutliers), np.delete(xPOS,IDXoutliers), np.delete(yPOS,IDXoutliers), np.delete(temperature,IDXoutliers), np.delete(qFLAG,IDXoutliers), np.delete(exposureTIME,IDXoutliers), np.delete(numberSTACKS,IDXoutliers)
    
    # Step 4. Remove flux outliers.
    # This is not implemented here, since it requires you to subtract any instrumental and physical signal from the lightcurve. As such, you need an iterative prewhitening code and / or a model for the lightcurve.
    
    
    # Step 5. Remove full satellite orbits not having a significant amount of measurements.
    timeBINmean, fluxBINmean, fluxBINstd, binNUMBERofELEMENTS, IDXoutliers = orbitBRITE.analyseNUMBERperORBIT(time, flux, numberSTACKS, satelliteORBITperiod)
    
    pl.figure()
    pl.plot(time, flux, 'kx', alpha=.4)
    pl.plot(time[IDXoutliers], flux[IDXoutliers], 'bs', ms=8)
    pl.xlabel('Time [d]'); pl.ylabel('Flux [adu]')
    pl.show()
    
    time, flux, xPOS, yPOS, temperature, qFLAG, exposureTIME, numberSTACKS = np.delete(time,IDXoutliers), np.delete(flux,IDXoutliers), np.delete(xPOS,IDXoutliers), np.delete(yPOS,IDXoutliers), np.delete(temperature,IDXoutliers), np.delete(qFLAG,IDXoutliers), np.delete(exposureTIME,IDXoutliers), np.delete(numberSTACKS,IDXoutliers)
    
    return time, flux, xPOS, yPOS, temperature, qFLAG, exposureTIME, numberSTACKS




if __name__ == '__main__':
    # This executes a routine. Normally, you will always change parameters here, unless someone defines an amazing routine to provide input from the command line. However, it is not recommended as you have way too many different options to run the different routines.
    
    # Some parameters you can change to specify your data
    fileIN = 'HD37742_6_OrionII-2014_BHr_setup7_APa2s5_R2.dat'
    fileCLIP = 'HD37742_OrionII_BHr_setup7_clipEXAMPLE.dat'
    
    pathEXAMPLE = 'example_data/'
    
    # Read in the setup file
    HJD, fluxRAW, xCCD, yCCD, tempCCD, JD, qFLAG, expTIME, nSTACKS, xRASTER, yRASTER, aperture = loadBRITE.loadSETUP(fileIN, pathEXAMPLE)
    print 'I have read in the file'
    
    # Convert the timing to the mid-exposure time
    HJDcorrected = adjust_timing(HJD, JD, expTIME, nSTACKS)
    print 'I have adjusted the timing'
    
    # Perform clipping
    HJD_clipped, fluxRAW_clipped, xCCD_clipped, yCCD_clipped, tempCCD_clipped, qFLAG_clipped, expTIME_clipped, nSTACKS_clipped = perform_clipping(HJDcorrected, fluxRAW, xCCD, yCCD, tempCCD, qFLAG, expTIME, nSTACKS)
    print 'I have rejected the outliers'
    
    # Save the clipped parameters
    saveBRITE.saveCLIPPED(fileCLIP, pathEXAMPLE, HJD_clipped, fluxRAW_clipped, xCCD_clipped, yCCD_clipped, tempCCD_clipped, expTIME_clipped, nSTACKS_clipped) # We do not save qFLAG anymore, since it should be = 1 (i.e. aperture fully rendered) be definition of the performed clipping.
    print 'I have saved the output'
    