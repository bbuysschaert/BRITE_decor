# -*- coding: utf-8 -*-
"""
This is an example routine to :
- show you how the time binning works using temperature (and / or gaps);
- show you how the PSF detrending works;
- get some actual work done.

The testing here is as following:
1. Read the clipped data;
2. Perform the binning of the time array using the temperature;
3. Correct the photometry within each bin for the temperature-dependent PSF changes;
4. Save the corrected data.

TODO add more commenting

Last update 21 March 2016


@author: Bram Buysschaert
"""

#===============================================================================
# 				Packages
#===============================================================================
import CheckMatplotlib # LOCAL ROUTINE, remove
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
import BRITE.inout.OUTdetrendTempPSF as savePSFdetrendBRITE

import BRITE.detrending.detrendTempPSF as PSFdetrendBRITE
import BRITE.detrending.binningWtemperature as TEMPbinBRITE

import BRITE.plotting.PLOTbinningWtemperature as plotTEMPbinBRITE
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
    filePSFcorr_diagnostic = 'HD37742_BHr_OrionII_decorr4PSF_BINS.dat'
    filePSFcorr_flux = 'HD37742_BHr_OrionII_decorr4PSF_FLUX.dat'
    
    pathEXAMPLE = '/home/bram/python/BRITE/exampleBRITEroutines'
    
    # Read in the clipped file
    HJD, fluxRAW, xCCD, yCCD, tempCCD, expTIME, nSTACKS = loadBRITE.loadCLIPPED(fileCLIP, pathEXAMPLE)
    print 'I have read in the file'
    
    
    # Do the temperature binning
    binENDindexes, binENDindexes_reason, tempCCDLONG = TEMPbinBRITE.openTIMEwithTEMPERATURE(HJD, tempCCD, TEMPcrit=1., GAPSinclude=True, GAPScrit=1., LOWESSfrac=0.30)
    print 'I have binned the data'
    
    
    # Plot the temperature binning    
    plotTEMPbinBRITE.PLOTopenTIMEwithTEMPERATURE(HJD, fluxRAW, tempCCD, tempCCDLONG, xCCD, yCCD, binENDindexes, binENDindexes_reason)
    pl.show()
    
    
    # Perform the detrending (and showing it, since show_ME=True)
    fluxPSFcorrected, correction, tckFIRSTstring_list, tckSECONDstring_list, diagnostic_list = np.array([]), np.array([]), [], [], []
    for ii in range(len(binENDindexes)):
      if ii == 0: # The first index
        HJD_bin, fluxRAW_bin, xCCD_bin, yCCD_bin = HJD[:binENDindexes[ii]+1], fluxRAW[:binENDindexes[ii]+1], xCCD[:binENDindexes[ii]+1], yCCD[:binENDindexes[ii]+1]
        
      elif ii == len(binENDindexes) -1: # The last index
        HJD_bin, fluxRAW_bin, xCCD_bin, yCCD_bin = HJD[binENDindexes[ii-1]+1:], fluxRAW[binENDindexes[ii-1]+1:], xCCD[binENDindexes[ii-1]+1:], yCCD[binENDindexes[ii-1]+1:]
        
      else: # The other indexes
	HJD_bin, fluxRAW_bin, xCCD_bin, yCCD_bin = HJD[binENDindexes[ii-1]+1:binENDindexes[ii]+1], fluxRAW[binENDindexes[ii-1]+1:binENDindexes[ii]+1], xCCD[binENDindexes[ii-1]+1:binENDindexes[ii]+1], yCCD[binENDindexes[ii-1]+1:binENDindexes[ii]+1]
      
      correction_bin, tckFIRSTstring_bin, tckSECONDstring_bin, diagnostic_bin = PSFdetrendBRITE.detrendTEMPpsfFULL(HJD_bin, fluxRAW_bin, xCCD_bin, yCCD_bin, show_ME=False, SPLINEknotpointsSPACING = np.array([0.1, 0.2,0.25,1./3.]), SPLINEorder = np.array([3,5], dtype='int32'))
      
      fluxPSFcorrected, correction = np.append(fluxPSFcorrected, fluxRAW_bin - correction_bin), np.append(correction, correction_bin)
      tckFIRSTstring_list.append(tckFIRSTstring_bin); tckSECONDstring_list.append(tckSECONDstring_bin);  diagnostic_list.append(diagnostic_bin)
    print 'I have corrected the flux'
   
    # Save the output
    savePSFdetrendBRITE.OUT_DIAGNOSTIC_detrendTEMPpsfFULL(filePSFcorr_diagnostic, pathEXAMPLE, HJD, binENDindexes, binENDindexes_reason, tckFIRSTstring_list, tckSECONDstring_list, diagnostic_list)
    savePSFdetrendBRITE.OUT_FLUX_detrendTEMPpsfFULL(filePSFcorr_flux, pathEXAMPLE, HJD, fluxRAW, fluxPSFcorrected, correction)   
    
    print 'I have saved the correction'
    
    
    
    