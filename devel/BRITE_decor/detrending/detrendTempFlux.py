# -*- coding: utf-8 -*-
"""
Routines to detrend for the temperature - flux instrumental effects
    
Last update 21 March 2016

@author: Bram Buysschaert
"""

#===============================================================================
# 				Packages
#===============================================================================
import numpy as np
import pylab as pl

import scipy.interpolate as scInterp

from BRITE_decor.fitting.splinefit import splineFIT, splineGOODNESSofFITandINFORMATIONCRITERION, reconvertTCKfromSTRING
from BRITE_decor.detrending.detrendPARAMFlux import detrendPARAMflux
from BRITE_decor.plotting.PLOTdetrendTempFlux import PLOTdetrendTEMPfluxFULL, PLOTdetrendTEMPfluxDIAGinformCRIT
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
def detrendTEMPflux(time, temperature, flux, **kwargs):
    """
    Routine to perform the detrending of the instrumental effect seen between the CCD temperature and BRITE flux.
    
    We call the general routine BRITE.detrending.detrendPARAMflux and use temperature as the param. Next, we update the kwarg dictionary, so appropriate settings are achieved, even when the user did not specify anything for the kwargs.
    
    NOTE see the following link on dictionary manipulation
    http://www.tutorialspoint.com/python/python_dictionary.htm
    
    Returns: The TCKoptimal for the most optimal fit (see also scipy.interpolate.splrep)
    
    @param time: time measurements [d]
    @type time: numpy array of length N
    @param temperature: temperature measurements [deg]
    @type temperature: numpy array of length N
    @param flux: flux measurements [adu]
    @type flux: numpy array of length N
    
    @return TCKoptimal: TCK of the most optimal spline representation
    @rtype: tuple
    
    @kwargs: SPLINEknotpointsSPACING: set of spacing for the different knotpoints - Default is np.array([1./3., 0.5, 1.0, 2.0]) [deg]
    @kwargs: SPLINEphaseSHIFT: value to consider for the phase shift for the same sets of knotpoints - Default is 0.01 [deg]
    @kwargs: SPLINEorder: order for the spline fits; does accept numpy arrays! - Default is np.array([3],dtype='int32')
    @kwargs: SPLINEtckLENGTH: number of bytes allocated for the to-string-converted TCK - Default is 2000 [integer]
    
    @kwargs: show_ME: Boolean to indicate if you want plotting at each possible step - Default is False [Boolean]
    @kwargs: show_FITS: Boolean to indicate if you want plotting after each bin fitting step - Default is False [Boolean]
    @kwargs: show_DIAG: Boolean to indicate if you want plotting when AIC != BIC - Default is False [Boolean]
    """
        
    # Reading in the kwargs and performing some minor checks, so we have everything in the correct input format
    SPLINEknotpointsSPACING = kwargs.get('SPLINEknotpointsSPACING', np.array([1./3., 0.5, 1.0, 2.0])) #[deg]
    SPLINEphaseSHIFT = kwargs.get('SPLINEphaseSHIFT', 0.01) #[deg]
    SPLINEorder = kwargs.get('SPLINEorder', np.array([3],dtype='int32')) #[integer]
    SPLINEtckLENGTH = np.int(kwargs.get('SPLINEstringLENGTH', 2000)) # [integer] (MINOR change needed for updating the directionary)
    show_ME = kwargs.get('show_ME', False); show_FITS = kwargs.get('show_FITS', False); show_DIAG = kwargs.get('show_DIAG', False)
    
    # Manipulating the kwarg dictonary, so we define the kwargs here for detrendPARAMflux in case the user did not specify anything. (REMEMBER: general routines require general fitting conditions, which are far from appropriate for *this* application, hence the mandatory update)
    kwargs['SPLINEknotpointsSPACING'] = SPLINEknotpointsSPACING; kwargs['SPLINEphaseSHIFT'] = SPLINEphaseSHIFT; kwargs['SPLINEorder'] = SPLINEorder; kwargs['SPLINEstringLENGTH'] = SPLINEtckLENGTH
    
    # Subtracted the mean of your flux, ensuring to not change any offsets to your flux when calculating the corrections.
    flux = flux - np.mean(flux)
    
    # Doing the fitting
    matrixAIC, matrixBIC, likelihoodMATRIX, matrixTCK = detrendPARAMflux(temperature, flux, **kwargs)
    
    # From here on, generic sequence of checking the fitting
    # ------------------------------------------------------
    
    # Looking for the optimal solution to apply the correction. To do so, we use the AIC and BIC. In doubt, we also resort to the likelihood. -- The nanmin is important, since there might be some NaNs in the matrix --
    AICmin = np.where(matrixAIC==np.nanmin(matrixAIC))
    BICmin = np.where(matrixBIC==np.nanmin(matrixBIC))
    
    # Determine the number of NaNs (sum over a matrix with ones + the np.isnan function) -- DEBUGGING at the moment
    NUMBERofNANAIC, NUMBERofNANBIC = np.sum(np.ones_like(matrixAIC[:,:,:])[np.isnan(matrixAIC[:,:,:])]), np.sum(np.ones_like(matrixBIC[:,:,:])[np.isnan(matrixBIC[:,:,:])])
    
    if NUMBERofNANAIC != 0 or  NUMBERofNANBIC != 0: #Not really used at present
      matrixAIC[np.isnan(matrixAIC)], matrixBIC[np.isnan(matrixBIC)] = 1.e75, 1.e75
      
      # Counting all values that are 1.e50, comparing them to the number of NaNs and the total number of entries
      matrixAICUNIQUE, matrixAICUNIQUEcounts = np.unique(matrixAIC, return_counts=True)
      matrixBICUNIQUE, matrixBICUNIQUEcounts = np.unique(matrixBIC, return_counts=True)
      
      NUMBERunusedAIC, NUMBERunusedBIC = matrixAICUNIQUEcounts[np.where(matrixAICUNIQUE == 1.e50)][0], matrixBICUNIQUEcounts[np.where(matrixBICUNIQUE == 1.e50)][0]
      
    # Ensure that both the AIC and BIC point towards the same optimal solution.
    if AICmin != BICmin:
      print bcolors.WARNING + '\t\tWARNING: AIC and BIC point towards a different set of knotpoints\n\t\tProviding more information now...' + bcolors.ENDC
      print bcolors.WARNING + '\t\tWARNING: AIC states {:1.3f} and BIC states {:1.3f} for the spacing of the knotpoints'.format(SPLINEknotpointsSPACING[AICmin[0][0]], SPLINEknotpointsSPACING[BICmin[0][0]]) + bcolors.ENDC
      print bcolors.WARNING + '\t\tWARNING: AIC states {:1.3f} and BIC states {:1.3f} for the phase shift of those sets of knotpoints'.format(AICmin[1][0]*SPLINEphaseSHIFT, BICmin[1][0]*SPLINEphaseSHIFT) + bcolors.ENDC
      print bcolors.WARNING + '\t\tWARNING: AIC states {:1.3f} and BIC states {:1.3f} for the order of the spline '.format(SPLINEorder[AICmin[2][0]], SPLINEorder[BICmin[2][0]]) + bcolors.ENDC
      print bcolors.WARNING + '\t\tWARNING: AIC states a minimum of {:.5e}, compared to {:.5e} for the BIC optimal solution'.format(matrixAIC[AICmin][0], matrixAIC[BICmin][0]) + bcolors.ENDC
      print bcolors.WARNING + '\t\tWARNING: BIC states a minimum of {:.5e}, compared to {:.5e} for the AIC optimal solution'.format(matrixBIC[BICmin][0], matrixBIC[AICmin][0]) + bcolors.ENDC
      print bcolors.WARNING + '\t\tWARNING: Likelihoods are {:.5e} and {:.5e}, for AIC and BIC respectively'.format(likelihoodMATRIX[AICmin][0], likelihoodMATRIX[BICmin][0]) + bcolors.ENDC
      
      # We follow the suggestion of the likelihood in case the AIC and BIC do not agree for the optimal solution.
      likelihoodAIC, likelihoodBIC = likelihoodMATRIX[AICmin][0], likelihoodMATRIX[BICmin][0] # The [0] is needed because of the np.where() function.
      if likelihoodAIC > likelihoodBIC:
        print bcolors.FAIL + '\t\tTaking the AIC as best solution, indicated by the LIKELIHOOD' + bcolors.ENDC
        print bcolors.FAIL + '\t\tTaking the AIC as best solution, indicated by the LIKELIHOOD' + bcolors.ENDC
        correctionPARAMS = AICmin
      elif likelihoodAIC < likelihoodBIC:
        print bcolors.FAIL + '\t\tTaking the **BIC** as best solution, indicated by the LIKELIHOOD' + bcolors.ENDC
        print bcolors.FAIL + '\t\tTaking the **BIC** as best solution, indicated by the LIKELIHOOD' + bcolors.ENDC
        correctionPARAMS = BICmin
        
      # Provide some diagnostics here for debugging, so we can actually study if the likelihood is a good suggestion.
      if show_ME or show_DIAG:
	TCKaic = reconvertTCKfromSTRING(matrixTCK[AICmin][0])
	TCKbic = reconvertTCKfromSTRING(matrixTCK[BICmin][0])
	PLOTdetrendTEMPfluxDIAGinformCRIT(flux, temperature, TCKaic, TCKbic)
	pl.show()
      
    elif AICmin == BICmin:
      correctionPARAMS = AICmin # Should not matter which one you take, since AICmin == BICmin
    
    # Perform the correction using the optimal settings for the fit   
    TCKoptimal = reconvertTCKfromSTRING(matrixTCK[correctionPARAMS][0])
    
    if show_ME or show_FITS:
      PLOTdetrendTEMPfluxFULL(time, flux, temperature, TCKoptimal)
      pl.show()
      
    return TCKoptimal  