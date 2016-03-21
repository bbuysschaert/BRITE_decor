# -*- coding: utf-8 -*-
"""
General fitting routine, used in the majority of the detrending scripts. Hence, param could be the CCD temperature, CCD position, or the satellite orbital phase.
    
Last update 21 March 2016

@author: Bram Buysschaert
"""

#===============================================================================
# 				Packages
#===============================================================================
import numpy as np

import scipy.interpolate as scInterp

from BRITE.fitting.splinefit import splineFIT, splineGOODNESSofFITandINFORMATIONCRITERION
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
def detrendPARAMflux(param, flux, **kwargs):
    """
    Routine to perform a general detrending of an instrumental effect seen between param and BRITE flux. It is assumed that the correction is fairly straightforward and requires only one fit.
    
    We permit the user to fit multiple different splines, having a range of knotpoint spacings (also being phase shifted), and different spline orders. Information criteria are then used to determine the most optimal representation.
    
    Returns: A matrix with the AIC, a matrix with the BIC, a matrix with the likelihood, and a matrix containing the TCKs of each performed fit.
    
    @param param: param measurements [???]
    @type param: numpy array of length N
    @param flux: flux measurements [adu]
    @type flux: numpy array of length N
    
    @return matrixAIClocal: matrix of the AIC of each fit
    @rtype: numpy 3D array of size KxPXO
    @return matrixBIClocal: matrix of the BIC of each fit
    @rtype: numpy 3D array of size KxPXO
    @return likelihoodMATRIXlocal: matrix of the likelihood of each fit
    @rtype: numpy 3D array of size KxPXO
    @return matrixTCKlocal: matrix of the likelihood of each fit
    @rtype: numpy 3D array of size KxPXO containing the TCK tuple converted to a string
    
    @kwargs: SPLINEknotpointsSPACING: set of spacing for the different knotpoints - Default is np.array([1./3., 0.5, 1.0]) [deg]
    @kwargs: SPLINEphaseSHIFT: value to consider for the phase shift for the same sets of knotpoints - Default is 0.01 [deg]
    @kwargs: SPLINEorder: order for the spline fits; does accept numpy arrays! - Default is np.array([3],dtype='int32')
    @kwargs: SPLINEtckLENGTH: number of bytes allocated for the to-string-converted TCK - Default is 2000 [integer]
    """
    # Reading in the kwargs and performing some minor checks, so we have everything in the correct input format
    SPLINEknotpointsSPACING = kwargs.get('SPLINEknotpointsSPACING', np.array([1./3., 0.5, 1.0])) #[deg]
    SPLINEphaseSHIFT = kwargs.get('SPLINEphaseSHIFT', 0.01) #[deg]
    SPLINEorder = kwargs.get('SPLINEorder', np.array([3],dtype='int32')) #[integer]
    SPLINEtckLENGTH = 'S' + str(np.int(kwargs.get('SPLINEstringLENGTH', 2000))) # [string] 
    # Making sure you gave a numpy array for SPLINEorder with integers. If not, change it
    try:
      len(SPLINEorder)
    except:
      SPLINEorder = np.array([SPLINEorder],dtype='int32')
    maxNUMBERphaseSHIFTS = np.int(np.max(SPLINEknotpointsSPACING)/SPLINEphaseSHIFT)
    
    # Setting up the local matrices, for which we store the output. We multiply everthing with 1.e50 since we want the minimum BIC / AIC, and in case nothing is calculated, we want to avoid it and being able to trace it. -- If 1.e50 is too small for your usage, you are doing something horribly wrong. -- NOTE that the likelihood should be maximised, thus np.zeros
    matrixAIClocal, matrixBIClocal, likelihoodMATRIXlocal, matrixTCKlocal = np.ones((len(SPLINEknotpointsSPACING),maxNUMBERphaseSHIFTS,len(SPLINEorder)))*1.e50, np.ones((len(SPLINEknotpointsSPACING),maxNUMBERphaseSHIFTS,len(SPLINEorder)))*1.e50, np.zeros((len(SPLINEknotpointsSPACING),maxNUMBERphaseSHIFTS,len(SPLINEorder))), np.ones((len(SPLINEknotpointsSPACING),maxNUMBERphaseSHIFTS,len(SPLINEorder)), dtype=SPLINEtckLENGTH)
    
    # Subtracted the mean of your flux, ensuring to not change any offsets to your flux when calculating the corrections.
    flux = flux - np.mean(flux)
    
    # Sorting the temperature and flux
    paramSORTED, fluxSORTED = (list(t) for t in zip(*sorted(zip(param, flux))))
    paramSORTED, fluxSORTED = np.array(paramSORTED), np.array(fluxSORTED)
    
    # Doing the fitting itself.
    for kk in range(len(SPLINEknotpointsSPACING)):	# Loop over the different sets of knotpoints
      NUMBERphaseSHIFTS = np.int(SPLINEknotpointsSPACING[kk]/SPLINEphaseSHIFT)
      for pp in range(NUMBERphaseSHIFTS): 		# Loop over the different phaseshifts for a given set of knotpoints
	if pp !=0: # A shift
	  paramKNOTPOINTS = np.arange(np.min(param) + SPLINEknotpointsSPACING[kk] - pp*SPLINEphaseSHIFT, np.max(param), SPLINEknotpointsSPACING[kk])
	else: # No shift shift
	  paramKNOTPOINTS = np.arange(np.min(param) + SPLINEknotpointsSPACING[kk], np.max(param), SPLINEknotpointsSPACING[kk])
	# Check whether the last element is too far away or over the maximum value of that position. This happens roughly in ~0.1% of the calculations
	if np.max(param) - paramKNOTPOINTS[-1] > SPLINEknotpointsSPACING[kk]:
	  paramKNOTPOINTS = np.append(paramKNOTPOINTS, paramKNOTPOINTS[-1] + SPLINEknotpointsSPACING[kk])
	elif np.max(param) - paramKNOTPOINTS[-1] <= 0:
	  paramKNOTPOINTS = np.delete(paramKNOTPOINTS, -1)
	
	for oo in range(len(SPLINEorder)): 		# Loop over the different orders of spline
	  NUMBERestimatedPARAMSparam = (len(paramKNOTPOINTS) + 1) * (SPLINEorder[oo] + 1)		# Preferred usage #NOTE +1 here, since the knotpoints DO NOT include beginning and ending
	  ###NUMBERestimatedPARAMSparam = np.int(1./SPLINEknotpointsSPACING[kk]) * (SPLINEorder[oo] + 1)	# Alternative, more prone to overfitting minor effects
	  TCKparam, TCKerror = splineFIT(paramSORTED, fluxSORTED, SPLINEgiveKNOTPOINTS = True, SPLINEknotpoints = paramKNOTPOINTS, SPLINEorder=SPLINEorder[oo])    
	  AICparam, BICparam, likelihoodPARAM = splineGOODNESSofFITandINFORMATIONCRITERION(param, flux, TCKparam, PARAMSdetermine=False, PARAMSestimated=NUMBERestimatedPARAMSparam) # You should provide the full coord array, no rebinned arrays. We provide the number of estimated parameters.
	  matrixAIClocal[kk,pp,oo], matrixBIClocal[kk,pp,oo], likelihoodMATRIXlocal[kk,pp,oo], matrixTCKlocal[kk,pp,oo] = AICparam, BICparam, likelihoodPARAM, str(TCKparam)
    return matrixAIClocal, matrixBIClocal, likelihoodMATRIXlocal, matrixTCKlocal