# -*- coding: utf-8 -*-
"""
Routines to detrend for the temperature dependent PSF changes, which produce instrumental effects between flux and position.

In general, you should call the routine detrendTEMPpsfFULL.

WARNING The routine will fail in case the positions are not well clipped. This can (and will) happen for even *one* bad datapoint. At present, I do not provide a routine to redo the fitting. Your best bet is to use a try-except statement, where you change the used spacings (make them larger!) for the knotpoints in the except window (or even ignore it at that moment). Remember, you can use the doSILENT = False to see why it happens. But it will only show you an errormessage=10, i.e. error on the input data (and NaNs in the whole AIC/BIC matrix; besides the 1.e50s).
    
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
from BRITE_decor.plotting.PLOTdetrendTempPSF import PLOTdetrendTEMPpsfFULL, PLOTdetrendTEMPpsfDIAGinformCRIT
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
def detrendTEMPpsfCOORD(flux, coord, **kwargs): # actually redundant and should be replaced by BRITE.detrending.detrendPARAMflux
    """
    Routine to perform the fitting between a CCD coordinate vs flux for the temperature dependent PSF of BRITE photometry, using the full position arrays. This routine is three times called during fitPSFpositionFULL. This saves some coding lines and makes everything easier to understand / clearer on what happens.
    
    Returns: A matrix with the AIC, a matrix with the BIC, a matrix with the likelihood, and a matrix containing the TCKs of each performed fit.
    
    @param flux: flux measurements [adu]
    @type flux: numpy array of length N
    @param coord: CCD position measurements along an axis [pixel]
    @type coord: numpy array of length N
    
    
    @return matrixAIClocal: matrix of the AIC of each fit
    @rtype: numpy 3D array of size KxPXO
    @return matrixBIClocal: matrix of the BIC of each fit
    @rtype: numpy 3D array of size KxPXO
    @return likelihoodMATRIXlocal: matrix of the likelihood of each fit
    @rtype: numpy 3D array of size KxPXO
    @return matrixTCKlocal: matrix of the likelihood of each fit
    @rtype: numpy 3D array of size KxPXO containing the TCK tuple converted to a string
    
    @kwargs: SPLINEknotpointsSPACING: set of spacing for the different knotpoints - default is np.array([0.2,0.25,1./3.,0.5]) [pixel]
    @kwargs: SPLINEphaseSHIFT: value to consider for the phase shift for the same sets of knotpoints - default is 0.01 [pixel]
    @kwargs: SPLINEtckLENGTH: number of bytes allocated for the to-string-converted TCK - default is 2000 [integer]
    @kwargs: doSILENT: silent the printing option of the function - Default is True [Boolean]
    """
    # Reading in the kwargs and performing some minor checks, so we have everything in the correct input format
    SPLINEknotpointsSPACING = kwargs.get('SPLINEknotpointsSPACING', np.array([0.2,0.25,1./3.])) #[pixel]
    SPLINEphaseSHIFT = kwargs.get('SPLINEphaseSHIFT', 0.01) #[pixel]
    SPLINEorder = kwargs.get('SPLINEorder', np.array([3],dtype='int32')) #[integer]
    SPLINEtckLENGTH = 'S' + str(np.int(kwargs.get('SPLINEstringLENGTH', 2000))) # [string] 
    silence = kwargs.get('doSILENT', True) #[boolean] --- diagnostic
    # Making sure you gave a numpy array for SPLINEorder with integers. If not, change it
    try:
      len(SPLINEorder)
    except:
      SPLINEorder = np.array([SPLINEorder],dtype='int32')
    maxNUMBERphaseSHIFTS = np.int(np.max(SPLINEknotpointsSPACING)/SPLINEphaseSHIFT)
    
    # Setting up the local matrices, for which we store the output. We multiply everthing with 1.e50 since we want the minimum BIC / AIC, and in case nothing is calculated, we want to avoid it and being able to trace it. -- If 1.e50 is too small for your usage, you are doing something horribly wrong. -- NOTE that the likelihood should be maximised, thus np.zeros
    matrixAIClocal, matrixBIClocal, likelihoodMATRIXlocal, matrixTCKlocal = np.ones((len(SPLINEknotpointsSPACING),maxNUMBERphaseSHIFTS,len(SPLINEorder)))*1.e50, np.ones((len(SPLINEknotpointsSPACING),maxNUMBERphaseSHIFTS,len(SPLINEorder)))*1.e50, np.zeros((len(SPLINEknotpointsSPACING),maxNUMBERphaseSHIFTS,len(SPLINEorder))), np.ones((len(SPLINEknotpointsSPACING),maxNUMBERphaseSHIFTS,len(SPLINEorder)), dtype=SPLINEtckLENGTH)
    
    # Sorting the positions
    coordSORTED, fluxSORTED = (list(t) for t in zip(*sorted(zip(coord, flux))))
    coordSORTED, fluxSORTED = np.array(coordSORTED), np.array(fluxSORTED)
    
    # Doing the fitting itself.
    for kk in range(len(SPLINEknotpointsSPACING)):	# Loop over the different sets of knotpoints
      NUMBERphaseSHIFTS = np.int(SPLINEknotpointsSPACING[kk]/SPLINEphaseSHIFT)
      for pp in range(NUMBERphaseSHIFTS): 		# Loop over the different phaseshifts for a given set of knotpoints
	#if pp != 0: # A shift
	  #coordKNOTPOINTS = np.arange(np.min(coord) + SPLINEknotpointsSPACING[kk] - pp*SPLINEphaseSHIFT, np.max(coord), SPLINEknotpointsSPACING[kk])
	#else: # No shift
	  #coordKNOTPOINTS = np.arange(np.min(coord) + SPLINEknotpointsSPACING[kk], np.max(coord), SPLINEknotpointsSPACING[kk])
	  
	## Check whether the last element is too far away or over the maximum value of that position. This happens roughly in ~0.1% of the calculations
	#if np.max(coord) - coordKNOTPOINTS[-1] > SPLINEknotpointsSPACING[kk]:
	  #coordKNOTPOINTS = np.append(coordKNOTPOINTS, coordKNOTPOINTS[-1] + SPLINEknotpointsSPACING[kk])
	#elif np.max(coord) - coordKNOTPOINTS[-1] <= 0:
	  #coordKNOTPOINTS = np.delete(coordKNOTPOINTS, -1)
	if pp !=0: # A shift
	  coordKNOTPOINTS = np.arange(np.min(coord) + SPLINEknotpointsSPACING[kk], np.max(coord), SPLINEknotpointsSPACING[kk]) + pp*SPLINEphaseSHIFT # New
	else: # No shift
	  coordKNOTPOINTS = np.arange(np.min(coord) + SPLINEknotpointsSPACING[kk], np.max(coord), SPLINEknotpointsSPACING[kk])
	# Check whether the last element is too far away or over the maximum value of that position. This happens roughly in ~50% of the calculations
	if np.max(coord) - coordKNOTPOINTS[-1] > SPLINEknotpointsSPACING[kk]:
	  coordKNOTPOINTS = np.append(coordKNOTPOINTS, coordKNOTPOINTS[-1] + SPLINEknotpointsSPACING[kk])
	elif np.max(coord) - coordKNOTPOINTS[-1] <= 0:
	  coordKNOTPOINTS = np.delete(coordKNOTPOINTS, -1)
	#  Check whether the first element is too far away from the minimum. This happens roughly in ~50% of the calculations
	if abs(np.min(coord) - coordKNOTPOINTS[0]) > SPLINEknotpointsSPACING[kk] + SPLINEphaseSHIFT:
	  coordKNOTPOINTS = np.append(coordKNOTPOINTS[0] - SPLINEknotpointsSPACING[kk], coordKNOTPOINTS)
	  
	for oo in range(len(SPLINEorder)): 		# Loop over the different orders of spline
	  NUMBERestimatedPARAMScoord = (len(coordKNOTPOINTS) + 1) * (SPLINEorder[oo] + 1)			# Preferred usage #NOTE +1 here, since the knotpoints DO NOT include beginning and ending
	  ###NUMBERestimatedPARAMScoord = np.int(1./SPLINEknotpointsSPACING[kk]) * (SPLINEorder[oo] + 1)	# Alternative, more prone to overfitting minor effects
	  TCKcoord, TCKerror = splineFIT(coordSORTED, fluxSORTED, SPLINEgiveKNOTPOINTS = True, SPLINEknotpoints = coordKNOTPOINTS, SPLINEorder = SPLINEorder[oo], doSILENT = silence)   
	  AICcoord, BICcoord, likelihoodCOORD = splineGOODNESSofFITandINFORMATIONCRITERION(coord, flux, TCKcoord, PARAMSdetermine=False, PARAMSestimated=NUMBERestimatedPARAMScoord, doSILENT = silence) # You should provide the full coord array, no rebinned arrays. We provide the number of estimated parameters.
	  matrixAIClocal[kk,pp,oo], matrixBIClocal[kk,pp,oo], likelihoodMATRIXlocal[kk,pp,oo], matrixTCKlocal[kk,pp,oo] = AICcoord, BICcoord, likelihoodCOORD, str(TCKcoord)
	  
    return matrixAIClocal, matrixBIClocal, likelihoodMATRIXlocal, matrixTCKlocal

def detrendTEMPpsfFULL(time, flux, xPOS, yPOS, **kwargs):
    """
    Master routine to perform the fitting between xPOS and yPOS vs flux for the temperature dependent PSF of BRITE photometry using the, full position arrays. The current technique is to perform a spline fit *without* the modulo 1 operation, for a large parameter space of knotpoints and / or spline orders.
    
    We provide a method to determine which set of knotpoints is actually the most appropriate, as this will influence your end result quite drastically. To do this, we use different information criteria and techniques to determine the goodness of fit. These criteratia are also used to determine which coordinate should be detrended first. - Hopefully it should not matter too much. - Since the position variations are most often 2D and thus cannot be represented by a PCA - I tried, we *ALWAYS* correct for both positions. In case there is no strong correlation, you shall see this in which group of knotpoints has been used.
    
    NOTE / UPDATE: we now call the routine detrendTEMPpsfCOORD to clean the coding up. As such, it will be easier to understand what is happening to a new user.
    
    NOTE We no longer assume that each pixel behaves similarly, hence no modulo 1 operation anymore along the positon axis - since it is actually a PSF effect instead of a pixel effect. Since missing datapoints, remaining outliers (bad clipping, go back a step!), ... will make your spline function do funny things, we will have to take additional care for the fitting and when we rebin the data.
    
    WARNING DANGER WARNING DANGER WARNING DANGER WARNING DANGER WARNING
    Because of the implementation of the spline fitting routines in scipy.interpolate (splrep function), you *need* to have at least sorted and ascending x values when performing the fitting process. 
    
    In addition, you can choose to have unique values along the x axis, but then you need to perform rebinning (=/= smoothing !!!). This is currently not implemented, but the routine rebinUNIQUE allows you to do this for two general parameter arrays. It provides an matrix as output, see the description of rebinUNIQUE, which should provide you with enough diagnostics to properly fit the spline to the rebinned values. At present, we only fit to the averages per bin, although one can implement an additional layer to decide if you should fit to the mean or the median value per position bin.
    WARNING DANGER WARNING DANGER WARNING DANGER WARNING DANGER WARNING
    
    NOTE Since we also want to compare the information criteria for x-position and y-position fits, we have to unify the number of fitting parameters somehow. If this is not done, the length of x (or y) will influence your BIC/AIC, since a larger lenght will give you more fitting parameters. As such, we use the spacing between consecutive knotpoints as a proxy for the number of fitting parameters, by calculating how many knotpoints we would have in a position region of length one. WARNING The phase shift might give you more knotpoints over the whole region. Unsure if we want to account for this! WARNING
    
    DANGER DANGER DANGER
    11/03/16: You need to subtract the mean of your flux before you apply any corrections. Otherwise, this is left in your 'correction flux' and will introduce mismatches between the different time bins, mess with your frequency diagrams, and introduce a whole range of effects you DO NOT want...
    
    15/03/16: Turned printing of found NaNs off at the moment. Too much spam, and you can now trace them using the diagCORRECTION.
    
    Returns: The correction you have to apply to the flux (which is a combination of the first and second correction), the TCK of the first correction, the TCK of the second correction, a diagnostic value tracing anything during the fitting process
    
    @param time: time measurements [d]
    @type time: numpy array of length N
    @param flux: flux measurements [adu]
    @type flux: numpy array of length N
    @param xPOS: CCD position measurements along x axis [pixel]
    @type xPOS: numpy array of length N
    @param yPOS: CCD position measurements along y axis [pixel]
    @type yPOS: numpy array of length N
    
    
    @return fluxCORRECTION: correction to apply to the flux
    @rtype: numpy array of length N
    @return: matrixTCKfirst[FIRSTcorrectionPARAMS][0]: TCK of the first correction, in string format
    @rtype: numpy.str_ (very annoying!)
    @return: matrixTCKsecond[SECONDcorrectionPARAMS][0]: TCK of the second correction, in string format
    @rtype: numpy.str_ (very annoying!)
    @return: diagCORRECTION: diagnostic to trace the corrections
    @rtype: numpy.int32
    
    @kwargs: SPLINEknotpointsSPACING: set of spacing for the different knotpoints - Default is np.array([0.2,0.25,1./3.]) [pixel]
    @kwargs: SPLINEphaseSHIFT: value to consider for the phase shift for the same sets of knotpoints - Default is 0.01 [pixel]
    @kwargs: SPLINEorder: order for the spline fits; does accept numpy arrays! - Default is np.array([3],dtype='int32')
    @kwargs: SPLINEtckLENGTH: number of bytes allocated for the to-string-converted TCK - Default is 2000 [integer]
    
    @kwargs: show_ME: Boolean to indicate if you want plotting at each possible step - Default is False [Boolean]
    @kwargs: show_FITS: Boolean to indicate if you want plotting after each bin fitting step - Default is False [Boolean]
    @kwargs: show_DIAG: Boolean to indicate if you want plotting when AIC != BIC - Default is False [Boolean]
    """
    
    # Reading in the kwargs and performing some minor checks, so we have everything in the correct input format
    SPLINEknotpointsSPACING = kwargs.get('SPLINEknotpointsSPACING', np.array([0.2,0.25,1./3.])) #[pixel]
    SPLINEphaseSHIFT = kwargs.get('SPLINEphaseSHIFT', 0.01) #[pixel]
    SPLINEorder = kwargs.get('SPLINEorder', np.array([3],dtype='int32')) #[integer]
    SPLINEtckLENGTH = 'S' + str(np.int(kwargs.get('SPLINEstringLENGTH', 2000))) # [string] 
    # Making sure you gave a numpy array for SPLINEorder with integers. If not, change it
    try:
      len(SPLINEorder)
    except:
      SPLINEorder = np.array([SPLINEorder],dtype='int32')
    maxNUMBERphaseSHIFTS = np.int(np.max(SPLINEknotpointsSPACING)/SPLINEphaseSHIFT)
    show_ME = kwargs.get('show_ME', False); show_FITS = kwargs.get('show_FITS', False); show_DIAG = kwargs.get('show_DIAG', False)

    # Subtracted the mean of your flux, ensuring to not change any offsets to your flux when calculating the corrections.
    flux = flux - np.mean(flux)
    
    # Preparing the two different diagnostic values
    diagCORRECTIONfirst, diagCORRECTIONsecond = np.int32(0), np.int32(0)
    
    """
	  First correction
    ----------------------------
    """
    
    # Setting up the master matrices, for which we store the output. We multiply everthing with 1.e50 since we want the minimum BIC / AIC, and in case nothing is calculated, we want to avoid it and being able to trace it. -- If 1.e50 is too small for your usage, you are doing something horribly wrong. -- NOTE that the likelihood should be maximised, thus np.zeros
    matrixAICfirst, matrixBICfirst, likelihoodMATRIXfirst, matrixTCKfirst = np.ones((2, len(SPLINEknotpointsSPACING),maxNUMBERphaseSHIFTS,len(SPLINEorder)))*1.e50, np.ones((2, len(SPLINEknotpointsSPACING),maxNUMBERphaseSHIFTS,len(SPLINEorder)))*1.e50, np.zeros((2, len(SPLINEknotpointsSPACING),maxNUMBERphaseSHIFTS,len(SPLINEorder))), np.ones((2, len(SPLINEknotpointsSPACING),maxNUMBERphaseSHIFTS,len(SPLINEorder)), dtype=SPLINEtckLENGTH)
    # Doing the correction using the x-position
    matrixAICfirst[0,:,:,:], matrixBICfirst[0,:,:,:], likelihoodMATRIXfirst[0,:,:,:], matrixTCKfirst[0,:,:,:] = detrendTEMPpsfCOORD(flux, xPOS, **kwargs)
    
    # Doing the correction using the y-position
    matrixAICfirst[1,:,:,:], matrixBICfirst[1,:,:,:], likelihoodMATRIXfirst[1,:,:,:], matrixTCKfirst[1,:,:,:] = detrendTEMPpsfCOORD(flux, yPOS, **kwargs)
    
    # Looking for the optimal solution to apply the correction. To do so, we use the AIC and BIC. In doubt, we also resort to the likelihood. -- The nanmin is important, since there might be some NaNs in the matrix --
    AICmin = np.where(matrixAICfirst==np.nanmin(matrixAICfirst))
    BICmin = np.where(matrixBICfirst==np.nanmin(matrixBICfirst))
    
    # Determine the number of NaNs (sum over a matrix with ones + the np.isnan function) -- DEBUGGING at the moment
    NUMBERofNANxAIC, NUMBERofNANxBIC = np.sum(np.ones_like(matrixAICfirst[0,:,:,:])[np.isnan(matrixAICfirst[0,:,:,:])]), np.sum(np.ones_like(matrixBICfirst[0,:,:,:])[np.isnan(matrixBICfirst[0,:,:,:])])
    NUMBERofNANyAIC, NUMBERofNANyBIC = np.sum(np.ones_like(matrixAICfirst[1,:,:,:])[np.isnan(matrixAICfirst[1,:,:,:])]), np.sum(np.ones_like(matrixBICfirst[1,:,:,:])[np.isnan(matrixBICfirst[1,:,:,:])])
    
    if NUMBERofNANxAIC != 0 or  NUMBERofNANxBIC != 0 or NUMBERofNANyAIC != 0 or  NUMBERofNANyBIC != 0:
      ####print '\tNumber of NaNs found for the x-axis in the AIC = {:4.0f} and in the BIC = {:4.0f}'.format(NUMBERofNANxAIC, NUMBERofNANxBIC)
      ####print '\tNumber of NaNs found for the y-axis in the AIC = {:4.0f} and in the BIC = {:4.0f}'.format(NUMBERofNANyAIC, NUMBERofNANyBIC)
      ####print '\tAlso providing the minima in the AIC and BIC matrix {:.5e} and {:.5e}, respectively'.format(np.nanmin(matrixAICfirst), np.nanmin(matrixBICfirst))
      
      # Changing all NaNs to 1.e75
      matrixAICfirst[np.isnan(matrixAICfirst)], matrixBICfirst[np.isnan(matrixBICfirst)] = 1.e75, 1.e75
      
      # Counting all values that are 1.e50, comparing them to the number of NaNs and the total number of entries
      matrixAICfirstUNIQUE, matrixAICfirstUNIQUEcounts = np.unique(matrixAICfirst, return_counts=True)
      matrixBICfirstUNIQUE, matrixBICfirstUNIQUEcounts = np.unique(matrixBICfirst, return_counts=True)
      
      NUMBERunusedAIC, NUMBERunusedBIC = matrixAICfirstUNIQUEcounts[np.where(matrixAICfirstUNIQUE == 1.e50)][0], matrixBICfirstUNIQUEcounts[np.where(matrixBICfirstUNIQUE == 1.e50)][0]
      
      ####print '\tFor the AIC:\n\t\tTotal entries = {:4.0f}\n\t\tUnused entries = {:4.0f}\n\t\tNaN entries = {:4.0f}'.format(matrixAICfirst.size, NUMBERunusedAIC, NUMBERofNANxAIC+NUMBERofNANyAIC)
      ####print '\tFor the BIC:\n\t\tTotal entries = {:4.0f}\n\t\tUnused entries = {:4.0f}\n\t\tNaN entries = {:4.0f}\n'.format(matrixBICfirst.size, NUMBERunusedBIC, NUMBERofNANxBIC+NUMBERofNANyBIC)
      
      # Diagnostic values
      diagCORRECTIONfirst += 10
      
    # Ensure that both the AIC and BIC point towards the same optimal solution.
    if AICmin != BICmin:
      print bcolors.WARNING + '\tWARNING: AIC and BIC point towards a different coordinate / set of knotpoints\n\tProviding more information now...' + bcolors.ENDC
      print bcolors.WARNING + '\tWARNING: AIC states {:1.0f} and BIC states {:1.0f} for the coordinate, where 0 = x and 1 = y'.format(AICmin[0][0], BICmin[0][0]) + bcolors.ENDC
      print bcolors.WARNING + '\tWARNING: AIC states {:1.3f} and BIC states {:1.3f} for the spacing of the knotpoints'.format(SPLINEknotpointsSPACING[AICmin[1][0]], SPLINEknotpointsSPACING[BICmin[1][0]]) + bcolors.ENDC
      print bcolors.WARNING + '\tWARNING: AIC states {:1.3f} and BIC states {:1.3f} for the phase shift of those sets of knotpoints'.format(AICmin[2][0]*SPLINEphaseSHIFT, BICmin[2][0]*SPLINEphaseSHIFT) + bcolors.ENDC
      print bcolors.WARNING + '\tWARNING: AIC states {:1.3f} and BIC states {:1.3f} for the order of the spline '.format(SPLINEorder[AICmin[3][0]], SPLINEorder[BICmin[3][0]]) + bcolors.ENDC
      print bcolors.WARNING + '\tWARNING: AIC states a minimum of {:.5e}, compared to {:.5e} for the BIC optimal solution'.format(matrixAICfirst[AICmin][0], matrixAICfirst[BICmin][0]) + bcolors.ENDC
      print bcolors.WARNING + '\tWARNING: BIC states a minimum of {:.5e}, compared to {:.5e} for the AIC optimal solution'.format(matrixBICfirst[BICmin][0], matrixBICfirst[AICmin][0]) + bcolors.ENDC
      print bcolors.WARNING + '\tWARNING: Likelihoods are {:.5e} and {:.5e}, for AIC and BIC respectively'.format(likelihoodMATRIXfirst[AICmin][0], likelihoodMATRIXfirst[BICmin][0]) + bcolors.ENDC
      
      # We follow the suggestion of the likelihood in case the AIC and BIC do not agree for the optimal solution.
      likelihoodAIC, likelihoodBIC = likelihoodMATRIXfirst[AICmin][0], likelihoodMATRIXfirst[BICmin][0] # The [0] is needed because of the np.where() function.
      if likelihoodAIC > likelihoodBIC:
        print bcolors.FAIL + '\tTaking the AIC as best solution, indicated by the LIKELIHOOD' + bcolors.ENDC
        print bcolors.FAIL + '\tTaking the AIC as best solution, indicated by the LIKELIHOOD' + bcolors.ENDC
        FIRSTcorrectionPARAMS = AICmin
        # Diagnostic values
        diagCORRECTIONfirst += 1
      elif likelihoodAIC < likelihoodBIC:
        print bcolors.FAIL + '\tTaking the **BIC** as best solution, indicated by the LIKELIHOOD' + bcolors.ENDC
        print bcolors.FAIL + '\tTaking the **BIC** as best solution, indicated by the LIKELIHOOD' + bcolors.ENDC
        FIRSTcorrectionPARAMS = BICmin
        # Diagnostic values
        diagCORRECTIONfirst += 2
        
      # Provide some diagnostics here for debugging, so we can actually study if the likelihood is a good suggestion.
      if show_ME or show_DIAG:
	TCKfirstAIC = reconvertTCKfromSTRING(matrixTCKfirst[AICmin][0])
	TCKfirstBIC = reconvertTCKfromSTRING(matrixTCKfirst[BICmin][0])
	
	if AICmin[0][0] == 0 & BICmin[0][0] == 0:
	  PLOTdetrendTEMPpsfDIAGinformCRIT(flux, xPOS, TCKfirstAIC, xPOS, TCKfirstBIC, 'xx')  
	elif AICmin[0][0] == 0 & BICmin[0][0] == 1:
	  PLOTdetrendTEMPpsfDIAGinformCRIT(flux, xPOS, TCKfirstAIC, yPOS, TCKfirstBIC, 'xy')
	elif AICmin[0][0] == 1 & BICmin[0][0] == 0:
	  PLOTdetrendTEMPpsfDIAGinformCRIT(flux, yPOS, TCKfirstAIC, xPOS, TCKfirstBIC, 'yx')
	elif AICmin[0][0] == 1 & BICmin[0][0] == 1:
	  PLOTdetrendTEMPpsfDIAGinformCRIT(flux, yPOS, TCKfirstAIC, yPOS, TCKfirstBIC, 'yy')   
	pl.show()

    elif AICmin == BICmin:
      FIRSTcorrectionPARAMS = AICmin # Should not matter which one you take, since AICmin == BICmin
    
    # Perform the first correction using the optimal settings for the fit
    if FIRSTcorrectionPARAMS[0][0] == 0: # x-position was chosen for the first correction
      TCKfirstCORRECTION = reconvertTCKfromSTRING(matrixTCKfirst[FIRSTcorrectionPARAMS][0])
      FLUXfirstCORRECTION = scInterp.splev(xPOS, TCKfirstCORRECTION)
      # Diagnostic values
      diagCORRECTIONfirst += 100; diagCORRECTIONsecond += 200   
    if FIRSTcorrectionPARAMS[0][0] == 1: # y-position was chosen for the first correction
      TCKfirstCORRECTION = reconvertTCKfromSTRING(matrixTCKfirst[FIRSTcorrectionPARAMS][0])
      FLUXfirstCORRECTION = scInterp.splev(yPOS, TCKfirstCORRECTION)    
      # Diagnostic values
      diagCORRECTIONfirst += 200; diagCORRECTIONsecond += 100 
    
    """
	  Second correction
    -----------------------------
    """
    # We do not set up the matrices for the output storage of the fitting, since these were already created by the detrendTEMPpsfCOORD (and are only 3D instead of 4D, since there is only one coordinate to fit)
    if FIRSTcorrectionPARAMS[0][0] == 0:# Doing the correction using the y-position, since the x-position was chosen for the first correction
      matrixAICsecond, matrixBICsecond, likelihoodMATRIXsecond, matrixTCKsecond = detrendTEMPpsfCOORD(flux - FLUXfirstCORRECTION, yPOS, **kwargs)
    if FIRSTcorrectionPARAMS[0][0] == 1:# Doing the correction using the y-position, since the x-position was chosen for the first correction
      matrixAICsecond, matrixBICsecond, likelihoodMATRIXsecond, matrixTCKsecond = detrendTEMPpsfCOORD(flux - FLUXfirstCORRECTION, xPOS, **kwargs)
    
    # Looking for the optimal solution to apply the correction. To do so, we use the AIC and BIC. In doubt, we also resort to the likelihood. -- The nanmin is important, since there might be some NaNs in the matrix --
    AICmin2 = np.where(matrixAICsecond==np.nanmin(matrixAICsecond))
    BICmin2 = np.where(matrixBICsecond==np.nanmin(matrixBICsecond))
    
    # Determine the number of NaNs (sum over a matrix with ones + the np.isnan function) -- DEBUGGING at the moment
    NUMBERofNANAIC, NUMBERofNANBIC = np.sum(np.ones_like(matrixAICsecond[:,:,:])[np.isnan(matrixAICsecond[:,:,:])]), np.sum(np.ones_like(matrixBICsecond[:,:,:])[np.isnan(matrixBICsecond[:,:,:])])
    
    if NUMBERofNANAIC != 0 or  NUMBERofNANBIC != 0:
      ####print '\t\tNumber of NaNs found for the second position axis in the AIC = {:4.0f} and in the BIC = {:4.0f}'.format(NUMBERofNANAIC, NUMBERofNANBIC)
      ####print '\t\tAlso providing the minima in the AIC and BIC matrix {:.5e} and {:.5e}, respectively'.format(np.nanmin(matrixAICsecond), np.nanmin(matrixBICsecond))
      # Changing all NaNs to 1.e75
      matrixAICsecond[np.isnan(matrixAICsecond)], matrixBICsecond[np.isnan(matrixBICsecond)] = 1.e75, 1.e75
      
      # Counting all values that are 1.e50, comparing them to the number of NaNs and the total number of entries
      matrixAICsecondUNIQUE, matrixAICsecondUNIQUEcounts = np.unique(matrixAICsecond, return_counts=True)
      matrixBICsecondUNIQUE, matrixBICsecondUNIQUEcounts = np.unique(matrixBICsecond, return_counts=True)
      
      NUMBERunusedAIC, NUMBERunusedBIC = matrixAICsecondUNIQUEcounts[np.where(matrixAICsecondUNIQUE == 1.e50)][0], matrixBICsecondUNIQUEcounts[np.where(matrixBICsecondUNIQUE == 1.e50)][0]
      
      ####print '\t\tFor the AIC:\n\t\t\tTotal entries = {:4.0f}\n\t\t\tUnused entries = {:4.0f}\n\t\t\tNaN entries = {:4.0f}'.format(matrixAICsecond.size, NUMBERunusedAIC, NUMBERofNANAIC)
      ####print '\t\tFor the BIC:\n\t\t\tTotal entries = {:4.0f}\n\t\t\tUnused entries = {:4.0f}\n\t\t\tNaN entries = {:4.0f}\n'.format(matrixBICsecond.size, NUMBERunusedBIC, NUMBERofNANBIC)
      
      # Diagnostic values
      diagCORRECTIONsecond += 10
    
    #DEBUGGING
    ###if len(AICmin2[0]) != 1:
      ###print matrixAICsecond[AICmin2], matrixBICsecond[BICmin2]
      ###print '\t\tFor the AIC:\n\t\t\tTotal entries = {:4.0f}\n\t\t\tUnused entries = {:4.0f}\n\t\t\tNaN entries = {:4.0f}'.format(matrixAICsecond.size, NUMBERunusedAIC, NUMBERofNANAIC)
      ###print '\t\tFor the BIC:\n\t\t\tTotal entries = {:4.0f}\n\t\t\tUnused entries = {:4.0f}\n\t\t\tNaN entries = {:4.0f}\n'.format(matrixBICsecond.size, NUMBERunusedBIC, NUMBERofNANBIC)
      
    
    # Ensure that both the AIC and BIC point towards the same optimal solution.
    if AICmin2 != BICmin2:
      print bcolors.WARNING + '\t\tWARNING: AIC and BIC point towards a different set of knotpoints\n\t\tProviding more information now...' + bcolors.ENDC
      print bcolors.WARNING + '\t\tWARNING: AIC states {:1.3f} and BIC states {:1.3f} for the spacing of the knotpoints'.format(SPLINEknotpointsSPACING[AICmin2[0][0]], SPLINEknotpointsSPACING[BICmin2[0][0]]) + bcolors.ENDC
      print bcolors.WARNING + '\t\tWARNING: AIC states {:1.3f} and BIC states {:1.3f} for the phase shift of those sets of knotpoints'.format(AICmin2[1][0]*SPLINEphaseSHIFT, BICmin2[1][0]*SPLINEphaseSHIFT) + bcolors.ENDC
      print bcolors.WARNING + '\t\tWARNING: AIC states {:1.3f} and BIC states {:1.3f} for the order of the spline '.format(SPLINEorder[AICmin2[2][0]], SPLINEorder[BICmin2[2][0]]) + bcolors.ENDC
      print bcolors.WARNING + '\t\tWARNING: AIC states a minimum of {:.5e}, compared to {:.5e} for the BIC optimal solution'.format(matrixAICsecond[AICmin2][0], matrixAICsecond[BICmin2][0]) + bcolors.ENDC
      print bcolors.WARNING + '\t\tWARNING: BIC states a minimum of {:.5e}, compared to {:.5e} for the AIC optimal solution'.format(matrixBICsecond[BICmin2][0], matrixBICsecond[AICmin2][0]) + bcolors.ENDC
      print bcolors.WARNING + '\t\tWARNING: Likelihoods are {:.5e} and {:.5e}, for AIC and BIC respectively'.format(likelihoodMATRIXsecond[AICmin2][0], likelihoodMATRIXsecond[BICmin2][0]) + bcolors.ENDC
      
      # We follow the suggestion of the likelihood in case the AIC and BIC do not agree for the optimal solution.
      likelihoodAIC, likelihoodBIC = likelihoodMATRIXsecond[AICmin2][0], likelihoodMATRIXsecond[BICmin2][0] # The [0] is needed because of the np.where() function.
      if likelihoodAIC > likelihoodBIC:
        print bcolors.FAIL + '\t\tTaking the AIC as best solution, indicated by the LIKELIHOOD' + bcolors.ENDC
        print bcolors.FAIL + '\t\tTaking the AIC as best solution, indicated by the LIKELIHOOD' + bcolors.ENDC
        SECONDcorrectionPARAMS = AICmin2
        # Diagnostic values
        diagCORRECTIONsecond += 1
      elif likelihoodAIC < likelihoodBIC:
        print bcolors.FAIL + '\t\tTaking the **BIC** as best solution, indicated by the LIKELIHOOD' + bcolors.ENDC
        print bcolors.FAIL + '\t\tTaking the **BIC** as best solution, indicated by the LIKELIHOOD' + bcolors.ENDC
        SECONDcorrectionPARAMS = BICmin2
        # Diagnostic values
        diagCORRECTIONsecond += 2
        
      # Provide some diagnostics here for debugging, so we can actually study if the likelihood is a good suggestion.
      if show_ME or show_DIAG:
	TCKsecondAIC = reconvertTCKfromSTRING(matrixTCKsecond[AICmin2][0])
	TCKsecondBIC = reconvertTCKfromSTRING(matrixTCKsecond[BICmin2][0])
	
	if FIRSTcorrectionPARAMS[0][0] == 0:
	  PLOTdetrendTEMPpsfDIAGinformCRIT(flux - FLUXfirstCORRECTION, yPOS, TCKsecondAIC, yPOS, TCKsecondBIC, 'yy')  
	elif FIRSTcorrectionPARAMS[0][0] == 1:
	  PLOTdetrendTEMPpsfDIAGinformCRIT(flux - FLUXfirstCORRECTION, xPOS, TCKsecondAIC, xPOS, TCKsecondBIC, 'xx')
	pl.show()  
      
    elif AICmin2 == BICmin2:
      SECONDcorrectionPARAMS = AICmin2 # Should not matter which one you take, since AICmin2 == BICmin2
    
    # Perform the second correction using the optimal settings for the fit
    if FIRSTcorrectionPARAMS[0][0] == 0: # y-position *has* to be taken for the second correction      
      TCKsecondCORRECTION = reconvertTCKfromSTRING(matrixTCKsecond[SECONDcorrectionPARAMS][0])
      FLUXsecondCORRECTION = scInterp.splev(yPOS, TCKsecondCORRECTION)    
      
      # Perform plotting
      if show_ME or show_FITS:
        PLOTdetrendTEMPpsfFULL(time, flux, xPOS, TCKfirstCORRECTION, yPOS, TCKsecondCORRECTION, 'xy', **kwargs) 
        pl.show()
    if FIRSTcorrectionPARAMS[0][0] == 1: # x-position *has* to be taken for the second correction
      TCKsecondCORRECTION = reconvertTCKfromSTRING(matrixTCKsecond[SECONDcorrectionPARAMS][0])
      FLUXsecondCORRECTION = scInterp.splev(xPOS, TCKsecondCORRECTION)
      
      # Perform plotting
      if show_ME or show_FITS:
        PLOTdetrendTEMPpsfFULL(time, flux, yPOS, TCKfirstCORRECTION, xPOS, TCKsecondCORRECTION, 'yx', **kwargs)
        pl.show()
    
    # Combining both diagnostic values to one integer. NOTE: At the moment of writing (15/03/16), each diagCORRECTIONxxx corresponds to three digits, hence we multiple diagCORRECTIONfirst with 1000.
    diagCORRECTION = diagCORRECTIONfirst * 1000 + diagCORRECTIONsecond
    
    print 'Minor debugging diagnostics provided for your convenience'
    print '1 - knotpoint spacing = {:1.2f}; phaseshift spacing = {:1.2f}'.format(SPLINEknotpointsSPACING[AICmin[1][0]], AICmin[2][0]*SPLINEphaseSHIFT)
    print '2 - knotpoint spacing = {:1.2f}; phaseshift spacing = {:1.2f}'.format(SPLINEknotpointsSPACING[AICmin2[0][0]], AICmin2[1][0]*SPLINEphaseSHIFT)
    
    return FLUXfirstCORRECTION + FLUXsecondCORRECTION, matrixTCKfirst[FIRSTcorrectionPARAMS][0], matrixTCKsecond[SECONDcorrectionPARAMS][0], diagCORRECTION    
