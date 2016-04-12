# -*- coding: utf-8 -*-
"""
General routines to perform a pca between two dimensions. At the moment, only a PCA routine between the x and y position is implemented. It is based on the work by Vanderburg and Johnson 2014, for the detrending of K2 data for the moving centroid position.

At present, I do not recommend to use this routine for the arclength. Moreover, the effect of the temperature dependent PSF changes on the position cannot be represented by a PCA!
    
Last update 17 March 2016

@author: Bram Buysschaert
"""

#===============================================================================
# 				Packages
#===============================================================================
import numpy as np

import scipy.interpolate as scInterp

from BRITE_decor.fitting.splinefit import splineFIT, splineGOODNESSofFITandINFORMATIONCRITERION
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


def positionPCA(time, xPOS, yPOS, **kwargs):
    """
    Routine to perform the PCA (principal component analysis) of the CCD positions. It converts a 2D motion to a dominant 1D motion. The covariance matrix is used to perform the PCA (WARNING a change in units of the position might have a strong influence on this matrix and thus the final result!).
    
    TODO
    Include a proper way to calculate the arclength along this 1D motion (instead of the 10th order polynomial now used in the K2 data reduction). To this end, the time measurments are already included as an input requirement.
    TODO
    
    NOTE: Currently (24/01/16) we are working with the numpy.cov function. This calculations the unnormalized covariance matrix. If you use the numpy.corrcoef routine, you would determine the normalized *correlation* matrix. It should (hopefully) not matter too much which routine you should call.
    DANGER: Within the same sidenote, we intend to warn the user that only a linear trend between the two parameters are properly determined by a correlation / covariance matrix. Anything else, like the periodic variations found here, should be treated with care.
    
    Returns: The dominant 1D motion and the less domination motion found by the PCA, unless ARCLENGTHoutput = True, than we return the arclength along the rotated motion.
    
    @param time: time measurements [d]
    @type time: numpy array of length N
    @param xPOS: CCD x position measurements [pixel]
    @type xPOS: numpy array of length N
    @param yPOS: CCD y position measurements [pixel]
    @type yPOS: numpy array of length N
        
    @return: dxPOSprime: the most dominant motion, recovered from the PCA
    @rtype: numpy array of length N
    @return: dyPOSprime: the non-dominant motion, recovered from the PCA
    @rtype: numpy array of length N
    @return: arclength: the arclength along the rotated position, recovered from the PCA (and a vague function)
    @rtype: numpy array of length N
    
    @kwargs: ARCLENGTHoutput: calculate the arclength for the dominant motion of the PCA, instead of this altered (rotated) coordinate - Default is False [Boolean]
    @kwargs: CLIPiteration: maximum number of iterations permitted to remove all *three*-sigma outliers - Default is 5 [Integer]
    """
    doARCLENGTH = kwargs.get('ARCLENGTHoutput', False)
    
    # Mean subtraction
    dxPOS, dyPOS = xPOS - np.mean(xPOS), yPOS - np.mean(yPOS)
    
    # Searching for the largest eigenvalue of the covariance matrix
    covarianceMATRIX = np.cov(dxPOS, dyPOS)
    covarianceMATRIXeigenVals, covarianceMATRIXeigenVecs = np.linalg.eig(covarianceMATRIX)
    ###print 'The two eigenvalues for the covariance matrix are\t{:1.5f} and {:1.5f}'.format(covarianceMATRIXeigenVals[0],covarianceMATRIXeigenVals[1])
  
    # Rotate the centroid positions to a new coordinate system, where x' and y' are along the eigenvectors
    #And x' is along the eigenvector with the largest eigenvalue
    if covarianceMATRIXeigenVals[0] > covarianceMATRIXeigenVals[1]:
      rotationMATRIX = covarianceMATRIXeigenVecs
    else:
      rotationMATRIX = covarianceMATRIXeigenVecs[:,::-1]
      
    dMATRIX = np.matrix((dxPOS, dyPOS))
    dMATRIXrotated = np.dot(dMATRIX.T,rotationMATRIX)
    
    dxPOSprime = np.zeros_like(dxPOS)
    dyPOSprime = np.zeros_like(dyPOS)
    
    #Issues with the matrix environment, so we copy elementwise the matrix to arrays. NOTE might need some better way to solve this or might need optimization.
    for ii in range(len(dxPOS)):
      dxPOSprime[ii] = dMATRIXrotated[ii,0]
      dyPOSprime[ii] = dMATRIXrotated[ii,1]
    
    if not(doARCLENGTH):
      return dxPOSprime, dyPOSprime
    else:
      # Performing the calculations for the arclength representation of the 'rotated' dataset.
      # We follow here the approach by Vanderburg+2014 (however, the fitting method you use for the relation between dxPOSprime and dyPOSprime is not mentioned)
      maxITER = kwargs.get('CLIPiteration',5)
      #10th order polynomial fit to dxPOSprime and dyPOSprime
      #outliers are taken into account
      dxPOSprimeToFit, dyPOSprimeToFit, idx3Sigma = np.copy(dxPOSprime), np.copy(dyPOSprime), np.array([],dtype='int32')  
      for ii in range(maxITER):
	poly10fit = np.polyfit(dxPOSprimeToFit,dyPOSprimeToFit,10)
	residualPoly10fit = dyPOSprimeToFit - np.polyval(poly10fit,dxPOSprimeToFit)
	residualPoly10fitMean, residualPoly10fitSTD = np.mean(residualPoly10fit), np.std(residualPoly10fit)
	
	idxToClip = np.where(np.abs(residualPoly10fit) > 3.*residualPoly10fitSTD)[0]
	for ll in range(len(idxToClip)):
	  idx3Sigma = np.append(idx3Sigma, np.where(dxPOSprime==dxPOSprimeToFit[idxToClip[ll]])[0]) #NOTE look in original array
	
	#Checking if there were *three*-sigma outliers
	if len(idxToClip) > 0: 
	  #Clipping out the outliers
	  dxPOSprimeToFit = np.delete(dxPOSprimeToFit, idxToClip)
	  dyPOSprimeToFit = np.delete(dyPOSprimeToFit, idxToClip)
	else:
	  break
	
      #calculate the derivative of the polynomial  
      poly10fitDerivative = np.polyder(poly10fit,m=1)
      
      #Sort dxPOSprime, in order to calculate the arclength
      if (np.__version__!='1.8.1') & (np.__version__!='1.9.2'): #OLD version
	sortRankdxPOSprime = np.argsort(dxPOSprime)
	dxPOSprimeSorted = dxPOSprime[:,sortRankdxPOSprime]
	dyPOSprimeSorted = dyPOSprime[:,sortRankdxPOSprime]
	timeSorted = time[:,sortRankdxPOSprime]
      else: #GENERAL version and most likely faster
	dxPOSprimeSorted, dyPOSprimeSorted, timeSorted = (list(t) for t in zip(*sorted(zip(dxPOSprime, dyPOSprime, time))))
	dxPOSprimeSorted, dyPOSprimeSorted, timeSorted = np.array(dxPOSprimeSorted), np.array(dyPOSprimeSorted), np.array(timeSorted)
      
      """
      Arclength calculation
      ---------------------
      """  
      #Make function of the integral for the arclength, in order to solve it as an explicit integral
      integral = lambda xx: np.sqrt(1. + (np.polyval(poly10fitDerivative,xx))**2.)
      
      #Calculate the integral for the arclength
      sSorted = np.zeros_like(dxPOSprimeSorted)
      for ii in range(len(sSorted)):
	sSorted[ii] = scInteg.quad(integral,dxPOSprimeSorted[0],dxPOSprimeSorted[ii])[0]
      
      #Resort the arclength to match the monotone increasing time array
      if (np.__version__!='1.8.1') & (np.__version__!='1.9.2'): #OLD version
	sortRankTime = np.argsort(timeSorted)
	arclength = sSorted[:,sortRankTime]
      else: #GENERAL version and most likely faster
	timeCopy, arclength = (list(t) for t in zip(*sorted(zip(timeSorted, sSorted))))
	timeCopy, arclength = np.array(timeCopy), np.array(arclength)
      
      
      return arclength