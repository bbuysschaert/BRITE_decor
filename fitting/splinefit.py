# -*- coding: utf-8 -*-
"""
Routines to perform general spline fitting and have some diagnostics
    
Last update 16 March 2016

@author: Bram Buysschaert
"""

#===============================================================================
# 				Packages
#===============================================================================
import numpy as np

import scipy.interpolate as scInterp
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
def reconvertTCKfromSTRING(TCKstring, **kwargs):
    """
    Routine to reconvert a TCK from the scInterp.splrep, which was saved as a string, back to an useful tuple.
    
    The following properties are known of the string:
    -The first 8 elements should be '(array(['
    -The second to last element should contain the order of the spline, i.e. "k)']". --> Can stop at element -3
    -Once the first numpy array is retrieved, we should find a '], array(['
    -There are '\n' saved in the string. --> Need to be removed
    
    
    NOTE: Any bug here, means that your predefined length for the string was too short (standard 2000 characters)
    
    Returns: The TCK, but as a tuple.
    
    @param TCKstring: string of the tck
    @type TCKstring: numpy.string_
    
    @return TCK: the reconverted tck
    @rtype: tuple
    """
    
    # Convert the (annoying) numpy.string_ to a 'classical' string
    # NOTE this has some extra effects
    TCKstring = str(TCKstring)
    
    # Get the rough arrays out already
    TCKfirstARRAYstr, TCKsecondARRAYstr = TCKstring.split('array([')[1:]
    TCKfirstARRAY, TCKsecondARRAY = [], []
    
    # Clean the first array
    TCKfirstARRAYstr = TCKfirstARRAYstr.split(']')[0].split(',')
    for ss in range(len(TCKfirstARRAYstr)):
      TCKarrayELEMENT = TCKfirstARRAYstr[ss].split(' ')
      for ee in range(len(TCKarrayELEMENT)):
	  try:
	    TCKfirstARRAY.append(float(TCKarrayELEMENT[ee]))
	  except:
	    pass
    # Clean the second array
    TCKsecondARRAYstr = TCKsecondARRAYstr.split(']')[0].split(',')
    for ss in range(len(TCKsecondARRAYstr)):
      TCKarrayELEMENT = TCKsecondARRAYstr[ss].split(' ')
      for ee in range(len(TCKarrayELEMENT)):
	  try:
	    TCKsecondARRAY.append(float(TCKarrayELEMENT[ee]))
	  except:
	    pass
    # Get the order of the spline
    TCKorder = int(TCKstring.split('array([')[2].split(']')[1].split(',')[1].split(')')[0])
    
    TCK = (np.array(TCKfirstARRAY), np.array(TCKsecondARRAY), TCKorder)
    return TCK
  
def splineFIT(param1, param2, **kwargs):
    """
    Spline representation for a given dataset. It is assumed that you can represent param2 as a function of param1. In the end, this is a simple decorator function to call the scipy.interpolate.splrep package to perform the spline representation.

    Returns: The optimal tck tuple from the spline fit ((t,c,k) is a tuple containing the vector of knots, the B-spline coefficients, and the degree of the spline)
    
    @param param1: parameter 1 values 
    @type param1: numpy array of length N
    @param param2: parameter 2 values 
    @type param2: numpy array of length N
    
    @return: tckSPLINE
    @rtype: tuple
    
    @kwargs: SPLINEgiveKNOTPOINTS: provide user specified knotpoints for the spline representation - Default is False [Boolean]
    @kwargs: SPLINEknotpointsSPACING: the user specified spacing between consecutive knotpoints - Default is 1.0 [units of param1]
    @kwargs: SPLINEknotpoints: the user specified knotpoints for the spline representation; should be used with SPLINEgiveKNOTPOINTS - Default is None (e.g. pop) [param1 units]
    @kwargs: SPLINEperiodic: state if the spline representation should be a periodic function or not - Default is False [Boolean]
    @kwargs: SPLINEorder: state which order the spline representation should be; preferibly odd order - default is 3 [integer]
    @kwargs: doSILENT: silent the printing option of the function - Default is True [Boolean]
    """
    
    manualKNOTPOINT = kwargs.get('SPLINEgiveKNOTPOINTS', False) #[Boolean]
    orderSPLINE = np.int(kwargs.get('SPLINEorder', 3)) #[integer]
    periodicSPLINE = kwargs.get('SPLINEperiodic', False) #[Boolean]
    doSILENT = kwargs.get('doSILENT', True) #[Boolean]
    
    if not(manualKNOTPOINT):
      spacingKNOTPOINTS = kwargs.get('SPLINEknotpointsSPACING', 1.0) #[units of param1]
      knotpointSPLINE = np.arange(min(param1), max(param1), spacingKNOTPOINTS) # [units of param1]
    else:
      knotpointSPLINE = kwargs.pop('SPLINEknotpoints')
    
    # WARNING your first element and last element of the knotpointSPLINE should be larger than the minimum and smaller than the maximum of param2, respectively. If not, splrep WILL complain.
    if not(periodicSPLINE):
      tckSPLINE, fpSPLINE, ierSPLINE, msgSPLINE  = scInterp.splrep(param1, param2, t=knotpointSPLINE[:], k=orderSPLINE, full_output=1)
    elif periodicSPLINE:
      tckSPLINE, fpSPLINE, ierSPLINE, msgSPLINE = scInterp.splrep(param1, param2, t=knotpointSPLINE[:], k=orderSPLINE, per=1, full_output=1)    
    
    if (ierSPLINE != 0) and not(doSILENT):
      print 'The spline fitting routine produced an error, which is the following:'
      print ierSPLINE, msgSPLINE
      if ierSPLINE == 10:
	print 'You should not use this fit, it is most likely FLAT ...'
    
    return tckSPLINE, ierSPLINE

def splineGOODNESSofFITandINFORMATIONCRITERION(param1, param2, tck, **kwargs): #Used  (09/03/16)
    """
    For a given spline fit which captures the behaviour of param2 with varying param1 do:
    - determine the goodness of fit using the residuals, a likelihood function
    - determine the information criterion for the fit (both Akaike and Bayesian information criterion) using the likelihood
    
    The used likelihood is slightly adapted for negative values (using the abs() function) from Duvall & Harvey 1986, Anderson+ 1990 and seen as Eq. (9) in Buysschaert, Beck, et al. 2016.
    
    WARNING We calculate the number of estimated parameters within a spline as followed:
    (number of fitting regions) * (order of spline + 1)
    
    NOTE:
    AIC is better in situations when a false negative finding would be considered more misleading than a false positive, and BIC is better in situations where a false positive is as misleading as, or more misleading than, a false negative
    NOTE:
    What if AIC and BIC do not point to the same model:
    http://stats.stackexchange.com/questions/577/is-there-any-reason-to-prefer-the-aic-or-bic-over-the-other
    TODO:
    Implement a Pearson's chi square
    
    REMEMBER
    The np.log is the natural logarithm.
    
    
    Returns: The AIC, BIC, and the likelihood of the model characterised by tck
    
    @param param1: parameter 1 values 
    @type param1: numpy array of length N
    @param param2: parameter 2 values 
    @type param2: numpy array of length N
    @param tck: tck of the scipy.interpolate.splrep function
    @type tck: tuple
    
    @return: AIC: Akaike information criterion
    @rtype: numpy.float
    @return: BIC: Bayesian information criterion
    @rtype: numpy.float
    @return: likelihood: likelihood of the fit
    @rtype: numpy.float
    
    @kwargs: fullOUTPUT: print all determined diagnostic values for the goodness of fit - Default is False [Boolean]
    @kwargs: PARAMSdetermine: use the tck to determine the number of estimated parameters in the fit - Default is True [Boolean]
    @kwargs: PARAMSestimated: provide the number of estimated parameters during the fit - is popped; so no default [float]
    """
    doSILENT = kwargs.get('doSILENT', True) #[Boolean]
    
    residuals = param2 - scInterp.splev(param1, tck)
    RSS = np.sum(residuals**2.)
    residualsSTD = np.std(residuals)
    likelihood = np.sum(np.log(np.abs(scInterp.splev(param1, tck))) + np.abs(param2 / scInterp.splev(param1, tck))) #NOTE absolute values of the fit, since param2 (and its representation) can be lower than 0, which is *not* ideal for a logarithm.
    
    
    if kwargs.get('PARAMSdetermine', True):
      estimatedPARAMS = (len(tck[0])/2 - 1) * (tck[-1] + 1) # (number of knotpoint REGIONS) * (order of the spline + 1) #NOTE -1 here, since the knotpoints include beginning and ending.
    else:
      estimatedPARAMS = kwargs.pop('PARAMSestimated')
    
    AIC = estimatedPARAMS * np.log(len(param1)) + len(param1) * np.log(likelihood)
    BIC = 2. * estimatedPARAMS + len(param1) * np.log(likelihood)
    
    if kwargs.get('fullOUTPUT', False):
      print '\tRSS\t= {:.5e}\n\tBIC\t= {:.5e}\n\tAIC\t= {:.5e}\n\tstd\t= {:.5e}\n\tlikelihood\t= {:.5e}\n'.format(RSS, BIC, AIC, residualsSTD, likelihood)
    
    if ((np.isnan(AIC)) or (np.isnan(BIC))) and not(doSILENT):
      print bcolors.WARNING + '\tWARNING WARNING\n\tYou have NaN values for either AIC or BIC. Printing out all diagnostic values below...\n\tRSS\t= {:.5e}\n\tBIC\t= {:.5e}\n\tAIC\t= {:.5e}\n\tstd\t= {:.5e}\n\tlikelihood\t= {:.5e}\n'.format(RSS, BIC, AIC, residualsSTD, likelihood) + bcolors.ENDC
      
    return AIC, BIC, likelihood