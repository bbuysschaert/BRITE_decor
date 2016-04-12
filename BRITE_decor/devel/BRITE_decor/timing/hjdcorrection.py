# -*- coding: utf-8 -*-
"""
Routines to correct the timing of the BRITE observations
    
Last update 16 March 2016

@author: Bram Buysschaert
"""

#===============================================================================
# 				Packages
#===============================================================================
import numpy as np

import scipy.interpolate as scInterp

#===============================================================================
# 				Code
#===============================================================================
def determineHELLCORR(JD, correction, **kwargs):
    """
    Use the JD and the calculated correction to determine a representation for the heliocentric correction. You can apply this to any other JD within the current JD interval. We calculate this representation using a spline fit (scipy.interpolate.splrep). You can then evaluate it for a given time using scipy.interpolate.splev(time, TCK).s
    
    WARNING this might give minor issues when you re-apply the fit if your last time stamp is outside the current JD interval
    
    Returns: the TCK of the fit (see scipy.interpolate.splrep)
    
    @param JD: Julian Dates [d]
    @type JD: numpy array of length N
    @param correction: a pre-calculated heliocentric correction
    @type correction: numpy array of length N
    
    @return tckJD: the TCK of the heliocentric correction.
    @rtype: TCK tuple (see scipy.interpolate.splrep)
    
    @kwargs: SPLINEknotpointsSPACING: spacing for knotpoints during the spline representation.- default is 5.0 [d]
    """
    SPLINEknotpointsSPACING = kwargs.get('SPLINEknotpointsSPACING', 5.0) #[d]
    
    
    knotpointsJD = np.arange(JD[0], JD[-1], SPLINEknotpointsSPACING)
    tckJD = scInterp.splrep(JD, correction, t=knotpointsJD[1:], k=3)
    
    return tckJD
  
def determineMIDobsTIME(JD, exposureTIME, stackTIME, numberSTACKS, **kwargs):
    """
    At present, all times given in the BRITE setup files, are beginning of the observation times. In case the exposure time is different for different setup files, you cannot just co-add two observations. Therefore, we provide a routine to determine the mid-exposure time.
    
    NOTE You have to apply the correction on the onboard JD, not the HJD!
    
    Return: the mid-exposure JDs
    
    @param JD: Julian Dates at the *beginning* of the observations [d]
    @type JD: numpy array of length N
    @param exposureTIME: exposure time for the observations [s]
    @type exposureTIME: numpy.float
    @param stackTIME: time it takes to stack an observation [s] - most often this will be 14s
    @type stackTIME: numpy.float
    @param numberSTACKS: number of stacked observations one observation represents [] - 1, 3, or 5
    @type numberSTACKS: numpy.float
    
    @return JD_MIDobs: Julian Dates at the *middle* of the observations [d]
    @rtype: numpy array of length N
    """
    
    # Conversion from seconds to day
    exposureTIME = exposureTIME / (24. * 3600.) #d
    stackTIME = stackTIME / (24. * 3600.) #d
    
    JD_MIDobs = JD + ((numberSTACKS * exposureTIME) + (numberSTACKS - 1.) * stackTIME) / 2.
    
    return JD_MIDobs