# -*- coding: utf-8 -*-
"""
Routines to study the data of BRITE observations per orbit
    
Last update 18 March 2016

@author: Bram Buysschaert
"""

#===============================================================================
# 				Packages
#===============================================================================
import numpy as np
#===============================================================================
# 				Code
#===============================================================================
def averageORBIT(time, flux, **kwargs):
    """
    Rebinning of the photometric measurements to bins of 1 orbit passage. The photometric data points are taken with almost a second cadence, which is too high.
    Yet, we can exploit this to our advantage. Rebinning allows us to calculate a mean, median and standard deviation (std) per bin.
    
    (still seems rather versatile and fast)
    
    WARNING     while loops	  WARNING
    WARNING beware of what you do WARNING
    
    Returns: a matrix with diagnostics about the time per orbital bin, a matrix with diagnostics about the flux per orbital bin, an array with the number of elements per orbital bin.
    
    @param time: time measurements [d]
    @type time: numpy array of length N
    @param flux: flux measurements [adu]
    @type flux: numpy array of length N    
    
    @return: timeBINmatrix: diagnostics about the time per orbital bin
    @rtype: numpy matrix of size 3xK
    @return: fluxBINmatrix:  diagnostics about the flux per orbital bin
    @rtype: numpy matrix of size 3xK
    @return: binNUMBERofELEMENTS: number of elements per orbital bin
    @rtype: numpy array of length K (dtype='int32')
    
    @kwargs: LIMIThard: impose a hard limit on the bins, if so specify the time difference in days - Default is not given, so the soft limit is the Default option [float]
    @kwargs: LIMITsoft: impose a soft limit on the bins, if so specity the multiple of the __median__ time difference in days - Default is 5 [float]
    """
    
    # Reading in the kwargs
    hard_limit = kwargs.get('LIMIThard') # When you did not specify LIMIThard, hard_limit == None.   [d]
    if hard_limit == None:
      soft_limit = kwargs.get('LIMITsoft', 5) #[float]
      time_difference = time[1:]-time[:-1]
      soft_limit = soft_limit * np.median(time_difference)
    
    # Setting up the differen lists, which will later be converted to a matrix.
    timeBINmean, timeBINmedian, timeBINstd = [], [], []
    fluxBINmean, fluxBINmedian, fluxBINstd = [], [], []
    binNUMBERofELEMENTS = []

    # Work with a copy. You might need the original at some later point.
    timeCOPY, fluxCOPY = np.copy(time), np.copy(flux)
        
    # Initialising the bins.
    timeBIN, fluxBIN = [], []
    while len(timeCOPY) > 0:
      if len(timeCOPY) > 2: # Normal usage.
	if hard_limit != None: #Hence, you specifed a hard limit
	  if timeCOPY[1]-timeCOPY[0] < hard_limit: # Difference is too small. Continue.
	    timeBIN.append(float(timeCOPY[0])); fluxBIN.append(float(fluxCOPY[0]))
	    #Remove first element to get out of the while loop
	    timeCOPY, fluxCOPY = np.delete(timeCOPY, 0), np.delete(fluxCOPY, 0)
	    
	  else: # Difference is too large. Stop and calculate the characteristics of the bin
	    #Do not forget the last item
	    timeBIN.append(float(timeCOPY[0])); fluxBIN.append(float(fluxCOPY[0]))
	    #Remove first element to get out of the while loop
	    timeCOPY, fluxCOPY = np.delete(timeCOPY, 0), np.delete(fluxCOPY, 0)
	    
	    # Time related storage
	    timeBINmean.append(float(np.mean(timeBIN))); timeBINmedian.append(float(np.median(timeBIN)));     timeBINstd.append(float(np.std(timeBIN)))
	    # Flux related storage
	    fluxBINmean.append(float(np.mean(fluxBIN))); fluxBINmedian.append(float(np.median(fluxBIN)));     fluxBINstd.append(float(np.std(fluxBIN)))
	    # Number of elementes / bin storage
	    binNUMBERofELEMENTS.append(int(len(timeBIN)))
	    
	    #Clearing the bin
	    timeBIN, fluxBIN = [], []
	else:
	  if timeCOPY[1]-timeCOPY[0] < soft_limit: # Difference is too small. Continue.
	    timeBIN.append(float(timeCOPY[0])); fluxBIN.append(float(fluxCOPY[0]))
	    #Remove first element to get out of the while loop
	    timeCOPY, fluxCOPY = np.delete(timeCOPY, 0), np.delete(fluxCOPY, 0)
	  else: # Difference is too large. Stop and calculate the characteristics of the bin
	    #Do not forget the last item
	    timeBIN.append(float(timeCOPY[0])); fluxBIN.append(float(fluxCOPY[0]))
	    #Remove first element to get out of the while loop
	    timeCOPY, fluxCOPY = np.delete(timeCOPY, 0), np.delete(fluxCOPY, 0)
	    
	    # Time related storage
	    timeBINmean.append(float(np.mean(timeBIN))); timeBINmedian.append(float(np.median(timeBIN)));     timeBINstd.append(float(np.std(timeBIN)))
	    # Flux related storage
	    fluxBINmean.append(float(np.mean(fluxBIN))); fluxBINmedian.append(float(np.median(fluxBIN)));     fluxBINstd.append(float(np.std(fluxBIN)))
	    # Number of elementes / bin storage
	    binNUMBERofELEMENTS.append(int(len(timeBIN)))
	    
	    #Clearing the bin
	    timeBIN, fluxBIN = [], []
	    
      else: # Less than three elements, indicate the end of your array and you put them all in the same bin. In case you did a proper clipping, you should not have any orbits with less than three elements.
	timeBIN.append(float(timeCOPY[0])); fluxBIN.append(float(fluxCOPY[0]))
	timeBIN.append(float(timeCOPY[1])); fluxBIN.append(float(fluxCOPY[1]))
	#Remove first element to get out of the while loop
	timeCOPY, fluxCOPY = np.delete(timeCOPY, [0,1]), np.delete(fluxCOPY, [0,1])
	
	# Time related storage
	timeBINmean.append(float(np.mean(timeBIN))); timeBINmedian.append(float(np.median(timeBIN)));     timeBINstd.append(float(np.std(timeBIN)))
	# Flux related storage
	fluxBINmean.append(float(np.mean(fluxBIN))); fluxBINmedian.append(float(np.median(fluxBIN)));     fluxBINstd.append(float(np.std(fluxBIN)))
	# Number of elementes / bin storage
	binNUMBERofELEMENTS.append(int(len(timeBIN)))
	
	#Clearing the bin (not needed, but could be used for checks if mandatory)
	timeBIN, fluxBIN = [], []
    
    # Convert the different lists to matrices.
    timeBINmatrix = np.array([timeBINmean, timeBINmedian, timeBINstd])
    fluxBINmatrix = np.array([fluxBINmean, fluxBINmedian, fluxBINstd])
    
    return timeBINmatrix, fluxBINmatrix, np.array(binNUMBERofELEMENTS, dtype='int32')

def phaseORBIT(time, Porbit, **kwargs):
    """
    Determine the orbital phase to a common epoch, which we fix throughout the analysis. This can be useful if you want to study data from the same satellite over different observing periods.
    
    DANGER: This is off course only valid if your satellite orbital period does not change over time.
    WARNING: it assumed that the Porbit is given in minutes
    
    Returns: the orbital phase, compared to the defined epoch.
    
    @param time: time measurements [d]
    @type time: numpy array of length N
    @param Porbit: orbital period of the satellite [min]NOTE
    @type Porbit: float
    
    @return: orbitPHASE: orbital phase of satellite compared to the specified epoch
    @rtype: numpy array of length N
    
    
    @kwargs: EPOCH: the specified epoch (date) you want to compare your observations to - Default is 2456293.5  [float]
    @kwargs: LIMITsoft: impose a soft limit on the bins, if so specity the multiple of the __median__ time difference in days - Default is 5 [float]
    
    The normal output will give you the number of the orbit (compared to the epoch) and the fraction (phase) of the orbit.
    The phase output will give you the phase of that particular orbit it is making.
    """
    # Reading in the kwargs
    epoch = kwargs.get('EPOCH', 2456293.5) #The standard epoch is the 1st of January 2013, midnight.
    
    # Converting the orbital period to days
    PorbitD = Porbit / (24. * 60.)
    
    return (time-epoch) % PorbitD / PorbitD

def timeORBIT(time, Porbit, **kwargs):
    """
    Determine the timing of the observations to a common epoch, which we fix throughout the analysis. This can be useful if you want to study data from the same satellite over different observing periods.
    
    DANGER: This is off course only valid if your satellite orbital period does not change over time.
    WARNING: it assumed that the Porbit is given in minutes
    
    Returns: integer = orbital number compared to the epoch, digits = fraction of the current orbit done
    
    @param time: time measurements [d]
    @type time: numpy array of length N
    @param Porbit: orbital period of the satellite [min]NOTE
    @type Porbit: float
    
    @return: orbitTIME: timing of the observations to the specified epoch
    @rtype: numpy array of length N
    
    
    @kwargs: EPOCH: the specified epoch (date) you want to compare your observations to - Default is 2456293.5  [float]
    @kwargs: LIMITsoft: impose a soft limit on the bins, if so specity the multiple of the __median__ time difference in days - Default is 5 [float]
    
    The normal output will give you the number of the orbit (compared to the epoch) and the fraction (phase) of the orbit.
    The phase output will give you the phase of that particular orbit it is making.
    """
    # Reading in the kwargs
    epoch = kwargs.get('EPOCH', 2456293.5) #The standard epoch is the 1st of January 2013, midnight.
    
    # Converting the orbital period to days
    PorbitD = Porbit / (24. * 60.)
    
    return (time-epoch)/PorbitD

def analyseNUMBERperORBIT(time, flux, numberSTACKS, Porbit, **kwargs):
    """
    Study the flux per orbit passage of the satellite. This is a variation on the averageORBIT routine. Moreover, this one also permits the user to mark orbits with not enough elements per bin as outliers.
    
    NOTE: Much more diagnostics are calculated in this routine, but are not outputted (yet).
    
    Returns: Diagnostics per orbital bin and the indexes of the datapoints, where not enough observations were taken during an orbit passage.
    
    @param time: time measurements [d]
    @type time: numpy array of length N
    @param flux: flux measurements [adu]
    @type flux: numpy array of length N    
    @param numberSTACKS: number of onboard stacked datapoints corresponding to one datapoint
    @type numberSTACKS: numpy array of length N    
    @param Porbit: orbital period of the satellite [min]NOTE
    @type Porbit: float
    
    @return: timeBINmean: mean time for the orbital passage
    @rtype: numpy array of length K
    @return: fluxBINmean: mean flux for the orbital passage
    @rtype: numpy array of length K
    @return: fluxBINstd: std on the flux for the orbital passage
    @rtype: numpy array of length K
    @return: binNUMBERofELEMENTS: number of observations for the orbital passage
    @rtype: numpy array of length K
    @return: idxBINoutliers: indexes of the the outliers
    @rtype: numpy array of length M
    
    @kwargs: ORBITminOBSERVATION: minimum number of observations one orbit passage has to have, to not be considered as an outlier - Default is 5 [float]
    @kwargs: ORBITphaseDIFF: difference in orbital phase, to be considered as two different satellite passages - Default is 0.4 [float]
    """
    # Reading in the kwargs
    minPerORBIT = kwargs.get('ORBITminOBSERVATION',5) # [Float]
    minSpacingORBIT = kwargs.get('ORBITphaseDIFF', 0.4) # [Float]
    
    # Getting the 'orbit time'. The epoch is passed through the kwargs
    orbitTIME = timeORBIT(time, Porbit, **kwargs)
    
    # Getting diagnostic values per orbit passage of the satellite
    # Passages with not enough datapoints are marked as outliers
    
    # Setting up the differen lists.
    timeBINmean, fluxBINmean, fluxBINstd, binNUMBERofELEMENTS = [], [], [], []
    idxBINoutliers = np.array([], dtype='int32')
    
    # Initialising the bins.
    timeBIN, fluxBIN, idxBIN, stackBIN = [], [], [], []

    # Loop over your time array.
    for tt in range(0,len(orbitTIME)-1,1):
      if np.abs(orbitTIME[tt+1] - orbitTIME[tt]) >= minSpacingORBIT: # End of your bin (=orbit)
	timeBIN.append(float(time[tt])); fluxBIN.append(float(flux[tt])); idxBIN.append(int(tt)); stackBIN.append(float(numberSTACKS[tt]))
	
	timeBINmean.append(float(np.mean(timeBIN))); fluxBINmean.append(float(np.mean(fluxBIN))); fluxBINstd.append(float(np.std(fluxBIN))); binNUMBERofELEMENTS.append(int(len(timeBIN)))
	# Checking if you did not have enough datapoints during this orbit (taking into account the onboard stacking)
	if len(timeBIN) * np.unique(stackBIN) <= minPerORBIT:
	  idxBINoutliers = np.append(idxBINoutliers, idxBIN)
	
	# Clearing the bins
	timeBIN, fluxBIN, idxBIN, stackBIN = [], [], [], []
	
      elif tt == len(orbitTIME)-2: # Last two elements of your array. Stop here.
	timeBIN.append(float(time[tt])); fluxBIN.append(float(flux[tt])); idxBIN.append(int(tt)); stackBIN.append(float(numberSTACKS[tt]))
	timeBIN.append(float(time[tt+1])); fluxBIN.append(float(flux[tt+1])); idxBIN.append(int(tt+1)); stackBIN.append(float(numberSTACKS[tt+1]))
	
	timeBINmean.append(float(np.mean(timeBIN))); fluxBINmean.append(float(np.mean(fluxBIN))); fluxBINstd.append(float(np.std(fluxBIN))); binNUMBERofELEMENTS.append(int(len(timeBIN)))
	# Checking if you did not have enough datapoints during this orbit (taking into account the onboard stacking)
	if len(timeBIN) * np.unique(stackBIN) <= minPerORBIT:
	  idxBINoutliers = np.append(idxBINoutliers, idxBIN)
	
	# Clearing the bin (not needed, but could be used for checks if mandatory)
	timeBIN, fluxBIN, idxBIN, stackBIN = [], [], [], []
	  
      else: # Normal usage
	timeBIN.append(float(time[tt])); fluxBIN.append(float(flux[tt])); idxBIN.append(int(tt)); stackBIN.append(float(numberSTACKS[tt]))
    
    
    timeBINmean, fluxBINmean, fluxBINstd, binNUMBERofELEMENTS = np.array(timeBINmean), np.array(fluxBINmean), np.array(fluxBINstd), np.array(binNUMBERofELEMENTS)
    
    
    return timeBINmean, fluxBINmean, fluxBINstd, binNUMBERofELEMENTS, idxBINoutliers