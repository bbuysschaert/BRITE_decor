# -*- coding: utf-8 -*-
"""
Routines which allow you to perform binning of the BRITE timeseries, using the longterm temperature variations. In addition, it is possibile to include the gapsize between datapoints for binning.
    
Last update 16 March 2016

@author: Bram Buysschaert
"""

#===============================================================================
# 				Packages
#===============================================================================
import numpy as np

import scipy.interpolate as scInterp

import statsmodels.api as sm
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
def convertSLOPEtoSIGN(slope, **kwargs):
    """
    A function that converts the value of the slope to a sign function (e.g. -1 or +1). We assume that a slope of zero has a positive sign and that we never have any NaN's
    
    Returns:
    The sign of the slope.
       
    @param slope: slope measurement
    @type time: numpy float
    
    @return: sign
    @rtype: numpy float
    
    @kwargs: SIGNinvert: invert the sign (might be useful) - default is False [Boolean]
    """
    invertSIGN = kwargs.get('SIGNinvert', False) #[Boolean]
    
    if not(invertSIGN):
      if slope > 0:
	return +1
      elif slope < 0:
	return -1
      elif slope == 0:
	print 'Congratulations, you have a slope with exactly zero!'
	return +1
    if invertSIGN:
      if slope < 0:
	return +1
      elif slope > 0:
	return -1
      elif slope == 0:
	print 'Congratulations, you have a slope with exactly zero!'
	return +1
    
def findENDofORBIT(time, index, **kwargs):
    """
    Routine to look for the nearest end of an orbit around a given index. 
    
    Returns:
    the index of the nearest end of an orbit
    
    @param time: time measurements [d]
    @type time: numpy array
    @param index: index within time array
    @type index: numpy integer
    
    @return: indexENDofORBIT
    @rtype: numpy integer
    
    @kwargs: IDXmaxSEARCH: maximum number of indexes you want to investigate to look for the end index of an orbit passage - default is 200 [idx]
    @kwargs: IDXsigma: integer indicating how large the difference has to be compared to the median time difference - default is 25 [integer]
    @kwargs: IDXforce: force to take a given index - default is False
    @kwargs: IDXtake: indicate which index is forced; needs to be used with forceIDX - default is indexUP.
    """
    # Reading the kwargs
    # ------------------
    maxIDXsearch = np.int(kwargs.get('IDXmaxSEARCH', 200)) #[idx]
    sigma = kwargs.get('IDXsigma', 25) # [integer]
    forceIDX = kwargs.get('IDXforce', False) #[Boolean]
    if forceIDX:
      takeIDX = kwargs.pop('IDXtake')
      if not takeIDX in ['IDXup', 'IDXdown']:
	raise ValueError, 'Please specify "IDXtake" while using "IDXforce", possible options are "IDXup" or "IDXdown"'
    
    # Performing the calculations
    # ---------------------------
    medianTIMEdiff = np.median(time[1:]-time[:-1])
    
    for ii in range(index, index+maxIDXsearch, 1):
      if time[ii+1] - time[ii] >= sigma*medianTIMEdiff:
	indexUP = ii  
	break
    for ii in range(index, index-maxIDXsearch, -1):
      if time[ii] - time[ii-1] >= sigma*medianTIMEdiff:
	indexDOWN = ii-1
	break
    
    # Check which one is the closest
    # ------------------------------
    # NOTE you want to check the time difference, not necessarily the index difference
    if not(forceIDX):
      if np.abs(time[index] - time[indexUP]) > np.abs(time[index] - time[indexDOWN]):
	return indexDOWN
      elif np.abs(time[index] - time[indexUP]) < np.abs(time[index] - time[indexDOWN]):
	return indexUP
      else:
	print bcolors.WARNING + '\t\tBoth time difference are equally large, taking the upper time difference.\n\t\tIndex \t\t= {:6.0f}\n\t\tIndexUP \t= {:6.0f}\t\ttimeUP \t= {:1.10e}\n\t\tIndexDOWN \t= {:6.0f}\t\ttimeDOWN \t= {:1.10e}'.format(index, indexUP, time[indexUP]-time[index], indexDOWN, time[indexDOWN]-time[index]) + bcolors.ENDC
	return indexUP
      
    elif takeIDX == 'IDXup':
      return indexUP
    elif takeIDX == 'IDXdown':
      return indexDOWN

def openTIMEwithTEMPERATURE(time, temperature, **kwargs):
    """
    Routine which will open the time base into smaller parts of the lightcurve, where each subset has a similar long term temperature. At the moment, up to two different techniques are used to determine the different timebins. i) the long term temperature difference and ii) the size of gaps within the data (if turned on; default is off).
    
    Diagnostics are provided on why a given index is considered as an endpoint for a time bin.
    
    NOTE
    At the moment we use a lowess filter to smoothen the temperature variations, determining the long term temperature variations.
    
    NOTE
    No minimum length for the time bins are specified. Thus, it is theoretically possible that one time bin has the length of one orbit. However, we provide a way to check if the length of the various time bins is too short compared to the specified minimum length. Yet, this is only a printing output and you should redefine your input parameters / start scripting if you want to change this.
    
    NOTE
    The slope of the long term variations is also considered and calculated within a window of slopeWINDOW, since nothing specifies that this is a monotoneous function. Whenever we pass a local extrema point, we change our temperature criterion from DeltaTEMP >= maxSPECIFIED to DeltaTEMP = 0. WARNING that slopeWINDOW might be able to significantly change your calculations. In case of onboard stacked data, this value should be significantly lower.
    
    Returns:
    list of indexes which should mark the end (include!) of the given time subset, a list of reasons why these indexes were selected, the long term temperature variations
    
    @param time: time measurements [d]
    @type time: numpy array of length N
    @param temperature: temperature measurements [deg]
    @type temperature: numpy array of length N
    
    @return: timeINDEXES
    @rtype: numpy array (dtype='int32') of length K
    @return: timeINDEXES_reason
    @rtype: list (dtype=string) of length K
    @return: temperatureLONG
    @rtype: numpy array of length N
    
    @kwargs: TEMPcrit: temperature difference between subsequent subset - default is 2.5 [deg]
    @kwargs: LOWESSfrac: fraction of the data you wish to use for the lowess filter - default is 0.2 (should be in range ~0.15 and ~0.35)
    @kwargs: GAPSinclude: stop your subsets when you reach a gap in the timing - default is False
    @kwargs: GAPScrit: size the gap has to be to stop the subset; needs to be used with stopGAPS - default is 0.3 [d]
    @kwargs: WINDOWslope: size in which you want to calculate the slope of the long term temperature variations - default is 100 [idx] (should be at least the approximate size of one orbit passage)
    @kwargs: BINminSIZE: the minimum timelength you would expect for a given temperature bin - default is 1.0 [d]
    
    TODO:
    -choose on how to handle extrema in your temperature (LATER)
    """
    
    dTEMP = kwargs.get('TEMPcrit', 2.5) # [deg]
    lowessFRAC = kwargs.get('LOWESSfrac', 0.2) # []
    stopGAPS = kwargs.get('GAPSinclude', False) # Boolean
    if stopGAPS:
      gapSIZE = kwargs.get('GAPScrit', 0.3) # [d]
    slopeWINDOW = np.int(kwargs.get('WINDOWslope', 100)) # [idx]
    binSIZEmin = kwargs.get('BINminSIZE', 1.0) # [d]
    
    #Determine the median amount of datapoints per day
    MEDIANdpointsDAY = 1./np.median(time[1:] - time[:-1])
    
    if MEDIANdpointsDAY < 3500: # You might have stacked data and need to adapt slopeWINDOW accordingly. This is only valid for BRITE data.
      print bcolors.WARNING + '\tWARNING\n\tYou might have onboard stacked data. Consider using a lower slopeWINDOW instead of {:3.0f}'.format(slopeWINDOW) + bcolors.ENDC
    
    """
        Determining the long term variations
     --------------------------------------------
    """
    # The long term temperature variations are calculated using a smoothing function. Since the time sampling between the different datapoints is not constant, you cannot use a convolution method. Therefore, a lowess (local linear regression fitting) filter is used, where the fraction of data considered is rather low (see lowessFRAC). If your setup is rather long, you might have to wait up to a few minutes for this step.
    
    print '\tPerforming lowess filtering ...'
    lowessLONG = sm.nonparametric.lowess(temperature, time, frac=lowessFRAC, delta=0.1)
    timeLONG, temperatureLONG = lowessLONG[:,0], lowessLONG[:,1]
    print bcolors.OKBLUE + '\t... all done.' + bcolors.ENDC
    if len(timeLONG) != len(time):
      print bcolors.FAIL + '\tERROR: The length of the longterm temperature array is different that the original array' + bcolors.ENDC
    
    temperatureMAX, temperatureMIN = np.max(temperatureLONG), np.min(temperatureLONG)
    # Check if the above are both at the beginning / end of the dataset. If not, inform the user.
    if (np.where(temperatureMAX==temperatureLONG)[0][0] != 0) and (np.where(temperatureMAX==temperatureLONG)[0][0] != len(temperatureLONG)-1):
      print bcolors.WARNING + '\tWARNING\n\tYour global *maximum* temperature is not at the beginning / end of the dataset. You will see this reflected in the subset bins.' + bcolors.ENDC
    if (np.where(temperatureMIN==temperatureLONG)[0][0] != 0) and (np.where(temperatureMIN==temperatureLONG)[0][0] != len(temperatureLONG)-1):
      print bcolors.WARNING + '\tWARNING\n\tYour global *minimum* temperature is not at the beginning / end of the dataset. You will see this reflected in the subset bins.' + bcolors.ENDC

    """
        Determining the temperature bins
     ---------------------------------------
    """
    # Defining the arrays for the indexes to define the different temperature bins
    timeINDEXES, timeINDEXES_reason = np.array([], dtype='int32'), []
    
    # While loop seems more appropriate, since it provides you with some control on when to skip elements. Although, you do need to define tt. Thus, initialising.
    # NOTE: not possible to determine the slope below this half (approximately one orbit for non-stacked data). Fully passed to the subset bin.
    # We use the longterm temperature behaviour to calculate the slope (linear flit within the window)
    temperatureINIT, indexINIT, slopeINIT = temperatureLONG[0], 0, np.polyfit(timeLONG[0:0+slopeWINDOW], temperatureLONG[0:0+slopeWINDOW], 1)[0]
    signINIT = convertSLOPEtoSIGN(slopeINIT)
    endBINtemperature = False
    tt = slopeWINDOW/2 + 1
    
    print '\n\tInitialising complete...'
    if not(stopGAPS):
      print '\tRemember, you have switched off to account for gaps while binning'
    
    while tt  < len(timeLONG):      
      # Not possible to determine the slope below this half (approximately one orbit). Fully passed to the subset bin.
      if tt >= len(timeLONG)-1 - slopeWINDOW/2 : # Last elements of your array. Stop here!
	timeINDEXES = np.append(timeINDEXES, len(timeLONG)-1); timeINDEXES_reason.append('end')
	try:
	  print '\tEnd of bin: ' + bcolors.BOLD + 'END' + bcolors.ENDC + '\t\t\tindex {:6.0f} (at time {:4.3f} d) with a binlength of {:4.3f} d.'.format(tt, time[tt]-time[0], time[len(timeLONG)-1]-time[timeINDEXES[-2]+1])
	  if time[len(timeLONG)-1]-time[timeINDEXES[-2]+1] < binSIZEmin:
	    print bcolors.WARNING + '\t\t\tWARNING: Timelength of this bin is smaller than the minimum length specified, i.e. {:1.3}d'.format(binSIZEmin) + bcolors.ENDC
	except:
	  print '\tEnd of bin: ' + bcolors.BOLD + 'END' + bcolors.ENDC + '\t\t\tindex {:6.0f} (at time {:4.3f} d) with a binlength of {:4.3f} d.'.format(tt, time[tt]-time[0], time[len(timeLONG)-1]-time[0])
	
	print bcolors.OKBLUE + '\t... found all endpoints for timebins.\n' + bcolors.ENDC
	break
      
      else: #Normal usage
	slopeTT = np.polyfit(timeLONG[tt-slopeWINDOW/2:tt+slopeWINDOW/2], temperatureLONG[tt-slopeWINDOW/2:tt+slopeWINDOW/2], 1)[0]
	signTT = convertSLOPEtoSIGN(slopeTT)
	
	# Deducing how the slope of your current datapoint compares to the first entry of your bin.
	if signINIT == +1 and signTT == +1 and temperatureLONG[tt] - temperatureINIT >= dTEMP:
	  endBINtemperature = True
	elif signINIT == -1 and signTT == -1 and temperatureINIT - temperatureLONG[tt] >= dTEMP:  
	  endBINtemperature = True
        elif signINIT == +1 and signTT == -1 and temperatureINIT >= temperatureLONG[tt]:
	  endBINtemperature = True
        elif signINIT == -1 and signTT == +1 and temperatureLONG[tt] >= temperatureINIT:  
	  endBINtemperature = True
        
        if endBINtemperature == True:        
	  try:
	    print '\tEnd of bin: ' + bcolors.BOLD + 'Temperature' + bcolors.ENDC + '\t\tindex {:6.0f} (at time {:4.3f} d) with a binlength of {:4.3f} d.'.format(tt, time[tt]-time[0], time[tt]-time[timeINDEXES[-1]+1])
	    if time[tt]-time[timeINDEXES[-1]+1] < binSIZEmin:
	      print bcolors.WARNING + '\t\t\tWARNING: Timelength of this bin is smaller than the minimum length specified, i.e. {:1.3}d'.format(binSIZEmin) + bcolors.ENDC
	  except:
	    print '\tEnd of bin: ' + bcolors.BOLD + 'Temperature' + bcolors.ENDC + '\t\tindex {:6.0f} (at time {:4.3f} d) with a binlength of {:4.3f} d.'.format(tt, time[tt]-time[0], time[tt]-time[0]) #This should only happen with the first bin
	    if time[tt]-time[0] < binSIZEmin:
	      print bcolors.WARNING + '\t\t\tWARNING: Timelength of this bin is smaller than the minimum length specified, i.e. {:1.3}d'.format(binSIZEmin) + bcolors.ENDC
	  print '\t\tTdiff = {:+2.2f}deg;\tdTEMP limit is {:+2.2f}deg.'.format(temperatureINIT - temperatureLONG[tt], dTEMP)
	  # Now take a more detailed look into which index to actually consider your final index.
	  if signTT != 0:
	    indexENDofORBIT = findENDofORBIT(time, tt, **kwargs)
	  else:
	    indexENDofORBIT = findENDofORBIT(time, tt, IDXforce=True, IDXtake='IDXup', **kwargs)
	  timeINDEXES = np.append(timeINDEXES, indexENDofORBIT); timeINDEXES_reason.append('temperature')
	  
	  # Make the new initial conditions for the next bin.
	  tt = indexENDofORBIT + 1
	  temperatureINIT, indexINIT = temperatureLONG[tt], tt
	  slopeINIT = np.polyfit(timeLONG[tt:tt+slopeWINDOW], temperatureLONG[tt:tt+slopeWINDOW], 1)[0]
	  signINIT = convertSLOPEtoSIGN(slopeINIT)
	  tt += 1; endBINtemperature = False
	  
	# Check if you did not pass a gap with a given gapSIZE. If so, take the previous index as your end of the bin.
	# NOTE this index should (normally) always be the end of an orbit, so no need to look into more detail.
	elif stopGAPS and time[tt+1] - time[tt] > gapSIZE: 
	  try:
	    print '\tEnd of bin: ' + bcolors.BOLD + 'Gap' + bcolors.ENDC + '\t\t\tindex {:6.0f} (at time {:4.3f} d) with a binlength of {:4.3f} d.'.format(tt, time[tt]-time[0], time[tt]-time[timeINDEXES[-1]+1])
	    if time[tt]-time[timeINDEXES[-1]+1] < binSIZEmin:
	      print bcolors.WARNING + '\t\t\tWARNING: Timelength of this bin is smaller than the minimum length specified, i.e. {:1.3}d'.format(binSIZEmin) + bcolors.ENDC
	  except:
	    print '\tEnd of bin: ' + bcolors.BOLD + 'Gap' + bcolors.ENDC + '\t\t\tindex {:6.0f} (at time {:4.3f} d) with a binlength of {:4.3f} d.'.format(tt, time[tt]-time[0], time[tt]-time[0]) #This should only happen with the first bin
	    if time[tt]-time[0] < binSIZEmin:
	      print bcolors.WARNING + '\t\t\tWARNING: Timelength of this bin is smaller than the minimum length specified, i.e. {:1.3}d'.format(binSIZEmin) + bcolors.ENDC
	  print '\t\ttdiff = {:2.3f}d;\t\tgapSIZE limit is {:2.3f}d.'.format(time[tt+1] - time[tt], gapSIZE)  
	  timeINDEXES = np.append(timeINDEXES, tt); timeINDEXES_reason.append('gap')
	
	  # Make the new initial conditions for the next bin.
	  tt = tt + 1
	  temperatureINIT, indexINIT = temperatureLONG[tt], tt
	  slopeINIT = np.polyfit(timeLONG[tt:tt+slopeWINDOW], temperatureLONG[tt:tt+slopeWINDOW], 1)[0]
	  signINIT = convertSLOPEtoSIGN(slopeINIT)
	
	# Temperature difference is not large enough to end your bin nor did you find a gap. So, you continue looping.
	else:
          tt += 1
    
    return timeINDEXES, timeINDEXES_reason, temperatureLONG
