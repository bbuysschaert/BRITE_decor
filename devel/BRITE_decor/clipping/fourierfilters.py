# -*- coding: utf-8 -*-
"""
Routines to perform filtering of a given flux array in Fourier space.
    
Last update 17 March 2016

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
def highPASSfilter(time, flux, cutOffFrequency, **kwargs):
  """
  Passes the lightcurve through a high-pass filter. This leaves only short-time variations in the lightcurve, while reducing any long-term effects.
  
  NOTE: An alternative way is to use a window averaging filter in the time domain. This behaves as a low pass filter. Dividing the low pass from the signal yields a high pass
  
  Returns: filtered flux measurements (centered around zero intensity units)
  
  @param time: time measurements [d]
  @type time: numpy array of length N
  @param flux: flux measurements [adu]
  @type flux: numpy array of length N
  @param cutOffFrequency: frequency limit for the highpass filter (needs inverse units of @param time)
  @type cutOffFrequency: integer
  
  @return: filteredFlux
  @rtype: numpy array of length N
  
  
  WARNING
  It is expected that your flux is in relative units. Nevertheless, the result from the filtered photometry results in photometry centered around zero.
  WARNING
  """
  """
  Equidistancy check
  ------------------
  """
  #The check is roughly in the order of seconds (if the time string is in days)
  meanTimeSpacing = np.mean(time[1:] - time[:-1])
  medianTimeSpacing = np.median(time[1:] - time[:-1])
  
  if np.round(meanTimeSpacing,4) == np.round(medianTimeSpacing,4):
    equidistant=True
  else:
    equidistant=False
  
  if equidistant == False:
    intFlux = scInterp.interp1d(time,flux) #This is a function
    equidistantTime = np.arange(time[0], time[-1]+np.round(medianTimeSpacing,5), np.round(medianTimeSpacing,5)) #Might be better, since it increases the length of the array
    equidistantTime[-1] = time[-1] #NOTE: this changes the equidistancy for the last pixels, yet it is necessary for a good backwards interpolation
    equidistantFlux = intFlux(equidistantTime)
  else:
    equidistantTime = time
    equidistantFlux = flux    
  arrayLength = len(equidistantTime)
  """
  Zero padding in order to have an even entry in the equidistant arrays
  -----------------
  """    
  if arrayLength % 2 == 0:
    paddingTime = np.arange(max(equidistantTime)+np.round(medianTimeSpacing,5), max(equidistantTime)+(arrayLength+1) * np.round(medianTimeSpacing,5), np.round(medianTimeSpacing,5))
  else:
    paddingTime = np.arange(max(equidistantTime)+np.round(medianTimeSpacing,5), max(equidistantTime)+(arrayLength) * np.round(medianTimeSpacing,5), np.round(medianTimeSpacing,5))
  equidistantTime = np.append(equidistantTime, paddingTime)
  equidistantFlux = np.append(equidistantFlux-np.mean(equidistantFlux),np.linspace(0,0,num=arrayLength))
  """
  Fourier transform
  -----------------
  """
  fourTransf = np.fft.rfft(equidistantFlux)
  fourierFreq = np.fft.fftfreq(len(equidistantTime),d=np.round(medianTimeSpacing,5)) #Needs implementation of the rfftfreq, np version 1.9
  
  maxIndexfourierFreq = np.int(len(fourTransf)-1) #Needed, since fftfreq gives you double the amount of elements (negative frequencies)
  idxToClip = np.where(fourierFreq[:maxIndexfourierFreq] < cutOffFrequency)[0]
  fourTransfFilter = np.copy(fourTransf)
  fourTransfFilter[idxToClip] = 0.  
  
  """
  Inverse fourier transform
  -------------------------
  """
  equidistantFluxFilter = np.fft.irfft(fourTransfFilter)

  if equidistant == False:
    intFluxFilter = scInterp.interp1d(equidistantTime,equidistantFluxFilter) #This is a function
    filteredFlux = intFluxFilter(time)
  else:
    filteredFlux = equidistantFluxFilter[:arrayLength+1]
  return filteredFlux

def lowPASSfilter(time, flux, cutOffFrequency, **kwargs):
  """
  Passes the lightcurve through a low-pass filter.
  This leaves only long-time variations in the lightcurve, while reducing any short-term effects.
  
  NOTE: An alternative way is to use a window averaging filter in the time domain. This behaves directly as a as a low pass filter. Dividing the low pass from the signal yields a high pass
  
  Returns: filtered flux measurements (centered around zero intensity units)
  
  @param time: time measurements [d]
  @type time: numpy array of length N
  @param flux: flux measurements [adu]
  @type flux: numpy array of length N
  @param cutOffFrequency: frequency limit for the lowpass filter (needs inverse units of @param time)
  @type cutOffFrequency: integer
  
  @return: filteredFlux
  @rtype: numpy array of length N
  
  
  WARNING
  It is expected that your flux is in relative units. Nevertheless, the result from the filtered photometry results in photometry centered around zero.
  WARNING
  """
  """
  Equidistancy check
  ------------------
  """
  #The check is roughly in the order of seconds (if the time string is in days)
  meanTimeSpacing = np.mean(time[1:] - time[:-1])
  medianTimeSpacing = np.median(time[1:] - time[:-1])
  
  if np.round(meanTimeSpacing,4) == np.round(medianTimeSpacing,4):
    equidistant=True
  else:
    equidistant=False
  
  if equidistant == False:
    intFlux = scInterp.interp1d(time,flux) #This is a function
    equidistantTime = np.arange(time[0], time[-1]+np.round(medianTimeSpacing,5), np.round(medianTimeSpacing,5)) #Might be better, since it increases the length of the array
    equidistantTime[-1] = time[-1] #NOTE: this changes the equidistancy for the last pixels, yet it is necessary for a good backwards interpolation
    equidistantFlux = intFlux(equidistantTime)
  else:
    equidistantTime = time
    equidistantFlux = flux    
  arrayLength = len(equidistantTime)
  """
  Zero padding in order to have an even entry in the equidistant arrays
  -----------------
  """
  if arrayLength % 2 == 0: #SOLVED an issue at some point, now it doesnt ...
    paddingTime = np.arange(max(equidistantTime)+np.round(medianTimeSpacing,5), max(equidistantTime)+(arrayLength) * np.round(medianTimeSpacing,5), np.round(medianTimeSpacing,5))
  else:
    paddingTime = np.arange(max(equidistantTime)+np.round(medianTimeSpacing,5), max(equidistantTime)+(arrayLength) * np.round(medianTimeSpacing,5), np.round(medianTimeSpacing,5))
  equidistantTime = np.append(equidistantTime, paddingTime)
  equidistantFlux = np.append(equidistantFlux-np.mean(equidistantFlux),np.linspace(0,0,num=arrayLength))
  print len(equidistantTime), len(equidistantFlux), len(paddingTime), arrayLength
  """
  Fourier transform
  -----------------
  """
  fourTransf = np.fft.rfft(equidistantFlux)
  fourierFreq = np.fft.fftfreq(len(equidistantTime),d=np.round(medianTimeSpacing,5)) #Needs implementation of the rfftfreq, np version 1.9
  
  maxIndexfourierFreq = np.int(len(fourTransf)-1) #Needed, since fftfreq gives you double the amount of elements (negative frequencies)
  idxToClip = np.where(fourierFreq[:maxIndexfourierFreq] > cutOffFrequency)[0] #NOTE Changed the sign here compared to the highPASSFIlter
  fourTransfFilter = np.copy(fourTransf)
  fourTransfFilter[idxToClip] = 0.  
  
  """
  Inverse fourier transform
  -------------------------
  """
  equidistantFluxFilter = np.fft.irfft(fourTransfFilter)

  if equidistant == False:
    intFluxFilter = scInterp.interp1d(equidistantTime,equidistantFluxFilter) #This is a function
    filteredFlux = intFluxFilter(time)
  else:
    filteredFlux = equidistantFluxFilter[:arrayLength+1]
  
  return filteredFlux

def multipleBANDPASSfilter(time, flux, centerOfBands, widthOfBands=None, invertBands=False, setNoise=False, **kwargs):
  """    
  Passes the lightcurve through a number of bandpass filters. The filtering is strictly done in Fourier space and matches an ideal filter.
  The only 'non-ideal' behaviour is made by the setting setNoise. When this setting is active, you don't set the fourier transform to 0+0j but to the noise in both the real and imaginary part.
  
  NOTE: the specified width for the frequency bands is the same for *all* bands
  
  Returns: filtered flux measurements (centered around zero intensity units)
  
  @param time: time measurements [d]
  @type time: numpy array of length N
  @param flux: flux measurements [adu]
  @type flux: numpy array of length N
  @param centerOfBands: central frequency of the filterbands (needs inverse units of @param time)
  @type centerOfBands: numpy array (can have length = 1)
  @param kwargs: optional arguments
  
  @return: filteredFlux
  @rtype: numpy array of length N
  
  @kwargs BANDSwidth: width of the filterbands - Default is 10*frequency resolution. [float]
  @kwargs BANDSonly: only take the flux which passes through the bands, instead of blocking it - Default is False [Boolean]
  @kwargs BANDSnoise: set the flux to an average noise level (in complex space) instead of putting int to 0+0j - Default is False [Boolean]
  """
  # Reading the kwargs
  widthOfBands = kwargs.get('BANDSwidth', 10 * 1./(time[-1]-time[0])) #[inverse units of @param time)]
  invertBands = kwargs.get('BANDSonly', False) #[Boolean]
  setNoise = kwargs.get('BANDSnoise', False) #[Boolean]
  
  # Trying to undersand the number of bands you are specifying
  try:
    numberOfBands = len(centerOfBands)
  except:
    if (isinstance(centerOfBands,int)) or (isinstance(centerOfBands,float)):
      numberOfBands = 1
    else:
      print "Data type for center of filterband(s) not understood"
      exit()
  """
  Equidistancy check
  ------------------
  """
  #The check is roughly in the order of seconds (if the time string is in days)
  meanTimeSpacing = np.mean(time[1:] - time[:-1])
  medianTimeSpacing = np.median(time[1:] - time[:-1])
  
  if np.round(meanTimeSpacing,4) == np.round(medianTimeSpacing,4):
    equidistant=True
  else:
    equidistant=False
  
  if equidistant == False:
    intFlux = scInterp.interp1d(time,flux) #This is a function
    equidistantTime = np.arange(time[0], time[-1]+np.round(medianTimeSpacing,5), np.round(medianTimeSpacing,5)) #Might be better, since it increases the length of the array
    equidistantTime[-1] = time[-1] #NOTE: this changes the equidistancy for the last pixels, yet it is necessary for a good backwards interpolation
    equidistantFlux = intFlux(equidistantTime)
  else:
    equidistantTime = time
    equidistantFlux = flux    
  arrayLength = len(equidistantTime)
  """
  Zero padding in order to have an even entry in the equidistant arrays
  -----------------
  """    
  if arrayLength % 2 == 0:
    paddingTime = np.arange(max(equidistantTime)+np.round(medianTimeSpacing,5), max(equidistantTime)+(arrayLength+1) * np.round(medianTimeSpacing,5), np.round(medianTimeSpacing,5))
  else:
    paddingTime = np.arange(max(equidistantTime)+np.round(medianTimeSpacing,5), max(equidistantTime)+(arrayLength) * np.round(medianTimeSpacing,5), np.round(medianTimeSpacing,5))
  equidistantTime = np.append(equidistantTime, paddingTime)
  equidistantFlux = np.append(equidistantFlux-np.mean(equidistantFlux),np.linspace(0,0,num=arrayLength))
  #print len(equidistantTime), len(equidistantFlux), len(paddingTime)
  """
  Fourier transform
  -----------------
  """
  fourTransf = np.fft.rfft(equidistantFlux)
  fourierFreq = np.fft.fftfreq(len(equidistantTime),d=np.round(medianTimeSpacing,5)) #Needs implementation of the rfftfreq, np version 1.9
  
  maxIndexfourierFreq = np.int(len(fourTransf)-1) #Needed, since fftfreq gives you double the amount of elements (negative frequencies)
  """
  Filtering
  ---------
  """
  idxToClipMax, idxToClipMin = [], []
  fourTransfFilter = np.copy(fourTransf)

  #Only one bandfilter
  if numberOfBands == 1:
    idxToClipMax.append(np.int64(np.where(fourierFreq[:maxIndexfourierFreq] < centerOfBands + widthOfBands/2.)[0][-1]))
    idxToClipMin.append(np.int64(np.where(fourierFreq[:maxIndexfourierFreq] > centerOfBands - widthOfBands/2.)[0][0]))
    
    if invertBands:
      if setNoise:
	fourTransfFilter[:idxToClipMin[0]] = np.median(fourTransf)
	fourTransfFilter[idxToClipMax[0]:] = np.median(fourTransf)
      else:
	fourTransfFilter[:idxToClipMin[0]] = 0
	fourTransfFilter[idxToClipMax[0]:] = 0
    else:
      if setNoise:
	fourTransfFilter[idxToClipMin[0]:idxToClipMax[0]+1] = np.median(fourTransf)
      else:
	fourTransfFilter[idxToClipMin[0]:idxToClipMax[0]+1] = 0
  #Multiple bandfilters
  else:
    for nn in range(numberOfBands):
      idxToClipMax.append(np.int64(np.where(fourierFreq[:maxIndexfourierFreq] < centerOfBands[nn] + widthOfBands/2.)[0][-1]))
      idxToClipMin.append(np.int64(np.where(fourierFreq[:maxIndexfourierFreq] > centerOfBands[nn] - widthOfBands/2.)[0][0]))
      
    if invertBands:
      for nn in range(numberOfBands):
	if setNoise:
	  if nn == 0:
	    fourTransfFilter[:idxToClipMin[nn]] = np.median(fourTransf)
	  elif nn == numberOfBands-1:
	    fourTransfFilter[idxToClipMax[nn]:] = np.median(fourTransf)
	    fourTransfFilter[idxToClipMax[nn-1]:idxToClipMin[nn]] = np.median(fourTransf)
	  else:
	    fourTransfFilter[idxToClipMax[nn-1]:idxToClipMin[nn]] = np.median(fourTransf)
	else:
	  if nn == 0:
	    fourTransfFilter[:idxToClipMin[nn]] = 0
	  elif nn == numberOfBands-1:
	    fourTransfFilter[idxToClipMax[nn]:] = 0
	    fourTransfFilter[idxToClipMax[nn-1]:idxToClipMin[nn]] = 0
	  else:
	    fourTransfFilter[idxToClipMax[nn-1]:idxToClipMin[nn]] = 0
    else:
      for nn in range(numberOfBands):
	if setNoise:
	  fourTransfFilter[idxToClipMin[nn]:idxToClipMax[nn]+1] = np.median(fourTransf)
	else:
	  fourTransfFilter[idxToClipMin[nn]:idxToClipMax[nn]+1] = 0
  """
  Inverse fourier transform
  -------------------------
  """
  equidistantFluxFilter = np.fft.irfft(fourTransfFilter)

  if equidistant == False:
    intFluxFilter = scInterp.interp1d(equidistantTime,equidistantFluxFilter) #This is a function
    filteredFlux = intFluxFilter(time)
  else:
    filteredFlux = equidistantFluxFilter[:arrayLength+1]

  return filteredFlux
