# -*- coding: utf-8 -*-
"""
Routines to load various files, which are either provided by the BRITE data reduction team or from your own analysis.
    
Last update 18 March 2016

@author: Bram Buysschaert
"""

#===============================================================================
# 				Packages
#===============================================================================
import numpy as np
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
def loadSETUP(fileIN, pathIN, **kwargs):
    """
    Routine to read the default setup files provided by the BRITE community. This routine does not only read data from the file, but also retrieves crucial information from the header of the file.
    
    Normally, it should work for both DR1 and DR2 output, although the user has to specify it when it is not DR2! (NOTE various DR2 reductions exist, but should all have the same output).
    
    Routine that loads the data files from BRITE per satellite. It can handle 'old' and 'new' data outputs, as seen during early October 2015
    
    Returns: all information from the setup file and from the header
    
    @param fileIN: name the setup file has
    @type fileIN: string
    @param pathIN: path to the location where the setup file is located
    @type pathIN: string
    
    @return HJD: heliocentric JD of the observations [d]
    @rtype HJD: numpy array of length N
    @return fluxRAW: raw flux measurements [adu]
    @rtype fluxRAW: numpy array of length N
    @return xPOS: CCD position measurements along the x position axis [pixel]
    @rtype xPOS: numpy array of length N
    @return yPOS: CCD position measurements along the x position axis [pixel]
    @rtype yPOS: numpy array of length N
    @return temperature: temperature measurements [deg]
    @rtype temperature: numpy array of length N
    @return JD: JD of the observations [d] != HJD
    @rtype JD: numpy array of length N
    @return qFLAG: quality flag for the aperture rendering [0, 1]
    @rtype qFLAG: numpy array of length N
    
    @return exposureTIME: exposure time for the observations [s]
    @rtype exposureTIME: numpy array of length N
    @return numberSTACKS: number of onboard stacked datapoints corresponding to this datapoint []
    @rtype numberSTACKS: numpy array of length N
    
    @return xRASTER: size of the CCD raster, along the x-axis, on which observations were taken [pixel]
    @rtype xRASTER: numpy float
    @return yRASTER: size of the CCD raster, along the x-axis, on which observations were taken [pixel]
    @rtype yRASTER: numpy float
    @return aperture: radius of the circular aperture, used to extract the photometry [pixel]
    @rtype aperture: numpy float
    
    @kwargs: DR: type of data reduction used to create the setup file - Default is DR2 [DR1, DR2]    
    """
    # Reading in the kwargs
    dataTYPE = kwargs.get('DR', 'DR2')
    if not dataTYPE in ['DR1', 'DR2']:
      raise ValueError, 'Please specify the "DR" properly as either "DR1" or "DR2".'
    # Checking if the last character of pathIN is an '/'
    if not(pathIN[-1] == '/'):
      pathIN += '/'
    # Checking if the suffix of the file is given
    if not fileIN[-4:] in ['.txt', '.dat']:
      fileIN += '.dat'  
    
    # Open the setup file
    if dataTYPE == 'DR1':
      HJD, fluxRAW, xPOS, yPOS, temperature, heliocentricCORRECTION = np.loadtxt(pathIN + fileIN, usecols=(0,1,2,3,4,5), unpack=True, comments='c')
      qFLAG = np.ones_like(HJD_set) #These were not calculated for DR1
      JD = HJD - heliocentricCORRECTION/(3600.*24.) #d
    if dataTYPE == 'DR2':  
      HJD, fluxRAW, xPOS, yPOS, temperature, JD, qFLAG = np.loadtxt(pathIN + fileIN, usecols=(0,1,2,3,4,5,6), unpack=True, comments='c')
      
    # Open the data file to read the header and get useful information out
    fileHEADER = open(pathIN + fileIN, 'r')
    for hh in range(100): # Should be enough to loop through the whole header
      lineHEADER = fileHEADER.readline()
      lineHEADER = lineHEADER.split(' ') # Split the string at spaces
      # It is the *second* element in the lineHEADER array that we are checking for the header keywords
      if lineHEADER[1] == 'ObsExpoT=':
	exposureTIME = np.float(lineHEADER[2]) / 1000. #s
      elif lineHEADER[1] == 'ObsStack=':
	numberSTACKS = np.float(lineHEADER[2])
      elif lineHEADER[1] == 'ROIxsiz':
	xRASTER = np.float(lineHEADER[3])
      elif lineHEADER[1] == 'ROIysiz':
	yRASTER = np.float(lineHEADER[3])  
      elif lineHEADER[1] == 'RedApert=':
	aperture = np.float(lineHEADER[2]) 
    fileHEADER.close() # Close the file to save memory
    
    return HJD, fluxRAW, xPOS, yPOS, temperature, JD, qFLAG, np.ones_like(HJD)*exposureTIME, np.ones_like(HJD)*numberSTACKS, xRASTER, yRASTER, aperture
  
def loadCLIPPED(fileIN, pathIN, **kwargs):
    """
    Routine to read the output of the clipping process, saved by saveCLIPPED.
    
    Returns: all information from the file.
    
    @param fileIN: name the setup file has
    @type fileIN: string
    @param pathIN: path to the location where the setup file is located
    @type pathIN: string
    
    @return time: heliocentric JD of the observations [d]
    @rtype time: numpy array of length N
    @return flux: flux measurements [adu]
    @rtype flux: numpy array of length N
    @return xPOS: CCD position measurements along the x position axis [pixel]
    @rtype xPOS: numpy array of length N
    @return yPOS: CCD position measurements along the x position axis [pixel]
    @rtype yPOS: numpy array of length N
    @return temperature: temperature measurements [deg]
    @rtype temperature: numpy array of length N
    
    @return exposureTIME: exposure time for the observations [s]
    @rtype exposureTIME: numpy array of length N
    @return numberSTACKS: number of onboard stacked datapoints corresponding to this datapoint []
    @rtype numberSTACKS: numpy array of length N
    """
    # Checking if the last character of pathIN is an '/'
    if not(pathIN[-1] == '/'):
      pathIN += '/'
    # Checking if the suffix of the file is given
    if not fileIN[-4:] in ['.txt', '.dat']:
      fileIN += '.dat'  
    
    time, flux, xPOS, yPOS, temperature, exposureTIME, numberSTACKS = np.loadtxt(pathIN + fileIN, unpack = True, comments='#')
    
    return time, flux, xPOS, yPOS, temperature, exposureTIME, numberSTACKS

def loadDETREND(fileIN, pathIN, **kwargs):
    """
    Routine to read the output of the detrending process, saved by saveDETREND.
    
    Returns: all information from the file.
    
    @param fileIN: name the setup file has
    @type fileIN: string
    @param pathIN: path to the location where the setup file is located
    @type pathIN: string
    
    @return time: heliocentric JD of the observations [d]
    @rtype time: numpy array of length N
    @return flux: flux measurements [adu]
    @rtype flux: numpy array of length N
    @return xPOS: CCD position measurements along the x position axis [pixel]
    @rtype xPOS: numpy array of length N
    @return yPOS: CCD position measurements along the x position axis [pixel]
    @rtype yPOS: numpy array of length N
    @return temperature: temperature measurements [deg]
    @rtype temperature: numpy array of length N
    
    @return exposureTIME: exposure time for the observations [s]
    @rtype exposureTIME: numpy array of length N
    @return numberSTACKS: number of onboard stacked datapoints corresponding to this datapoint []
    @rtype numberSTACKS: numpy array of length N
    
    @return correction: correction for the flux to detrend it [adu]
    @rtype correction: numpy array of length N
    """
    # Checking if the last character of pathIN is an '/'
    if not(pathIN[-1] == '/'):
      pathIN += '/'
    # Checking if the suffix of the file is given
    if not fileIN[-4:] in ['.txt', '.dat']:
      fileIN += '.dat'    
      
    time, flux, xPOS, yPOS, temperature, exposureTIME, numberSTACKS, correction = np.loadtxt(pathIN + fileIN, unpack = True)
    
    return time, flux, xPOS, yPOS, temperature, exposureTIME, numberSTACKS, correction
  
def loadPSFdetrend(fileIN, pathIN, **kwargs):
    """
    Routine to read the output of the PSF detrending process, saved by OUT_FLUX_detrendTEMPpsfFULL.
    
    Returns: all information from the file.
    
    @param fileIN: name the setup file has
    @type fileIN: string
    @param pathIN: path to the location where the setup file is located
    @type pathIN: string
    
    @return time: heliocentric JD of the observations [d]
    @rtype time: numpy array of length N
    @return flux: flux measurements [adu]
    @rtype flux: numpy array of length N
    @return fluxCORRECTED = flux measurements corrected for the temperature dependent PSF changes[adu]
    @rtype fluxCORRECTED: numpy array of length N    
    @return correction: correction for the flux to detrend it [adu]
    @rtype correction: numpy array of length N
    """
    # Checking if the last character of pathIN is an '/'
    if not(pathIN[-1] == '/'):
      pathIN += '/'
    # Checking if the suffix of the file is given
    if not fileIN[-4:] in ['.txt', '.dat']:
      fileIN += '.dat'    
      
    time, flux, fluxCORRECTED, correction = np.loadtxt(pathIN + fileIN, unpack = True)
    
    return time, flux, fluxCORRECTED, correction