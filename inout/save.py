# -*- coding: utf-8 -*-
"""
Routines to save output to various files.
    
Last update 18 March 2016

@author: Bram Buysschaert
"""

#===============================================================================
# 				Packages
#===============================================================================
import numpy as np

from time import strftime
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
def saveCLIPPED(fileOUT, pathOUT, time, flux, xPOS, yPOS, temperature, exposureTIME, numberSTACKS, **kwargs):
    """
    Routine to save the output of the clipping process.
    
    Returns: Nothing, but saves a file named fileOUT at the specified location pathOUT.
    
    @param fileOUT: name you wish to give to the output file (.txt NO need to be specified)
    @type fileOUT: string
    @param pathOUT: path to the location where you wish to save the output
    @type pathOUT: string
    
    @param time: heliocentric JD of the observations [d]
    @type time: numpy array of length N
    @param flux: flux measurements [adu]
    @type flux: numpy array of length N
    @param xPOS: CCD position measurements along the x position axis [pixel]
    @type xPOS: numpy array of length N
    @param yPOS: CCD position measurements along the x position axis [pixel]
    @type yPOS: numpy array of length N
    @param temperature: temperature measurements [deg]
    @type temperature: numpy array of length N
    
    @param exposureTIME: exposure time for the observations [s]
    @type exposureTIME: numpy array of length N
    @param numberSTACKS: number of onboard stacked datapoints corresponding to this datapoint []
    @type numberSTACKS: numpy array of length N
    """
    # Checking if the last character of pathOUT is an '/'
    if not(pathOUT[-1] == '/'):
      pathOUT += '/'
    # Checking if the suffix of the file is given
    if not fileOUT[-4:] in ['.txt', '.dat']:
      fileOUT += '.dat'  
      
    # Preparing the header of the output file
    headerSTRING = 'BRITE photometry, which was clipped for outliers on ' + strftime("%Y-%m-%d %H:%M:%s") + '.'
    headerSTRING +='\n----------------------------------------'
    headerSTRING +='\nColumn1: time measurements [d]'
    headerSTRING +='\nColumn2: flux [adu]'
    headerSTRING +='\nColumn3: CCD centroid x-position [pixel]'
    headerSTRING +='\nColumn4: CCD centroid y-position [pixel]'
    headerSTRING +='\nColumn5: CCD temperature [deg]'
    headerSTRING +='\nColumn6: exposure time of the observations [s]'
    headerSTRING +='\nColumn7: number of stacked observations corresponding to one datapoint []'
    headerSTRING +='\n----------------------------------------'
    
    # Constructing the matrix
    dtOUT = np.dtype([('time', np.float), ('flux', np.float), ('xPOS', np.float), ('yPOS', np.float), ('temperature', np.float), ('exposureTIME', np.float), ('numberSTACKS', np.float)])
    matrixOUT = np.zeros(len(time), dtype=dtOUT)
    matrixOUT['time'] = time; matrixOUT['flux'] = flux; matrixOUT['xPOS'] = xPOS; matrixOUT['yPOS'] = yPOS; matrixOUT['temperature'] = temperature; matrixOUT['exposureTIME'] = exposureTIME; matrixOUT['numberSTACKS'] = numberSTACKS
    
    # The actual saving using a numpy.savetxt    
    np.savetxt(pathOUT + fileOUT, matrixOUT, fmt=('%.12e %.7f %.4f %.4f %.4f %.2f %i'), delimiter=' ', header=headerSTRING, comments='#')

def saveDETREND(fileOUT, pathOUT, time, fluxRAW, xPOS, yPOS, temperature, exposureTIME, numberSTACKS, correction, **kwargs):
    """
    Routine to save the output of the detrending process.
    
    Returns: Nothing, but saves a file named fileOUT at the specified location pathOUT.
    
    @param fileOUT: name you wish to give to the output file (.txt NO need to be specified)
    @type fileOUT: string
    @param pathOUT: path to the location where you wish to save the output
    @type pathOUT: string
    
    @param time: heliocentric JD of the observations [d]
    @type time: numpy array of length N
    @param flux: raw, uncorrected flux measurements [adu]
    @type flux: numpy array of length N
    @param xPOS: CCD position measurements along the x position axis [pixel]
    @type xPOS: numpy array of length N
    @param yPOS: CCD position measurements along the x position axis [pixel]
    @type yPOS: numpy array of length N
    @param temperature: temperature measurements [deg]
    @type temperature: numpy array of length N
    
    @param exposureTIME: exposure time for the observations [s]
    @type exposureTIME: numpy array of length N
    @param numberSTACKS: number of onboard stacked datapoints corresponding to this datapoint []
    @type numberSTACKS: numpy array of length N
    
    @param correction: correction for the flux to detrend it [adu]
    @type correction: numpy array of length N
    """
    # Checking if the last character of pathOUT is an '/'
    if not(pathOUT[-1] == '/'):
      pathOUT += '/'
    # Checking if the suffix of the file is given
    if not fileOUT[-4:] in ['.txt', '.dat']:
      fileOUT += '.dat' 
      
    # Preparing the header of the output file
    headerSTRING = 'BRITE photometry, which was detrended for instrumental effects on ' + strftime("%Y-%m-%d %H:%M:%s") + '.'
    headerSTRING +='\n----------------------------------------'
    headerSTRING +='\nColumn1: time measurements [d]'
    headerSTRING +='\nColumn2: raw, uncorrected flux [adu]'
    headerSTRING +='\nColumn3: CCD centroid x-position [pixel]'
    headerSTRING +='\nColumn4: CCD centroid y-position [pixel]'
    headerSTRING +='\nColumn5: CCD temperature [deg]'
    headerSTRING +='\nColumn6: exposure time of the observations [s]'
    headerSTRING +='\nColumn7: number of stacked observations corresponding to one datapoint []'
    headerSTRING +='\nColumn8: final correction to apply to the raw flux [adu]'
    headerSTRING +='\n----------------------------------------'
    
    # Constructing the matrix
    dtOUT = np.dtype([('time', np.float), ('flux', np.float), ('xPOS', np.float), ('yPOS', np.float), ('temperature', np.float), ('exposureTIME', np.float), ('numberSTACKS', np.float), ('correction', np.float)])
    matrixOUT = np.zeros(len(time), dtype=dtOUT)
    matrixOUT['time'] = time; matrixOUT['flux'] = fluxRAW; matrixOUT['xPOS'] = xPOS; matrixOUT['yPOS'] = yPOS; matrixOUT['temperature'] = temperature; matrixOUT['exposureTIME'] = exposureTIME; matrixOUT['numberSTACKS'] = numberSTACKS; matrixOUT['correction'] = correction
    
    # The actual saving using a numpy.savetxt    
    np.savetxt(pathOUT + fileOUT, matrixOUT, fmt=('%.12e %.7f %.4f %.4f %.4f %.2f %i %.7f'), delimiter=' ', header=headerSTRING, comments='#')    
