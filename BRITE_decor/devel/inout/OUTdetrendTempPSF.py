# -*- coding: utf-8 -*-
"""
Routines to save the output of the detrending process for the temperature dependent PSF changes, which produce instrumental effects between flux and position.
    
Last update 16 March 2016

TODO: Make it mandatory to provide the used kwargs, saving them in the header for future reference and / or comparison.

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
def cleanTCKstring(TCKstring, **kwargs):
    """
    Routine to reconvert a TCK from the scInterp.splrep, which was saved as a string, back to an useful tuple.
    
    The following properties are known of the string:
    -The first 8 elements should be '(array(['
    -The second to last element should contain the order of the spline, i.e. "k)']". --> Can stop at element -3
    -Once the first numpy array is retrieved, we should find a '], array(['
    -There are '\n' saved in the string. --> Need to be removed
    
    Returns: The TCK, but as a tuple.
    
    @param TCKstring: string of the tck
    @type TCKstring: numpy.string_
    
    @return TCKstringCLEAN: the cleaned tck
    @rtype: string
    """
    
    # Convert the (annoying) numpy.string_ to a 'classical' string
    # NOTE this has some extra effects
    TCKstring = str(TCKstring)
    
    TCKstringCLEAN = TCKstring.replace(' ', '').replace('\\n', '').replace('\n', '')
    return TCKstringCLEAN

def OUT_DIAGNOSTIC_detrendTEMPpsfFULL(fileOUT, pathOUT, time, endINDEXbin, endINDEXbin_reason, tckFIRST, tckSECOND, diagnosticCORRECTION, **kwargs):
    """
    Routine to save the output of the detrending and correction process of the temperature-dependent PSF changes. As such, you can keep track of what the different fits were, which issues you might have had, and most importantly, what the indexes of your different time bins are.
    
    WARNING: These indexes will refer to the current length N of your detrended flux / time arrays. Be sure to update these, if you need them, whenever you perform clipping in a future step. To help you with this, we also include the time[endINDEXbin] in the file
    
    NOTE for a more simple example on how to save strings and floats to a text file, look on:
    http://stackoverflow.com/questions/16621351/how-to-use-python-numpy-savetxt-to-write-strings-and-float-number-to-an-ascii-fi (the non-accepted answer!)
    
    NOTE You will need a dedicated routine to open the output of this routine. Most likely, you will have to work with a file.open('r') and a file.readline(), because of the string / TCK.
    
    TODO Get your kwargs and dump them into the header.
        
    Returns: Nothing, but saves a file named fileOUT at the specified location pathOUT.
    
    @param fileOUT: name you wish to give to the output file (.txt NO need to be specified)
    @type fileOUT: string
    @param pathOUT: path to the location where you wish to save the output
    @type pathOUT: string
    @param time: time measurements [d]
    @type time: numpy array of length N
    @param endINDEXbin: end index (of your time array for that bin)
    @type endINDEXbin: numpy array (dtype='int32') of length K
    @param endINDEXbin_reason: reason to end your temperature bin
    @type endINDEXbin_reason: list of length K (containing strings)
    @param tckFIRST: TCK of the first correction, in string format
    @type tckFIRST: list with numpy.str_ of length K
    @param tckSECOND: TCK of the second correction, in string format
    @type tckSECOND: list with numpy.str_ of length K
    @param diagnosticCORRECTION: diagnostic to trace the corrections
    @type diagnosticCORRECTION: numpy array (dtype='int32') of length K
    
    @kwargs: SPLINEtckLENGTH: number of bytes allocated for the to-string-converted TCK - default is 2000 [integer]
    """
    # Reading in the kwargs
    SPLINEtckLENGTH = 'S' + str(np.int(kwargs.get('SPLINEstringLENGTH', 2000))) # [string]
    
    
    # Checking if the last character of pathOUT is an '/'
    if not(pathOUT[-1] == '/'):
      pathOUT += '/'
    # Checking if the suffix of the file is given
    if not fileOUT[-4:] in ['.txt', '.dat']:
      fileOUT += '.dat'  
    
    # Preparing the header of the output file
    headerSTRING = 'Detrending for the temperature dependent PSF, observed in BRITE photometry, done on ' + strftime("%Y-%m-%d %H:%M:%s") + '.'
    headerSTRING +='\nThis file contains the information of the temperature bins and the different fits used during the detrending process.'
    headerSTRING +='\nThe adopted technique is to:'
    headerSTRING +='\n1. split the observations in smaller time bins, which have been selected according to the long term temperature variations and (possibly) too long gaps.'
    headerSTRING +='\n2. perform 1D fits within each bin for flux vs xPOS and flux vs yPOS -- !!always both!! --.'
    headerSTRING +='\nNow, on to the header for each columnn:'
    headerSTRING +='\n----------------------------------------'
    headerSTRING +='\nColumn1: end index of the complete array, for that bin (include index!) [integer]'
    headerSTRING +='\nColumn2: reason for the ending of that bin (include index!) [string]; tem = temperature difference; gap = gapsize too large; end = end of the file'
    headerSTRING +='\nColumn3: end time of that bin [d]'
    headerSTRING +='\nColumn4: time span of that bin [d]'
    headerSTRING +='\nColumn5: diagnostic value of the fitting process within that bin [integer]'
    headerSTRING +='\nColumn6: full TCK (scipy.interpolate.splrep) of the most optimal fit for the first coordinate [string]'
    headerSTRING +='\nColumn7: full TCK (scipy.interpolate.splrep) of the most optimal fit for the second coordinate [string]'
    headerSTRING +='\n----------------------------------------'
    
    
    # Calculating the time span of each bin
    timeSPAN = []
    for kk in range(len(endINDEXbin)):
      if kk == 0:
	timeSPAN.append(float(time[endINDEXbin[kk]] - time[0]))
      else:
	timeSPAN.append(float(time[endINDEXbin[kk]] - time[endINDEXbin[kk-1]+1]))
    
    
    # Cleaning the TCK strings (mainly getting rid of the '\n', in addition, we remove all spaces'
    tckFIRSTclean, tckSECONDclean = [], []
    for kk in range(len(endINDEXbin)):
      tckFIRSTclean.append(cleanTCKstring(tckFIRST[kk]))
      tckSECONDclean.append(cleanTCKstring(tckSECOND[kk]))
    
    # Constructing the matrix
    dtOUT = np.dtype([('endINDEXbin', np.int32), ('endINDEXbin_reason', 'S3'), ('timeENDbin', np.float), ('timeSPAN', np.float), ('diagnosticCORRECTION', np.int32), ('tckFIRST', SPLINEtckLENGTH), ('tckSECOND', SPLINEtckLENGTH)])
    matrixOUT = np.zeros(len(endINDEXbin), dtype=dtOUT)
    matrixOUT['endINDEXbin'] = endINDEXbin; matrixOUT['endINDEXbin_reason'] = endINDEXbin_reason; matrixOUT['timeENDbin'] = time[endINDEXbin]; matrixOUT['timeSPAN'] = timeSPAN; matrixOUT['diagnosticCORRECTION'] = diagnosticCORRECTION; matrixOUT['tckFIRST'] = tckFIRSTclean; matrixOUT['tckSECOND'] = tckSECONDclean
    
    # The actual saving using a numpy.savetxt    
    np.savetxt(pathOUT + fileOUT, matrixOUT, fmt=('%i %3s %.12e %3.4f %i %s %s'), delimiter=' ', header=headerSTRING, comments='#')

    return
    
def OUT_FLUX_detrendTEMPpsfFULL(fileOUT, pathOUT, time, flux, fluxCORRECTED, CORRECTION, **kwargs):
    """
    Routine to save the output of the detrending and correction process of the temperature-dependent PSF changes. More particularly, the corrected flux. This can then be passed on to further detrending routines.
    
    
    Returns: Nothing, but saves a file named fileOUT at the specified location pathOUT.
    
    @param fileOUT: name you wish to give to the output file (.txt NO need to be specified)
    @type fileOUT: string
    @param pathOUT: path to the location where you wish to save the output
    @type pathOUT: string
    @param time: time measurements [d]
    @type time: numpy array of length N
    @param flux: flux measurements [adu]
    @type flux: numpy array of length N
    @param fluxCORRECTED: corrected photometry [adu]
    @type fluxCORRECTED: numpy array of length N
    @param CORRECTION: correction to apply to the flux [adu]
    @type CORRECTION: numpy array of length N
    """
    # Checking if the last character of pathOUT is an '/'
    if not(pathOUT[-1] == '/'):
      pathOUT += '/'
    # Checking if the suffix of the file is given
    if not fileOUT[-4:] in ['.txt', '.dat']:
      fileOUT += '.dat'  
    
    # Preparing the header of the output file
    headerSTRING = 'Detrending for the temperature dependent PSF, observed in BRITE photometry, done on ' + strftime("%Y-%m-%d %H:%M:%s") + '.'
    headerSTRING +='\nThis file contains the corrected flux and the correction which has been applied in order to deduce it.'
    headerSTRING +='\nThe adopted technique is to:'
    headerSTRING +='\n1. split the observations in smaller time bins, which have been selected according to the long term temperature variations and (possibly) too long gaps.'
    headerSTRING +='\n2. perform 1D fits within each bin for flux vs xPOS and flux vs yPOS -- !!always both!! --.'
    headerSTRING +='\nNow, on to the header for each columnn:'
    headerSTRING +='\n----------------------------------------'
    headerSTRING +='\nColumn1: time measurements [d]'
    headerSTRING +='\nColumn2: original, uncorrected flux [adu]'
    headerSTRING +='\nColumn3: corrected flux [adu]'
    headerSTRING +='\nColumn4: derived correction, applied to the flux [adu]'
    headerSTRING +='\n----------------------------------------'
    
    # Constructing the matrix
    dtOUT = np.dtype([('time', np.float), ('flux', np.float), ('fluxCORRECTED', np.float), ('CORRECTION', np.float)])
    matrixOUT = np.zeros(len(time), dtype=dtOUT)
    matrixOUT['time'] = time; matrixOUT['flux'] = flux; matrixOUT['fluxCORRECTED'] = fluxCORRECTED; matrixOUT['CORRECTION'] = CORRECTION
    
    # The actual saving using a numpy.savetxt    
    np.savetxt(pathOUT + fileOUT, matrixOUT, fmt=('%.12e %.7e %.7e %.7e'), delimiter=' ', header=headerSTRING, comments='#')
    
    return