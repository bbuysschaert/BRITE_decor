# -*- coding: utf-8 -*-
"""
Routine to cut a given setup file into smaller chunks. This can be useful if you want to do a proper PSFdetrend.

WARNING: This is on average still a messy routine, since I make global parameters to pass them through the different routines. This is to somewhat circumvent the pl.show(), where it becomes difficult to pass things through when it is called.

DANGER: Not much commenting is done for these routines, until they are properly tested.

TODO: enable a kwarg to choose between a lowess filter and a scipy.interpolate.splrep for the determination of the 'trend'
    
Last update 06 April 2016


@author: Bram Buysschaert
"""

#===============================================================================
# 				Packages
#===============================================================================
import numpy as np

import pylab as pl

import scipy.interpolate as scInterp

from matplotlib.ticker import NullFormatter, MultipleLocator, FormatStrFormatter, NullLocator
from matplotlib import patches, rc, rcParams, pyplot
import matplotlib.gridspec as gridspec

import statsmodels.api as sm

import BRITE_decor.inout.save as saveBRITE
#===============================================================================
# 				Settings
#===============================================================================
from string import whitespace as ws
import re
rc('xtick', labelsize=20) 
rc('ytick', labelsize=20) 
rc("lines", markeredgewidth=2.0)
rc("axes", linewidth=2.0)
#rc('text', usetex=True)
rcParams["font.size"] = 25

p = re.compile('(%s)' % ('|'.join([c for c in ws])))
s = " \ nu "
nu = p.sub('',s)
#===============================================================================
# 				Code
#===============================================================================
def ontypeEVENT(event):
    """
    Defining the routine which handles the ontype events.
    """
    global fileOUT_global, pathOUT_global, time_global, fluxRAW_global, xPOS_global, yPOS_global, temperature_global, exposureTIME_global, numberSTACKS_global
    global timeLONG_global, fluxLONG_global, temperatureLONG_global, xPOSLONG_global, yPOSLONG_global
    global figCUTTING
    global idxCUT, idxEXPORTED, kwargs_global
    
    if event.inaxes is  None:
      print 'You are in the wrong figure.'
      return
    else:
      timeEVENT = event.xdata
      
    if event.key == 'H' or event.key == 'h':
      print 'This is the help function of the interactiveCUTsetup routine.'
      print 'This routine permits you to cut up a setup file into smaller chunks.'
      print '\th: help function -- this'
      print '\t1: define the beginning of the region you wish to cut out'
      print '\t2: define the end of the region you wish to cut out'
      print '\te: extract the region you wish to cut out'
      #print '\tr: reset to the original setup files' # Not really possible.
      print '\tq: quit'
      return
    
    if event.key == '1':
      indexFIRST, timeFIRST = min(enumerate(time_global), key=lambda z: abs(z[1]-timeEVENT))
      idxCUT[0] = indexFIRST
      
      print 'Adding index {:.0f} as first index of the cut'.format(indexFIRST)
      
      event.canvas.figure.clear()
      
      axFLUX = figCUTTING.add_subplot(411)
      axTEMP = figCUTTING.add_subplot(412, sharex=axFLUX)
      axxPOS = figCUTTING.add_subplot(413, sharex=axFLUX)
      axyPOS = figCUTTING.add_subplot(414, sharex=axFLUX)
      
      axFLUX.plot(time_global, fluxRAW_global, 'k.', ms=6, alpha=.3); axFLUX.plot(timeLONG_global, fluxLONG_global, 'b.', ms=6, alpha=.6); axFLUX.grid(color='r', lw=3, alpha=.3)
      axTEMP.plot(time_global, temperature_global, 'k.', ms=6, alpha=.3); axTEMP.plot(timeLONG_global, temperatureLONG_global, 'b.', ms=6, alpha=.6); axTEMP.grid(color='r', lw=3, alpha=.3)
      axxPOS.plot(time_global, xPOS_global, 'k.', ms=6, alpha=.3); axxPOS.plot(timeLONG_global, xPOSLONG_global, 'b.', ms=6, alpha=.6); axxPOS.grid(color='r', lw=3, alpha=.3)
      axyPOS.plot(time_global, yPOS_global, 'k.', ms=6, alpha=.3); axyPOS.plot(timeLONG_global, yPOSLONG_global, 'b.', ms=6, alpha=.6); axyPOS.grid(color='r', lw=3, alpha=.3)
      try:
	axFLUX.plot(time_global[idxEXPORTED], fluxRAW_global[idxEXPORTED], 'rx', ms=6, alpha=.3); axTEMP.plot(time_global[idxEXPORTED], temperature_global[idxEXPORTED], 'rx', ms=6, alpha=.3); axxPOS.plot(time_global[idxEXPORTED], xPOS_global[idxEXPORTED], 'rx', ms=6, alpha=.3); axyPOS.plot(time_global[idxEXPORTED], yPOS_global[idxEXPORTED], 'rx', ms=6, alpha=.3)
      except:
	pass
      
      axFLUX.axvline(x=time_global[idxCUT[0]], color ='g', lw=5, alpha=.4); axTEMP.axvline(x=time_global[idxCUT[0]], color ='g', lw=5, alpha=.4); axxPOS.axvline(x=time_global[idxCUT[0]], color ='g', lw=5, alpha=.4);  axyPOS.axvline(x=time_global[idxCUT[0]], color ='g', lw=5, alpha=.4)
      axFLUX.axvline(x=time_global[idxCUT[1]], color ='m', lw=5, alpha=.4); axTEMP.axvline(x=time_global[idxCUT[1]], color ='m', lw=5, alpha=.4); axxPOS.axvline(x=time_global[idxCUT[1]], color ='m', lw=5, alpha=.4);  axyPOS.axvline(x=time_global[idxCUT[1]], color ='m', lw=5, alpha=.4)
      
      axFLUX.set_ylabel('Flux [adu]'); axTEMP.set_ylabel('Temp [deg]'); axxPOS.set_ylabel('xPOS [pix]'); axyPOS.set_ylabel('yPOS [pix]'); axyPOS.set_xlabel('Time [d]')
      event.canvas.draw()
      return
      
    if event.key == '2':
      indexSECOND, timeSECOND = min(enumerate(time_global), key=lambda z: abs(z[1]-timeEVENT))
      idxCUT[1] = indexSECOND
      
      print 'Adding index {:.0f} as second index of the cut'.format(indexSECOND)
      
      event.canvas.figure.clear()
      
      axFLUX = figCUTTING.add_subplot(411)
      axTEMP = figCUTTING.add_subplot(412, sharex=axFLUX)
      axxPOS = figCUTTING.add_subplot(413, sharex=axFLUX)
      axyPOS = figCUTTING.add_subplot(414, sharex=axFLUX)
      
      axFLUX.plot(time_global, fluxRAW_global, 'k.', ms=6, alpha=.3); axFLUX.plot(timeLONG_global, fluxLONG_global, 'b.', ms=6, alpha=.6); axFLUX.grid(color='r', lw=3, alpha=.3)
      axTEMP.plot(time_global, temperature_global, 'k.', ms=6, alpha=.3); axTEMP.plot(timeLONG_global, temperatureLONG_global, 'b.', ms=6, alpha=.6); axTEMP.grid(color='r', lw=3, alpha=.3)
      axxPOS.plot(time_global, xPOS_global, 'k.', ms=6, alpha=.3); axxPOS.plot(timeLONG_global, xPOSLONG_global, 'b.', ms=6, alpha=.6); axxPOS.grid(color='r', lw=3, alpha=.3)
      axyPOS.plot(time_global, yPOS_global, 'k.', ms=6, alpha=.3); axyPOS.plot(timeLONG_global, yPOSLONG_global, 'b.', ms=6, alpha=.6); axyPOS.grid(color='r', lw=3, alpha=.3)
      axFLUX.axvline(x=time_global[idxCUT[0]], color ='g', lw=5, alpha=.4); axTEMP.axvline(x=time_global[idxCUT[0]], color ='g', lw=5, alpha=.4); axxPOS.axvline(x=time_global[idxCUT[0]], color ='g', lw=5, alpha=.4);  axyPOS.axvline(x=time_global[idxCUT[0]], color ='g', lw=5, alpha=.4)
      axFLUX.axvline(x=time_global[idxCUT[1]], color ='m', lw=5, alpha=.4); axTEMP.axvline(x=time_global[idxCUT[1]], color ='m', lw=5, alpha=.4); axxPOS.axvline(x=time_global[idxCUT[1]], color ='m', lw=5, alpha=.4);  axyPOS.axvline(x=time_global[idxCUT[1]], color ='m', lw=5, alpha=.4)
      try:
	axFLUX.plot(time_global[idxEXPORTED], fluxRAW_global[idxEXPORTED], 'rx', ms=6, alpha=.3); axTEMP.plot(time_global[idxEXPORTED], temperature_global[idxEXPORTED], 'rx', ms=6, alpha=.3); axxPOS.plot(time_global[idxEXPORTED], xPOS_global[idxEXPORTED], 'rx', ms=6, alpha=.3); axyPOS.plot(time_global[idxEXPORTED], yPOS_global[idxEXPORTED], 'rx', ms=6, alpha=.3)
      except:
	pass
      
      axFLUX.set_ylabel('Flux [adu]'); axTEMP.set_ylabel('Temp [deg]'); axxPOS.set_ylabel('xPOS [pix]'); axyPOS.set_ylabel('yPOS [pix]'); axyPOS.set_xlabel('Time [d]')
      event.canvas.draw()  
      return
      
      
    if event.key == 'E' or event.key == 'e':
      indexTOexport = np.arange(idxCUT[0], idxCUT[1], 1, dtype='int32')
      fileOUT = fileOUT_global + '_{:3.3f}_'.format(time_global[idxCUT[0]] - time_global[0]).replace('.', 'p').replace(' ', '0') + '{:3.3f}'.format(time_global[idxCUT[1]] - time_global[0]).replace('.', 'p').replace(' ', '0') + '.dat'
      
      idxEXPORTED = np.append(idxEXPORTED, indexTOexport)
      idxEXPORTED = np.unique(idxEXPORTED) # just in case
      
      print 'Trying to export your cut to ' + fileOUT
      
      dataTYPE = kwargs_global.get('TYPEdata', 'clipped')
      if not dataTYPE in ['detrend', 'clipped']:
	raise ValueError, 'dataTYPE kwarg not understood...'
      
      if dataTYPE == 'clipped':
	saveBRITE.saveCLIPPED(fileOUT, pathOUT_global, time_global[indexTOexport], fluxRAW_global[indexTOexport], xPOS_global[indexTOexport], yPOS_global[indexTOexport], temperature_global[indexTOexport], exposureTIME_global[indexTOexport], numberSTACKS_global[indexTOexport])
      elif dataTYPE == 'detrend':
	fluxCORRECTION = kwargs_global.get('CORRECTIONflux')
	saveBRITE.saveDETREND(fileOUT, pathOUT_global, time_global[indexTOexport], fluxRAW_global[indexTOexport], xPOS_global[indexTOexport], yPOS_global[indexTOexport], temperature_global[indexTOexport], exposureTIME_global[indexTOexport], numberSTACKS_global[indexTOexport], fluxCORRECTION[indexTOexport])
      
      idxCUT = np.array([0,0], dtype='int32')
      
      event.canvas.figure.clear()
      
      axFLUX = figCUTTING.add_subplot(411)
      axTEMP = figCUTTING.add_subplot(412, sharex=axFLUX)
      axxPOS = figCUTTING.add_subplot(413, sharex=axFLUX)
      axyPOS = figCUTTING.add_subplot(414, sharex=axFLUX)
      
      axFLUX.plot(time_global, fluxRAW_global, 'k.', ms=6, alpha=.3); axFLUX.plot(timeLONG_global, fluxLONG_global, 'b.', ms=6, alpha=.6); axFLUX.grid(color='r', lw=3, alpha=.3)
      axTEMP.plot(time_global, temperature_global, 'k.', ms=6, alpha=.3); axTEMP.plot(timeLONG_global, temperatureLONG_global, 'b.', ms=6, alpha=.6); axTEMP.grid(color='r', lw=3, alpha=.3)
      axxPOS.plot(time_global, xPOS_global, 'k.', ms=6, alpha=.3); axxPOS.plot(timeLONG_global, xPOSLONG_global, 'b.', ms=6, alpha=.6); axxPOS.grid(color='r', lw=3, alpha=.3)
      axyPOS.plot(time_global, yPOS_global, 'k.', ms=6, alpha=.3); axyPOS.plot(timeLONG_global, yPOSLONG_global, 'b.', ms=6, alpha=.6); axyPOS.grid(color='r', lw=3, alpha=.3)
      try:
	axFLUX.plot(time_global[idxEXPORTED], fluxRAW_global[idxEXPORTED], 'rx', ms=6, alpha=.3); axTEMP.plot(time_global[idxEXPORTED], temperature_global[idxEXPORTED], 'rx', ms=6, alpha=.3); axxPOS.plot(time_global[idxEXPORTED], xPOS_global[idxEXPORTED], 'rx', ms=6, alpha=.3); axyPOS.plot(time_global[idxEXPORTED], yPOS_global[idxEXPORTED], 'rx', ms=6, alpha=.3)
      except:
	pass
      
      axFLUX.set_ylabel('Flux [adu]'); axTEMP.set_ylabel('Temp [deg]'); axxPOS.set_ylabel('xPOS [pix]'); axyPOS.set_ylabel('yPOS [pix]'); axyPOS.set_xlabel('Time [d]')
      event.canvas.draw()
      return
    
    if event.key == 'Q' or event.key == 'q':	
      pl.close(figCUTTING)
      return     
    
def interactiveCUTsetup(fileOUT, pathOUT, time, fluxRAW, xPOS, yPOS, temperature, exposureTIME, numberSTACKS, **kwargs):
    """
    This routine permits you to cut up a setup file into smaller chunks. It is made as such, that it should handle both the data from a loadCLIPPED and a loadDETREND.
    NOTE: You cannot pass kwargs through mpl_connect
    """     
    
    # Declaring the global values
    global fileOUT_global, pathOUT_global, time_global, fluxRAW_global, xPOS_global, yPOS_global, temperature_global, exposureTIME_global, numberSTACKS_global
    global timeLONG_global, fluxLONG_global, temperatureLONG_global, xPOSLONG_global, yPOSLONG_global
    global figCUTTING
    global idxCUT, idxEXPORTED, kwargs_global
    idxCUT, idxEXPORTED = np.array([0,0],dtype='int32'), np.array([], dtype='int32')
    kwargs_global = kwargs
    
    fileOUT_global, pathOUT_global, time_global, fluxRAW_global, xPOS_global, yPOS_global, temperature_global, exposureTIME_global, numberSTACKS_global = fileOUT, pathOUT, time, fluxRAW, xPOS, yPOS, temperature, exposureTIME, numberSTACKS
    
    # Doing the lowess filtering (for guidance)
    lowessFRAC = kwargs.get('LOWESSfrac', 0.05) # []
    fluxLOWESS = sm.nonparametric.lowess(fluxRAW_global, time_global, frac=lowessFRAC)
    temperatureLOWESS = sm.nonparametric.lowess(temperature_global, time_global, frac=lowessFRAC)
    xPOSLOWESS = sm.nonparametric.lowess(xPOS_global, time_global, frac=lowessFRAC)
    yPOSLOWESS = sm.nonparametric.lowess(yPOS_global, time_global, frac=lowessFRAC)
    timeLONG_global, fluxLONG_global, temperatureLONG_global, xPOSLONG_global, yPOSLONG_global = fluxLOWESS[:,0], fluxLOWESS[:,1], temperatureLOWESS[:,1], xPOSLOWESS[:,1], yPOSLOWESS[:,1]
    
    
    # Making the figure
    figCUTTING = pl.figure(figsize=(16,16))
    axFLUX = figCUTTING.add_subplot(411)
    axTEMP = figCUTTING.add_subplot(412, sharex=axFLUX)
    axxPOS = figCUTTING.add_subplot(413, sharex=axFLUX)
    axyPOS = figCUTTING.add_subplot(414, sharex=axFLUX)
    
    axFLUX.plot(time_global, fluxRAW_global, 'k.', ms=6, alpha=.3); axFLUX.plot(timeLONG_global, fluxLONG_global, 'b.', ms=6, alpha=.6); axFLUX.grid(color='r', lw=3, alpha=.3)
    axTEMP.plot(time_global, temperature_global, 'k.', ms=6, alpha=.3); axTEMP.plot(timeLONG_global, temperatureLONG_global, 'b.', ms=6, alpha=.6); axTEMP.grid(color='r', lw=3, alpha=.3)
    axxPOS.plot(time_global, xPOS_global, 'k.', ms=6, alpha=.3); axxPOS.plot(timeLONG_global, xPOSLONG_global, 'b.', ms=6, alpha=.6); axxPOS.grid(color='r', lw=3, alpha=.3)
    axyPOS.plot(time_global, yPOS_global, 'k.', ms=6, alpha=.3); axyPOS.plot(timeLONG_global, yPOSLONG_global, 'b.', ms=6, alpha=.6); axyPOS.grid(color='r', lw=3, alpha=.3)
    
    axFLUX.set_ylabel('Flux [adu]'); axTEMP.set_ylabel('Temp [deg]'); axxPOS.set_ylabel('xPOS [pix]'); axyPOS.set_ylabel('yPOS [pix]'); axyPOS.set_xlabel('Time [d]')
    
    
    figCUTTING.canvas.mpl_connect('key_press_event',ontypeEVENT)
    
    
    pl.show() #DANGER pl.show is implied here. Doing it any later will most likely not work, since the return statement for idx_OUTLIERS will be before the pl.show
    
    return