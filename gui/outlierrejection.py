# -*- coding: utf-8 -*-
"""
Routines to do an interactive outlier rejection for a given parameter array, varying over time. This is not fully generalised, since you can also delete all parameters within (roughly) one orbit passage. However, it should provide you with an example on how to make matplotlib GUIs.

WARNING: This is on average still a messy routine, since I make global parameters to pass them through the different routines. This is to somewhat circumvent the pl.show(), where it becomes difficult to pass things through when it is called.

DANGER: Not much commenting is done for these routines, until they are properly tested.

TODO: enable a kwarg to choose between a lowess filter and a scipy.interpolate.splrep for the determination of the 'trend'
    
Last update 18 March 2016


@author: Bram Buysschaert
"""

#===============================================================================
# 				Packages
#===============================================================================
# import CheckMatplotlib # LOCAL ROUTINE, remove
import numpy as np

import pylab as pl

import scipy.interpolate as scInterp

from matplotlib.ticker import NullFormatter, MultipleLocator, FormatStrFormatter, NullLocator
from matplotlib import patches, rc, rcParams, pyplot
import matplotlib.gridspec as gridspec
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
def lookwithinORBIT(time, idxVALUE, **kwargs):
      """
      Look for all the values within the same orbit, assuming that it would be closer than 50min from the current datapoint.
      Return all the indexes of those values.
      """
      halfPERIOD = 50#min - rough value DANGER DANGER should be made a kwarg (unsure of mpl_connect passes them through though)
      halfPERIOD = halfPERIOD/ (24. * 60.)
      
      timeVALUE = time[idxVALUE]
      timeMIN, timeMAX = timeVALUE-halfPERIOD, timeVALUE+halfPERIOD
      #idxORBIT = np.where(((time>=timeMIN) and (time<=timeMAX)))[0] #should give you an array of elements within that orbit #UNCLEAR WHY IT DOES NOT WORK ISSUES WITH ALL / ANY
      idxORBIT = []
      for ii in range(idxVALUE-200,idxVALUE+200,1):
	if ii >= len(time):
	  break
	elif ii < 0:
	  break
	elif (time[ii] >=timeMIN) and (time[ii] <= timeMAX):
	  idxORBIT.append(int(ii))
      idxORBIT = np.array(idxORBIT)
      return idxORBIT
    
def ontypeEVENT(event, **kwargs):
    """
    Defining the routine which handles the ontype events.
    """
    global timeGLOBAL, parameterARRAYGLOBAL, idx_OUTLIERS, figINToutlier, workORBIT, parameterSPLINE
    if event.inaxes is  None:
      print 'You are in the wrong figure.'
      return
    else:
      #Choosing the correct subfigure
      axEVENT = event.inaxes
      timeEVENT = event.xdata
      if axEVENT.is_last_row() == True:
	print 'Please work in the parameter figure, not the residuals.'
	return
      
    if event.key == 'H' or event.key == 'h':
      print 'This is the help function of the interactiveOUTLIER routine.'
      print 'It allows you to mark datapoints of the current parameter array as outliers and use these in the outlier rejection'
      print '\th: help function -- this'
      print '\to: work in full orbit or not; currently you have a ' + str(workORBIT)
      print '\ta: add a given datapoint as an outlier'
      print '\td: delete a given datapoint as an outlier'
      print '\tp: print all the identified outliers'
      print '\tq: quit'
      return
    if event.key == 'O' or event.key == 'o':
      if workORBIT == False:
	print 'You are now editing full orbits'
	workORBIT = True
      elif workORBIT == True:
	print 'You are no longer editing full orbits'
	workORBIT = False
	return
    if event.key == 'A' or event.key == 'a':
      indexADD, timeADD = min(enumerate(timeGLOBAL), key=lambda z: abs(z[1]-timeEVENT))
      print 'The closest value I found next to your cursor would be {:5.5f}'.format(timeADD)
      if workORBIT:
	indexORBIT = lookwithinORBIT(timeGLOBAL, indexADD)
	indexADD = indexORBIT
      idx_OUTLIERS = np.append(idx_OUTLIERS, indexADD)
      idx_OUTLIERS = np.unique(idx_OUTLIERS) # catches any multiple additions
      #Updating the figure
      event.canvas.figure.clear()
      axPARAMETER = figINToutlier.add_subplot(211)
      axRESIDUAL =  figINToutlier.add_subplot(212, sharex=axPARAMETER)
      axPARAMETER.plot(timeGLOBAL, parameterARRAYGLOBAL, 'k.', ms=5, alpha=.4)
      axPARAMETER.plot(timeGLOBAL, parameterSPLINE, 'b-', lw=3, alpha=.4)
      axRESIDUAL.plot(timeGLOBAL, parameterARRAYGLOBAL-parameterSPLINE, 'k.', ms=5, alpha=.4)
      try:
	axPARAMETER.plot(timeGLOBAL[idx_OUTLIERS], parameterARRAYGLOBAL[idx_OUTLIERS], 'rx', ms=15, alpha=.6)
	axRESIDUAL.plot(timeGLOBAL[idx_OUTLIERS], parameterARRAYGLOBAL[idx_OUTLIERS]-parameterSPLINE[idx_OUTLIERS], 'rx', ms=15, alpha=.6)
      except:
	pass
      
      event.canvas.draw()
      return
    if event.key == 'D' or event.key == 'd':
      indexDELETE, timeDELETE = min(enumerate(timeGLOBAL), key=lambda z: abs(z[1]-timeEVENT))
      #Check if it is within a certain tolerance
      if abs(timeDELETE-timeEVENT) > 0.1: #tolerance value
	print 'The closest value I found to *delete* would be {:5.5f}, but it is too far away from your cursor at {:5.5f}'.format(timeDELETE, timeEVENT)
	return
      else:
	if workORBIT:
	  indexORBIT = lookwithinORBIT(timeGLOBAL, indexDELETE)
	  indexDELETE = indexORBIT
	idxREMOVEoutliers = []
	try:
	  for oo in range(len(indexDELETE)):
	    idx, value = min(enumerate(idx_OUTLIERS), key=lambda z: abs(z[1]-indexDELETE[oo]))
	    idxREMOVEoutliers.append(int(idx))
	except:
	  idx, value = min(enumerate(idx_OUTLIERS), key=lambda z: abs(z[1]-indexDELETE)) # no [oo]
	  idxREMOVEoutliers.append(int(idx))
	idxREMOVEoutliers = np.array(idxREMOVEoutliers)
	idx_OUTLIERS = np.delete(idx_OUTLIERS, idxREMOVEoutliers)
	idx_OUTLIERS = np.unique(idx_OUTLIERS) # catches any multiple additions
	#Updating the figure
	event.canvas.figure.clear()
	axPARAMETER = figINToutlier.add_subplot(211)
	axRESIDUAL =  figINToutlier.add_subplot(212, sharex=axPARAMETER)
	axPARAMETER.plot(timeGLOBAL, parameterARRAYGLOBAL, 'k.', ms=5, alpha=.4)
	axPARAMETER.plot(timeGLOBAL, parameterSPLINE, 'b-', lw=3, alpha=.4)
	axRESIDUAL.plot(timeGLOBAL, parameterARRAYGLOBAL-parameterSPLINE, 'k.', ms=5, alpha=.4)
	try:
	  axPARAMETER.plot(timeGLOBAL[idx_OUTLIERS], parameterARRAYGLOBAL[idx_OUTLIERS], 'rx', ms=15, alpha=.6)
	  axRESIDUAL.plot(timeGLOBAL[idx_OUTLIERS], parameterARRAYGLOBAL[idx_OUTLIERS]-parameterSPLINE[idx_OUTLIERS], 'rx', ms=15, alpha=.6)
	except:
	  pass
	event.canvas.draw()
	return
    if event.key == 'P' or event.key == 'p':	
      print idx_OUTLIERS
    if event.key == 'Q' or event.key == 'q':	
      pl.close(figINToutlier)
      return     
    
def interactiveOUTLIER(time, parameterARRAY, **kwargs):
    """
    This routine allows you to perform an interactive outlier selection of a certain parameter array.
    The figure will show you the parameter varying over time, a simple spline fit to it and the residuals compared to this fit.
    You will be able to select individual outlier points or mark full orbits as an outlier
    """     
    #global time, parameterARRAY, idx_OUTLIERS, figINToutlier, workORBIT, parameterSPLINE
    global timeGLOBAL, parameterARRAYGLOBAL, idx_OUTLIERS, figINToutlier, workORBIT, parameterSPLINE
    timeGLOBAL, parameterARRAYGLOBAL, workORBIT, idx_OUTLIERS = np.copy(time - time[0]), np.array(parameterARRAY), False, np.array([],dtype='int32')
    
    timeGLOBAL -= timeGLOBAL[0]
    #REST OF THE ROUTINE
    timeKNOTPOINTS = kwargs.get('knotpoints', np.arange(min(timeGLOBAL), max(timeGLOBAL), 5.))
    tckPARAMETER = scInterp.splrep(timeGLOBAL, parameterARRAY, t=timeKNOTPOINTS[1:-1], k=3)
    parameterSPLINE = scInterp.splev(timeGLOBAL,tckPARAMETER)
    
    figINToutlier = pl.figure(figsize=(16,16))
    axPARAMETER = figINToutlier.add_subplot(211)
    axRESIDUAL =  figINToutlier.add_subplot(212, sharex=axPARAMETER)
    axPARAMETER.plot(timeGLOBAL, parameterARRAY, 'k.', ms=5, alpha=.4)
    axPARAMETER.plot(timeGLOBAL, parameterSPLINE, 'b-', lw=3, alpha=.4)
    try:
      axPARAMETER.plot(HJD_MIDOBS[idx_OUTLIERS], parameterARRAY[idx_OUTLIERS], 'rx', ms=15, alpha=.6)
    except:
      pass
    axRESIDUAL.plot(timeGLOBAL, parameterARRAY-parameterSPLINE, 'k.', ms=5, alpha=.4)
    figINToutlier.canvas.mpl_connect('key_press_event',ontypeEVENT)
    
    
    pl.show() #DANGER pl.show is implied here. Doing it any later will most likely not work, since the return statement for idx_OUTLIERS will be before the pl.show
    
    return idx_OUTLIERS
