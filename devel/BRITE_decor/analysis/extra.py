"""

These functions are random and supplementary to the BRITE_decor routine, but they do not have a specific place to go quite yet, so they are all stuck here.

Last update 22 June 2016

@author = Bert Pablo

"""

#===============================================================================#      Packages
#===============================================================================

#Basic Packages

import numpy as np
import matplotlib.pyplot as plt
import copy
import scipy.interpolate as scInterp
from glob import glob
#BRITE_decor packages

import BRITE_decor.detrending.detrendPositionFlux as POSdetrendBRITE
import BRITE_decor.detrending.detrendTempFlux as TEMPdetrendBRITE
import BRITE_decor.detrending.detrendOrbitFlux as ORBdetrendBRITE

import BRITE_decor.clipping.percentageclipping as percentageclipBRITE
#===============================================================================
#   Functions
#===============================================================================
	
def map(flux, xph, yph, sizex=10, sizey=10):
    """
    make a 2D map of pixel with number of bins in each direction given sizex and size y

    @flux - flux of star
    @xph - x pixel position of center mod 1
    @yph - y pixel position of center mod 1
    @sizex - number of x bins for map
    @sizey - number of y bins for map

    """


    sizex = int(sizex)
    sizey = int(sizey)
    data = np.column_stack((flux,xph,yph))
    mask = np.zeros((sizex,sizey))
    for x in range(sizex):
#		print x
        for y in range(sizey):
            mask[x,y] = np.sum(data[:,0][(np.floor(data[:,1])==x) & (np.floor(data[:,2])==y)]/len(data[:,0][(np.floor(data[:,1])==x) & (np.floor(data[:,2])==y)]))
            if  len(data[:,0][(np.floor(data[:,1])==x) & (np.floor(data[:,2])==y)])==0:
                mask[x,y] = np.mean(flux)

    return mask




def intrem(flux, xpos, ypos):

    """
    create a x, y map of a pixel and remove it from the flux. 

    @flux- flux array
    @xpos - array of x center positions for the star psf
    @ypos - array of y center positions for the star psf
 
    @return intrapixel corrected flux
    @rtype: numpy array

    """

    fluxn = np.array([])
    xcph = xpos % 1 * 10
    ycph = ypos % 1 * 10
    mask = map(flux, xcph, ycph)
    mask = mask -np.median(mask)
    rem = np.column_stack((np.floor(xcph).astype(int), np.floor(ycph).astype(int)))
    nrem = [(x[0],x[1]) for x in rem]
#	print "NREM", nrem
#	print nrem[0]
#	print mask[nrem[0]]
    for x in range(len(flux)):
        fluxn = np.append(fluxn, flux[x] - mask[nrem[x]])
#    for x in range(len(flux)):
#        print "x1", x
#        fluxn = np.append(fluxn, flux[x] - mask[nrem[x]])
#    for x in range(len(flux)):

#        print "x2", x
#        fluxn = np.append(fluxn, flux[x] - mask[nrem[x]])
#        print x    
#    for x in range(len(flux)):
#        print x, len(flux)
##        print x, nrem[x], flux[x]
#        print fluxn
#        fluxn = np.append(fluxn, flux[x] - mask[nrem[x]])
#        print fluxn
#        print len(flux)
#	new = np.column_stack((data[:,0], fluxn, data[:,[2,-3,-2,-1]]))

    return fluxn






def binnedINTRAPIXEL(flux, xpos, ypos, binENDindexes, plot=False):

    """
    This function will loop over all the individual defined data chunks and perform an intrapixel correction on each one.
    @flux- flux array
    @xpos - array of x center positions for the star psf
    @ypos - array of y center positions for the star psf
    @binENDindexes - indexes corresponding to the end of each data chunk used to split the full array.
    @return intrapixel corrected flux
    @rtype: numpy array

    """
    init_ind = 0 # starting index
    flux_corrected = []
    for x in range(len(binENDindexes)):

        
        fin_ind = binENDindexes[x]+1 # ending index
        
        # define data chunk
        fluxchunk = flux[init_ind:fin_ind]
        xchunk = xpos[init_ind:fin_ind]
        ychunk = ypos[init_ind:fin_ind]
        
        # apply intrapixel correction
            
        fluxchunknew =  intrem(fluxchunk, xchunk, ychunk)   
       
        #plot if that's what you want to do 
        if plot:
            f, (ax1, ax2) = plt.subplots(2,1)
            ax1.set_xlabel('X Position')
            ax1.set_ylabel('Flux')
            ax2.set_ylabel('Flux')
            ax2.set_xlabel('Y Position')

            ax1.plot(xchunk, fluxchunk, 'k.', label='original')
            ax1.plot(xchunk, fluxchunknew, 'r.', ms=2, label='corrected')
            ax1.legend(prop={'size':7})
            ax2.plot(ychunk, fluxchunk, 'k.', label='original')
            ax2.plot(ychunk, fluxchunknew, 'r.', ms=2, label='corrected')
            ax2.legend(prop={'size':7})
            plt.show()                        
        # append timechunk to an array

        flux_corrected.append(fluxchunknew)
        
        # redefine init_ind

        init_ind = binENDindexes[x]+1

    flux_corrected = np.concatenate(flux_corrected)

    # flux is the only value which has changed and therefore the only value returned

    return flux_corrected


def phase_shift(phase):

    """
    Shift in phase to find longest contiguous piece. This is done so there are not large gaps in the phased data which can make problems when fitting.

    @phase- phase array (array)
    @flux- flux array (array)

    @return- phase shift, phshift which you must subtract

    """

        # sort data
    y = np.argsort(phase)
    data = phase[y]

    # double phase data so it goes from 0 to 2
    ph2 = data.copy()
    ph2 = ph2+1.0
    datan = np.concatenate((data,ph2))

    #find the largest gap in phase

    gaps = datan[1:] - datan[:-1] 
    indices = np.where(gaps==max(gaps))[0]
    
    #take data which avoids this gap if gap is greater than 1 then take original data

    if datan[indices[0]] >= 1.0:

        dataf = data
    
    elif len(indices) == 1:

        dataf = data

    else:
         
        dataf = datan[(datan>= datan[indices[0]+1]) & (datan < datan[indices[1]+1.0])]  


    phshift = min(dataf)-min(data)

    return phshift

def corr_strength(data, col, names = None, plot=False):

    """

    Determine the strength of each parameter with flux.

    @data - data (array)
    @col - data columns to be considered (Flux must be the first column). (list)
    @columns - columns dictionary for providing names to the correlations (optional). Otherwise just the columns will be printed. (dict)
    @return correlations and corresponding data column numbers in descending order
    """

    datan = data[:,col]
    datan = datan-np.mean(datan, axis=0)
    datan = datan/np.std(datan, axis=0) 
    cov_matrix = np.cov(datan.T)
    correlations = cov_matrix[0]
#    correlations = np.delete(correlations, 0)
#    names2 = np.copy(names)
    y = np.argsort(abs(correlations[1::]))[::-1]
    print y
    print col
    colmax = col[1::][y[0]]
    cormax = abs(correlations[1::][y[0]])
    y = np.argsort(abs(correlations))[::-1]
    for x in y:
        if x != 0:
            if names:
                print names[x], correlations[x]
            else:
                print 'Column '+str(x), correlations[x]
            if plot:
                plt.figure()
                plt.title(names[x])
                plt.plot(datan[:,x], datan[:,0], 'k.', ms=0.5)
    if plot:
        plt.show() 
#    print cov_matrix
    
    return colmax, cormax


def ind_split(time, gapsize=0.5):

    """
    Find indices on which data gap is larger than required.

    @time -must be in chronological order (array)
    @gapsize - minimum gap size required to split data

    """
    dataa = time[:-1]
    datab = time[1:]
    diff = datab-dataa
#    elements = np.linspace(0,len(diff)-1, len(diff)).astype(int) 
#    cuts = elements[diff>gap]
    cuts = np.where(diff>gapsize)[0]
    cuts = np.append(cuts, len(time)-1)
    return cuts


def array_split(data, indices):

    """
    Split mxn array at given indices in n.
    @data - mxn data array
    @indices- indices where data should be cut
 
    """
#    indices = np.append(indices, int(len(data)-1))
#    cuts = np.append(cuts, int(len(data)-1))
    datasn = []

    init = 0
    for x in indices:
        datan = data[init:x+1]
        datasn.append(datan)
        init = x+1

    return datasn


def plot_trends(datasets, cols, title=None):

    """
    plot each lc for each dataset in consecutive frames in top panel and
    the correlation of a given parameter with time in the bottom column.

    datasets - list of numpy arrays
    cols - cols in dataset being considered (time must be first, followed by flux, and then whatever parameter you would like)

    """

    nrows = 2
    ncols = int(len(datasets))
    
    f, axarr = plt.subplots(nrows, ncols, sharey = True)
    plt.locator_params(nbins=5)
    temp = np.vstack(datasets)
#    plt.locator_params(axis='y',nbins=3)
#    plt.locator_params(axis='x',nbins=3)
    if title:
        f.suptitle(title+' Correlation')
    for x in range(len(datasets)):
        datas = datasets[x]
        mint = min(datas[:,cols[0]]-datasets[0][:,cols[0]][0])
        maxt = max(datas[:,cols[0]]-datasets[0][:,cols[0]][0])
        ts = np.round((maxt-mint)/4., 0)
        print mint, maxt, ts
        if (maxt - mint) <= 2.0:
            ts = 0.5
        axarr[0,x].xaxis.set_ticks(np.arange(np.floor(mint)-1,np.ceil(maxt)+1,ts))
        axarr[0,x].plot(datas[:,cols[0]]-datasets[0][:,cols[0]][0], datas[:,cols[1]], 'k.')
#        axarr[1,x].set_ylim(minfl, maxfl)
#        axarr[1,x].plot(datas[:,cols[2]], datas[:,cols[1]],  'k.')
        minte =  min(datas[:,cols[2]])
        maxte =  max(datas[:,cols[2]])
        ts1 = np.round((maxte-minte)/4., 0)
        print maxte-minte
        if (maxte - minte) <= 1.0:
            minte=1
            maxte=0
            ts1=4

        axarr[1,x].xaxis.set_ticks(np.arange(np.floor(minte)-1,np.ceil(maxte)+1,ts1))
        axarr[1,x].plot(datas[:,cols[2]], datas[:,cols[1]],  'k.')
    plt.show()
    return


def param_detrend(time, flux, param, name, **kwargs):

    """
    Detrend a given parameter with flux. While the same fitting routines are accessed each time, this will provide fitting parameters e.g. knotspacing more appropriate to each correction.
    @flux - flux array
    @param - array of given parameter.
    @name - name of parameter: 'XCEN', 'YCEN', 'Phase', 'CCDT' or 'temp', 'Orbtemp' or orbtemp
    @return - flux correction
    """

    if name == 'XCEN' or name == 'YCEN':

        #correct for position
        try:
            tckPOScorrection = POSdetrendBRITE.detrendPOSITIONflux(time, param, flux,SPLINEknotpointsSPACING = SPLINEknotpointsSPACING,  **kwargs)

            correction = scInterp.splev(param, tckPOScorrection)
        except:
            #sometimes the knotpointspacing is larger than the lenght of the parameters
            SPLINEknotpointsSPACING = np.linspace((max(param)-min(param))/4.0, max(param)-min(param)-0.1, 4)   
            tckPOScorrection = POSdetrendBRITE.detrendPOSITIONflux(time, param, flux,SPLINEknotpointsSPACING = SPLINEknotpointsSPACING,  **kwargs)

            correction = scInterp.splev(param, tckPOScorrection)
   
    elif name == 'Phase' or name == 'phase':
        
        #correct for phase. This also has the requirement of being periodic so flux at 0 must equal flux at 1.
        # create realistic knot spacings. If the phase space is not covered well, then the defaults will not be adequate
        parmin = min(param)
        parmax = max(param)
        diff = parmax-parmin
        SPLINEknotpointsSPACING =np.array([diff/6., diff/5., diff/4., diff/3., diff/2., diff/1.5])
        #if there are big shifts then we need to shift in phase to avoid
        phshift = phase_shift(param)
        if phshift > 0.05:
            param = (param-phshift)%1
            parmin = min(param)
            parmax = max(param)
            diff = parmax-parmin
            SPLINEknotpointsSPACING =np.array([diff/4., diff/3., diff/2., diff/1.5])

            
            
        tckPhasecorrection = ORBdetrendBRITE.detrendORBITflux(time, param, flux,
SPLINEknotpointsSPACING=SPLINEknotpointsSPACING, **kwargs)
        correction = scInterp.splev(param, tckPhasecorrection)
        
        
    elif name == 'Orbtemp' or name == 'orbtemp':

        #this parameter has yet to be useful so I'm not sure what the best route to take here is.
        try:
            tckTEMPcorrection = TEMPdetrendBRITE.detrendTEMPflux(time, param, flux, **kwargs)
            correction = scInterp.splev(param, tckTEMPcorrection)
#        tckPOScorrection = POSdetrendBRITE.detrendPOSITIONflux(time, param, flux, **kwargs)

#        correction = scInterp.splev(param, tckPOScorrection) 
        except:
            #sometimes the knotpointspacing is larger than the lenght of the parameters
            SPLINEknotpointsSPACING = np.linspace((max(param)-min(param))/4.0, max(param)-min(param)-0.1, 4)   
            tckTEMPcorrection = TEMPdetrendBRITE.detrendTEMPflux(time, param, flux,SPLINEknotpointsSPACING = SPLINEknotpointsSPACING,  **kwargs)

            correction = scInterp.splev(param, tckTEMPcorrection)


    elif name == 'CCDT' or name == 'Temp' or name == 'temp':

        try:
            tckTEMPcorrection = TEMPdetrendBRITE.detrendTEMPflux(time, param, flux, **kwargs)
            correction = scInterp.splev(param, tckTEMPcorrection)
        except:
            SPLINEknotpointsSPACING = np.linspace((max(param)-min(param))/4.0, max(param)-min(param)-0.1, 4)   
            tckTEMPcorrection = TEMPdetrendBRITE.detrendTEMPflux(time, param, flux,SPLINEknotpointsSPACING = SPLINEknotpointsSPACING,  **kwargs)

            correction = scInterp.splev(param, tckTEMPcorrection)



    else:

        SPLINEknotpointsSPACING = np.linspace((max(param)-min(param))/4.0, max(param)-min(param), 4)   
        tckPOScorrection = POSdetrendBRITE.detrendPOSITIONflux(time, param, flux,SPLINEknotpointsSPACING = SPLINEknotpointsSPACING,  **kwargs)

#        raise ValueError(str(name)+' is not a valid option. you have either made an error or this parameter has not been implemented yet.')


    #TODO insert additional parameters for newer data releases that include chopping 

    return correction

def auto_detrend(data, cols, columns, mincorr=0.1, plot=False, **kwargs):

    """
    for a given set of data. Determine correlation, If correlation is over mincorr, correct and iterate until no correlations over mincorr exist. One correction maximum per parameter is done.
    @data - data array.
    @cols - columns to consider for correlations
    @columns - dictionary of columns
    @mincorr - minimum correlation necessary in order to detrend on a given parameter.
    """ 

    inv_columns = {v: k for k, v in columns.items()}
    maxcorr = mincorr
    names = [ inv_columns[k] for k in cols]
    time = data[:,columns['HJD']]
    flux = data[:,columns['FLUX']]
    coln = copy.copy(cols)
    datan = np.copy(data)
    if plot:
        plot_names = []
        plot_fluxes = []
        plot_params = []
        plot_corrections = []

    while (maxcorr >= mincorr) and (len(coln) > 1):

        colnum, maxcorr = corr_strength(datan, coln, names = names, plot=False)
        print colnum, maxcorr
        if maxcorr >= mincorr:
            param = data[:,colnum]
            name = inv_columns[colnum]
            print("Correcting for "+str(name))
            correction = param_detrend(time, flux, param, name, **kwargs)
            datan[:,columns['FLUX']] = flux-correction
            if plot:
                plot_names.append(name)
                plot_fluxes.append(flux)
                plot_params.append(param)
                plot_corrections.append(correction)
            
            
            flux = flux - correction
            coln.remove(colnum)
            names.remove(inv_columns[colnum])

    if plot == True:

        nrows = 2
        ncols = int(len(plot_names))
        print ncols
        if ncols == 0:
            print("I didn't correct anything so there is nothing to plot")
        elif ncols == 1:
            f, axarr = plt.subplots(nrows, ncols, sharey = True)
            x= 0
            axarr[0].set_title(plot_names[x])
            axarr[0].plot(plot_params[x], plot_fluxes[x]-np.mean(plot_fluxes[x]), 'k.', ms=1.0)
            axarr[0].plot(plot_params[x], plot_corrections[x], 'r.')
            axarr[1].plot(plot_params[x], plot_fluxes[x]-np.mean(plot_fluxes[x]) -plot_corrections[x], 'k.', ms=1.0)
            axarr[1].set_title('Residiuals')
            plt.show()

        else: 
            f, axarr = plt.subplots(nrows, ncols, sharey = True)
            for x in range(len(plot_names)):
                print "shape", axarr.shape
                print plot_names[x], x
                axarr[0,x].set_title(plot_names[x])
                print len(plot_params[x]), len(plot_fluxes[x]), np.mean(plot_fluxes[x]) 
                axarr[0,x].plot(plot_params[x], plot_fluxes[x]-np.mean(plot_fluxes[x]), 'k.', ms=1.0)
                axarr[0,x].plot(plot_params[x], plot_corrections[x], 'r.')
                axarr[1,x].plot(plot_params[x], plot_fluxes[x]-np.mean(plot_fluxes[x]) -plot_corrections[x], 'k.', ms=1.0)
                axarr[1,x].set_title('Residuals')
            plt.show()

    return flux

def orb_err(flux):

    """
    Calculate the RMS error of a set of flux points.
    @flux - array of fluxes
    @return - rms error
    
    """

    sig = abs(np.std(flux))
    errmean = sig/np.sqrt(len(flux))

    return errmean

def data_err(time, flux, oper):

    """
    Calculate the RMS error of an entire dataset.
    @time - time array
    @flux - flux array
    @oper - orbital period of telescope (for binning on the orbit)
    @return - RMS error of dataset (one value)

    """

    dmax = oper/(2*1440.)

    tstart = time[0]
    tfin = tstart+dmax 

    indices = ind_split(time, oper/2.0)
    datasplit = array_split(flux, indices)
    rms_error = []

    for x in datasplit:
        rms_error.append(orb_err(x))

    tot_err = np.sqrt(np.sum(np.array(rms_error)**2.0)/float(len(rms_error)))

    return tot_err


def bin_brite(data, tcol=0 , tel=None, per=None, error = False, return_dict =False):

    """
    
    Bin data on the  BRITE satellite orbital period. This will bin every quantity in the given array AND will provide an rms column for the flux if error = True.
    @data - data array. Flux is assumed to be the second column (array)
    @tcol - data column corresponding to the time array (integer)
    @tel - BRITE telescope ID. This is necessary unless period is provided (string)
    @per - orbital period...will accept telescope first if given (float)
    @error - calculate error in fluxes (boolean)
    @return_dict - when calculating errors this returns dictionary entries for added column
    @return - binned data array with flux errors if error == True and 2 separate dictionary entries if return_dict=True.
    
    """

    if len(data.shape) > 1:

        time = data[:,tcol]

    else:
        time = data[0]

    if tel:

        oper = {'UBr':100.3708/1440., 'BAb':100.3617/1440., 'BTr':98.2428/1440., 'BLb':99.6651/1440., 'BHr':97.0972/1440.} 
        per = oper[tel]

    if per == 'None':

        raise ValueError('you must provide a telescope ID OR your own period')


    datan = copy.copy(data)
    dmax = per/(2*1440.)

    tstart = time[0]
    tfin = tstart+dmax 

    indices = ind_split(time, per/2.0)
    datasplit = array_split(data, indices)
    binned = []
    if error:
        rms_error=[]
    for x in datasplit:
        
        binned.append(np.mean(x, axis=0))
        
        if error:
            
            #calculate rms error here
            rms_error.append(orb_err(x[:,1]))
            colnew = {}
            parnew = {}
    datab = np.vstack((binned))

    
    if error:
        rms_error = np.array(rms_error)
        
        #append rms error to datab array as last column
        datab = np.column_stack((datab, rms_error))
        colnew['Error'] = len(datab[0])-1
        parnew['column'+str(len(datab[0])-1)] = 'Error'

    if return_dict and error:
        return datab, colnew, parnew

    else:
        return datab   
    
def outlier_id(param, time, percent_low, percent_high, frac = 0.15, **kwargs):

    lowess = sm.nonparametric.lowess(param, time, frac=frac)
    time_l, param_l = lowess[:,0], lowess[:,1]
    param_corr = param-param_l
    toutliers, tmask = percentageclipBRITE.percentageFILTER(param_corr, percent_low, percent_high, full_output=True)

    return toutliers, tmask

def quick_cut(xcen, ycen,percent_low, percent_high):

    xoutliers, xmask =  percentageclipBRITE.percentageFILTER(xcen, percent_low, percent_high, full_output=True)
    youtliers, ymask =  percentageclipBRITE.percentageFILTER(xcen, percent_low, percent_high, full_output=True)

    xyoutliers = np.concatenate((xoutliers, youtliers))
    xyoutliers = set(list(xyoutliers))
    xyoutliers = sorted(xyoutliers)
    mask = np.full(len(xcen), True, dtype=bool)
    mask[xyoutliers]=False

    return mask

def orb_clip(data, per, iterative=True, sig = 3.5):

    """
    Bin on the orbit and do a sigma clip of each orbit.  
    
    data = data array. Time and Flux must be the first and second column respectively. (array)
    per = orbital period (float)
    iterative = do and iterative clip until, repeating until no values are above the give threshold
    sig = sigma threshold (float)
    
    """

    indices = ind_split(data[:,0], per/2.0)
    datasplit = array_split(data, indices)
    datasplitNEW = []
    
    for x in datasplit:
        threshold = sig*np.std(x[:,1])
        absflux = np.abs(x[:,1]-np.mean(x[:,1]))
        xnew = x[absflux < threshold]

        if iterative == True:
            threshold = sig*np.std(xnew[:,1])
            absflux = np.abs(xnew[:,1]-np.mean(xnew[:,1]))

            while np.any(absflux > threshold):
                
                xnew = xnew[absflux < threshold]
                threshold = sig*np.std(xnew[:,1])
                absflux = np.abs(xnew[:,1]-np.mean(xnew[:,1]))
        
        datasplitNEW.append(xnew)

    dataNEW = np.vstack((datasplitNEW))

    return dataNEW

def orb_cut(data, orbper, sig = 4.0, minpoints=4.0):

    """
    Cut the orbits which have an error larger than the threshold 

    @data - data array: time and flux must be the 1st and second column respectively (2D array)
    @sig - sigma threshold
    @minpoints - minimum number of points acceptable in a given orbit.
    @return - filtered data (2D array)
    """

    indices = ind_split(data[:,0], orbper/2.0)
    datasplit = array_split(data, indices)
    datasplitn = []
    sigs = []
    for x in range(len(datasplit)):
        if len(datasplit[x]) > 4:
            datasplitn.append(datasplit[x])
            sigs.append(np.std(datasplit[x][:,1]))
        

    threshold = np.mean(sigs)+sig*np.std(sigs)
    results = sigs < threshold

    while np.any(results == False):
        datasplitn = [ item for item, flag in zip( datasplitn, results ) if flag == True ]
        sigs = sigs < threshold
        results = sigs < threshold

    datan = np.vstack((datasplitn))

    return datan

def comb_datasets(tel, root, pathIN, obsid=None):

    """
    Combine all Brite datasets of a given telescope and obsid

    @tel - telescope id (string)
    @obsid - observation id (string)
    @pathIN -path to directory where files are kept (str)
    @filt - combine all files of a given filter. two valid choices: r, b. will override telescope if selected
    @return None. A file is created in the pathIN directory of the combined data
    """


    if not pathIN[-1] == '/':
        pathIN = pathIN+'/'
    search_string = pathIN
    if obsid:
        search_string = search_string+'*'+str(obsid)
    if tel == 'red' or tel == 'blue':

        search_string = search_string+'*'+str(tel[0])+'_'
#        comb_list = glob(pathIN+'*'+str(obsid)+'*'+str(tel)+'*'+str(root))
    else:
        search_string = search_string+'*'+str(tel)
#        comb_list = glob(pathIN+'*'+str(obsid)+'*'+str(filt)+'_*'+str(root))
    search_string = search_string+'*'+str(root)
    comb_list = glob(search_string)
    datas = []
    
    for x in comb_list:

        data = np.loadtxt(x, comments='c')
        datas.append(data)

    
    datan = np.vstack((datas))
    y = np.argsort(datan[:,0])
    datan = datan[y]
    targ = comb_list[0].split('/')[-1].split('_')[0]

    if not obsid:
        obsid = 'all'
    root = root.split('.')[0]
#    if 'b.' in  comb_list[0]:
    filename = str(pathIN)+str(targ)+'_'+str(obsid)+'_'+str(tel)+'_'+str(root)+'_comb.dat'
    np.savetxt(filename, datan)

    return

def duty_cycle(time, tel=None, orbper=None):

    """
    Calculate the duty cycle for a given set of BRITE data. Specifically, this calculates the number of orbits where data with respect to the number of orbits available.
    
    @time = array of observation times (array)
    @orbper = orbital period of satellite (float)

    """
    if tel:
        oper = {'UBr':100.3708/1440., 'BAb':100.3617/1440., 'BTr':98.2428/1440., 'BLb':99.6651/1440., 'BHr':97.0972/1440.} 
        orbper = oper[tel]
    
    obs_poss = int((max(time)-min(time))/orbper)
    #fake binning out because it requires two columns
    data = np.column_stack((time,time))
    binned = bin_brite(data, tcol=0 , per=orbper, error = False, return_dict =False)  

    num_obs = len(binned)

    duty = (num_obs/float(obs_poss))*100

    return duty
