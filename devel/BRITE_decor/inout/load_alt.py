"""
Routines to load various files from BRITE reduction

Last update: 23 June 2016

@author: Bert Pablo

"""


import numpy as np

def load_header(filename, comment='c'):

#    f=open(filename)
    header = []
    with open(filename) as fp:
        for line in fp:

            if line[0] == comment:
                header.append(line)
            else:
                break 

    return header           

"""
get_params() Return a dictionary with  all the parameters and their values for a BRITE header file.


"""

def get_params(header):

    # pull out all lines with a parameter value
    
    headern = [i for i in header if '=' in i]
    params = {}
    for x in range(len(headern)):

        #remove unwanted info from string
        stuff = headern[x].strip('c').strip().split('/')[0].split('=')
        # load into dictionary
        params[stuff[0].strip()] = stuff[1].strip()


    return params

def create_filename(targ, tel, obnum, obsid, setup, release, typ='stare', position=None):

    """
    create_filename: Use all the parameters of your data set to create the filename for the BRITE input file. This is purely a convenience function.

    """

    #define filename
    if typ == 'stare':
        base = '_APa2s5_R'+str(release)+'.dat'
    elif typ == 'chop':
        if position == None:
            base = '_APa3s2chop_R'+str(release)+'.dat'
        else:
            base = '_APa3s2chop_R'+str(release)+'_pos'+str(position)+'.dat'

    filename = str(targ)+'_'+str(obnum)+'_'+str(obsid)+'_'+str(tel)+'_setup'+str(setup)+str(base)
    
    
    return filename


def load_dataset(filename, pathIN=None):

    """
    Given a handful of parameters about the data file, load the file as 2D array, along with a dictionary which specifies the location of each column, and another dictionary which gives all the parameters in the file header.


    """



    # define period for making phase column
    per = {'UBr':100.3708/1440., 'BAb':100.3617/1440., 'BTr':98.2428/1440., 'BLb':99.6651/1440., 'BHr':97.0972/1440.}

#    oper = per[tel]
    
    #if a filename is not defined retrieve one

    


    if pathIN is not None:
        if not(pathIN[-1] == '/'):
            pathIN += '/'
        filename = pathIN+filename   

    #load data default comments for data = c

    data = np.loadtxt(filename, comments='c')  
    # remove NaNs

    #add later if necessary

    #get header and parameter values

    header = load_header(filename)
    params = get_params(header)
    
    #create extra columns with information on important values

    stacks = np.ones(len(data[:,0]))*np.int(params['ObsStack']) #number of stacks
    exptime = np.ones(len(data[:,0]))*np.float(params['ObsExpoT'])/1000 #exposure time divided by 1000 because it's given in milliseconds


    # add orbital period to header if not already there
    try:
        params['orbper']

    except:
        tel = params['SatellID']
        orbper = per[tel]
        params['orbper'] = orbper


    # Now create a dictionary which gives the column name and location

    columns = {}
    a = 0
    for x in params:
        if 'column' in x:
            
            col = params[x] 
            columns[col] = int(x[-1])-1
            a=a+1
    try:
        columns['stack']
    except:
        columns['stack'] = a
        data = np.column_stack((data, stacks))
        params['column'+str(a+1)] = 'stack'
    try:
        columns['exptime']
    except:
        columns['exptime'] = a+1
        data = np.column_stack((data, exptime))
        params['column'+str(a+2)] = 'exptime'

    return data, columns, params

"""
make a header file using a dictionary of keyword values. This will follow the same basic format as the original BRITE ascii, files, except the definitions of each keyword are not kept

"""


def make_header(params, comments='c'):

    header = ' start header ------------------------------------------------------------------------------------\n'
    for x in params:
         header+=' '+str(x)+'= '+str(params[x])+'\n'
    header+=' end header ---------------------------------------------------------------------------------------'
    return header

def save_dataset(filename, data, columns=None, params=None, pathOUT=None):
    header = None
    c={}
    r={}
    if pathOUT is not None:
        if not(pathOUT[-1] == '/'):
            pathOUT += '/'
        filename = pathOUT+filename   

    if columns != None:

        try:

            data = np.delete(data, columns['FLAG'], 1)
            r = dict(params)
            c = dict(columns)
            del r['column'+str(int(columns['FLAG'])+1)]
            del c['FLAG']
            if params != None:
                
                #get rid of all the columns
                for x in r.keys():
                    if 'column' in x:
                        del r[x]

                #add back all the columns that still exist
                for x in c:
                    if c[x] < columns['FLAG']:
                        entry = 'column'+str(c[x]+1)
                    else:
                        entry = 'column'+str(c[x])
                    r[entry] = x
                    
                r['OFilName'] = filename

                header = make_header(r)

        except:
            
            print('There is no column "FLAG" to delete')



    #make header if there is info to do it

    if params != None:

        if not header: 
            r = dict(params)
            r['OFilName'] = filename
            header = make_header(r)
            
        np.savetxt(filename, data, header=header, comments='c')
    
    else:
        np.savetxt(filename, data)

    return
