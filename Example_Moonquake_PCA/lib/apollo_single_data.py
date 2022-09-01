#File names and paths are required for pre-processing functions. 

import obspy
import scipy
import numpy as np
import os,glob
import datetime
from math import log10,floor

def read_single_moonquakes(lis):
    fil =open(lis,"r")  
    st = obspy.Stream()
    lineSTORE = [] 
    for line in fil: 
        sti = obspy.read(line.split(',')[0]) 
        #print(line.split(',')[0],line.split(',')[1])
        sti[0].stats.picks = obspy.core.AttribDict() 
        sti[0].stats.picks.cor1 = line.split(',')[1]
        sti[0].stats.picks.cor1 = int(float(line.split(',')[1][:-1]))
        lineSTORE.append(line.split('\n')[0])
        st=st+sti
    return st


def write_single_moonquakes_info(lis):
    """
    :param     catalog: which set of moonquakes stacks to grab
                         (default: 'Bulow')
    """

    # read the data
    fil =open(lis,"r")  
    st = obspy.Stream()
    lineSTORE = [] 
    for line in fil: 
        sti = obspy.read(line.split('\n')[0]) 
        lineSTORE.append(line.split('\n')[0])
        st=st+sti

    #d = 

    # time window that should include the maximum
    twin=[ datetime.timedelta(minutes=20),datetime.timedelta(minutes=40)]
    
    # check that there's data for each
    iok=np.zeros(len(st),dtype=bool)

    tmax=[obspy.UTCDateTime(1970,1,1,1)]*len(iok)
    tspick=[obspy.UTCDateTime(1970,1,1,1)]*len(iok)
    for k in range(0,len(st)):
        tr=st[k]
        # check that they're not all zero
        iok[k]=np.diff(minmax(tr.data))>0

        # find maximum
        if iok[k]:
            find_maxenv(tr,pk='t4',flm=[0.75,1.75],f_low=1/60,twin=twin)
            tmax[k]=tr.stats.t4

            # also set a tentative pick time for now
            tspick[k]=tmax[k]-100


    # write to file
    fname = '/home/users/lina3509/Documents/Planets/S16_peaked_pickes_MHN.txt'
    fl=open(fname,'w')
    for k in range(0,len(st)):
        tr=st[k]
        #print(lineSTORE[k])
        fl.write(str(lineSTORE[k])+","+str(tspick[k])+"\n")
    fl.close()



def write_stack_info(catalog='moonquake_trial'):
    """
    :param     catalog: which set of moonquakes stacks to grab
                         (default: 'Bulow')
    """

    # read the data
    #st=read_stacks(catalog=catalog)
    st = read_single_moonquakes(lis='')  

    for k in range(0,len(st)):
        tr=st[k]
        # check that they're not all zero
        iok[k]=np.diff(minmax(tr.data))>0

        # find maximum
        if iok[k]:
            #twin = [tr.stats.starttime + 10, tr.stats.starttime + 17000]
            find_maxenv(tr,pk='t4',flm=[0.75,1.75],f_low=1/60,twin=twin)
            tmax[k]=tr.stats.starttime+tr.stats.t4

            # also set a tentative pick time for now
            tspick[k]=tmax[k]-100

    # directory of interest
    fdir='ADD PATH TO INPUT LIST HERE'
    #if catalog in ['moonquake_trial']:
        #f#dir2='moonquake_trial'
    #fdir=os.path.join(fdir,fdir2)

    # write to file
    fname = 'ADD PATH TO OUTPUT FILE HERE'
    fl=open(fname,'w')
    for k in range(0,len(st)):
        tr=st[k]
        tstr=tmax[k].strftime('%Y-%b-%d-%H-%M-%S')
        tstrp=tspick[k].strftime('%Y-%b-%d-%H-%M-%S')
        fl.write(tr.get_id()+',{:d},{:s},{:s}\n'.format(iok[k],tstr,tstrp))
    fl.close()



def find_maxenv(tr,pk='t4',flm=[0.25,2],f_low=1/60,twin=None):
    """
    find the location with the max amplitude and mark it
    :param           tr: the waveforms or traces
    :param           pk: the pick to mark
    :param          flm: frequency band to filter to first
    :param         twin: a time window to allow (default: None, ignored)
    """

    if isinstance(tr,obspy.Stream):
        for tri in tr:
            find_maxenv(tri,pk=pk,flm=flm,twin=twin)
    else:
        tri=tr.copy()
        tri.filter('bandpass',freqmin=flm[0],freqmax=flm[1],zerophase=True)
        if twin is not None:
            tri.trim(starttime = tr.stats.starttime + twin[0],endtime= tr.stats.starttime + twin[1])
        tre=tri.copy()
        env=obspy.signal.filter.envelope(tre.data)
        tre.data=env
        tre.taper(side='both',max_percentage=0.2,max_length=1/f_low*2)
        tre.filter('lowpass',freq=f_low,zerophase=True)
        imax=np.argmax(tre.data)
        tmax=tre.stats.starttime+imax*tr.stats.delta
        #print(tmax)
        
        tr.stats[pk]=tmax-tr.stats.starttime

def read_stack_picks(st,catalog='moonquake_trial',pk='t3',remnopick=True,bystation=True):

    """
    :param          st: a set of stacks 
    #:param     minqual: minimum quality to accept
    :param     catalog: which set of stacks to consider
    :param   remnopick: remove events that don't have a pick
    :return        sts: a subset of the stacks, as selected
    """

    # directory of interest
    fdir='ADDD PATH TO DIRECTORY WITH SCRIPTS'
    if catalog in ['moonquake_trial']:
        fdir2='moonquake_trial'
    fdir=os.path.join(fdir,fdir2)
    fname=os.path.join(fdir,'saved_picks')

    # read the data
    vls=np.loadtxt(fname,delimiter=',',dtype=str)
    ids=vls[:,0]
    tmax=[obspy.UTCDateTime(datetime.datetime.strptime(tm,'%Y-%b-%d-%H-%M-%S'))
          for tm in vls[:,1]]
    tmax=dict([(ids[k],tmax[k]) for k in range(0,len(ids))])


    # do by pick
    if bystation:
        tmaxv=np.array(list(tmax.values()))
        tmaxk=np.array(list(tmax.keys()))
        sevk=np.array(['.'.join(ky.split('.')[1:3]) for ky in tmaxk])
        sev,ix=np.unique(sevk,return_inverse=True)
        tmaxs={}
        for k in range(0,sev.shape[0]):
            vls=tmaxv[ix==k]
            vls=vls[0]+np.median(vls-vls[0])
            tmaxs[sev[k]]=vls
            for ky in tmaxk[ix==k]:
                tmax[ky]=vls

    # assign values
    for tr in st:
        if pk in tr.stats.__dict__.keys():
            tpk=tmax.get(tr.get_id(),obspy.UTCDateTime(1900,1,1))
            if tpk>obspy.UTCDateTime(1950,1,1):
                tr.stats[pk]=tpk-tr.stats.starttime
            else:
                if remnopick:
                    st.remove(tr)

            
def read_stack_tmax(st,catalog='moonquake_trial',pk='t4'):

    """
    :param          st: a set of stacks 
    #:param     minqual: minimum quality to accept
    :param     catalog: which set of stacks to consider
    :return        sts: a subset of the stacks, as selected
    """

    # directory of interest
    fdir='ADD PATH TO DIRECTORY WITH SCRIPTS'
    if catalog in ['moonquake_trial']:
        fdir2='moonquake_trial'
    fdir=os.path.join(fdir,fdir2)
    fname=os.path.join(fdir,'stack_quality')

    # read the data
    vls=np.loadtxt(fname,delimiter=',',dtype=str)
    ids=vls[:,0]
    tmax=[obspy.UTCDateTime(datetime.datetime.strptime(tm,'%Y-%b-%d-%H-%M-%S'))
          for tm in vls[:,2]]
    tmax = dict(zip(ids, tmax))
    #tmax=dict([(ids[k],tmax[k]) for k in range(0,len(ids))])

    # assign values
    for tr in st:
        tr.stats[pk]=tmax[tr.get_id()]-tr.stats.starttime

def select_stacks(st,catalog='moonquake_trial'):
    """
    :param          st: a set of stacks 
    #:param     minqual: minimum quality to accept
    :param     catalog: which set of stacks to consider
    :return        sts: a subset of the stacks, as selected
    """

    # directory of interest
    fdir='ADD PATH TO DIRECTORY WITH SCRIPTS'
    if catalog in ['moonquake_trial']:
        fdir2='moonquake_trial'
    fdir=os.path.join(fdir,fdir2)
    fname=os.path.join(fdir,'stack_quality')

    # read the data
    vls=np.loadtxt(fname,delimiter=',',dtype=str)
    ids=vls[:,0]
    quals=vls[:,1].astype(int)
    dqual=dict([(ids[k],quals[k]) for k in range(0,len(ids))])

    # a new stream
    sts=obspy.Stream()
    for tr in st:
        #if dqual[tr.get_id()]>=minqual:
        sts.append(tr)

    return sts
    
def minmax(x,bfr=1.):
    """
    :param      x:   set of values
    :param    bfr:   how much to multiply the limits by (default: 1.)
    :return   lms:   limits
    """

    # minmax
    lms = np.array([np.min(x),np.max(x)])

    if bfr!=1.:
        lms = np.mean(lms)+np.diff(lms)[0]*bfr*np.array([-.5,.5])

    return lms
    
 
def closest(xvals,x):
    """
    :param      xvals:    set of sorted values
    :param          x:    values of interest
    :return        ix:    index of closest value
    """

    xvals = np.atleast_1d(xvals)
    x = np.atleast_1d(x)

    # index before and after
    ix = np.searchsorted(xvals,x,'left')

    # in range
    ix = np.maximum(ix,1)
    ix = np.minimum(ix,len(xvals)-1)
    
    # before or after?
    dx1 = np.abs(xvals[ix-1]-x)
    dx2 = np.abs(xvals[ix]-x)
    dx1 = dx1<dx2

    ix[dx1] = ix[dx1]-1
    ix = np.maximum(ix,0)

    return ix
    
    
def prepfiltmask(st,tmask=3.):
    """
    :param        st: waveforms or trace
    :param     tmask: time window to mask within some interval or endpoint, 
                         in seconds
    :return      msk: a set of waveforms with masks
    """

    if isinstance(st,obspy.Trace):

        # allowable values
        try:
            ms=np.logical_or(st.data.mask,np.isnan(st.data.data))
            data=st.data.data
        except:
            ms = np.isnan(st.data)
            data=st.data

        # times
        tm=np.arange(0.,data.size)

        # interpolate
        if np.sum(~ms):
            data[ms]=scipy.interp(tm[ms],tm[~ms],data[~ms])
        else:
            data[:]=0.
        
        # and copy data
        st.data = data    

        # to mask
        nwin = int(tmask/st.stats.delta)
        nwin = np.maximum(nwin,1)
        win = scipy.signal.boxcar(nwin*2+1)
        ms = ms.astype(float)
        ms = scipy.signal.convolve(ms,win,mode='same')

        # also the beginning and end
        if tmask != 0.:
            ms[0:nwin+1]=1.
            ms[-nwin:]=1.

        # place in trace
        ms = np.minimum(ms,1.)
        msk = st.copy()
        msk.data = ms

    elif isinstance(st,obspy.Stream):
        msk = obspy.Stream()

        for tr in st:
            mski = prepfiltmask(tr,tmask=tmask)
            msk.append(mski)

    return msk
    
def addfiltmask(st,msk):
    """
    :param        st: waveforms or trace
    :param       msk: a set of waveforms with masks
    """
    
    if isinstance(st,obspy.Trace):
        # the mask
        #ms = msk.data.astype(bool)
        ms = msk.data > 0.1
        try:
            # if there's already a mask, combine them
            st.data.mask=np.logical_or(st.data.mask,ms)
        except:
            st.data=np.ma.masked_array(st.data,mask=ms)

    elif isinstance(st,obspy.Stream):

        for tr in st:
            # select mask
            ms = msk.select(id=tr.id)[0]

            # add filter
            addfiltmask(tr,ms) 
    
