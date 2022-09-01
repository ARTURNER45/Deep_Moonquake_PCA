import obspy
import waveformdb
from waveformtable import Waveform
import responsecorr
import seisproc
import numpy as np
import os,glob
import general
import datetime

#def read_single_moonquakes(catalog='moonquake_trial',stats=['s15','s16'],chans=['MH1','MH2'],nsts = ['A01'] ):
def read_single_moonquakes(lis):
    fil =open(lis,"r")  
    st = obspy.Stream()
    lineSTORE = [] 
    for line in fil: 
        #print(line) 
        #print(line.split(',')[0]) 
        #print(line.split(',')[1]) 
        #dire='/home/users/lina3509/Documents/Planets/moonquake_trial/A01/S12/'
        #a = line.split(',')[0] 
        sti = obspy.read(line.split(',')[0]) 
        #print(line.split(',')[0],line.split(',')[1])
        sti[0].stats.picks = obspy.core.AttribDict() 
        sti[0].stats.picks.cor1 = line.split(',')[1]
        sti[0].stats.picks.cor1 = int(float(line.split(',')[1][:-1]))
        lineSTORE.append(line.split('\n')[0])
        st=st+sti
    return st
#,lineSTORE
    # note the stations and channels
    #if isinstance(chans,str):
        #chans=[chans]
    #chans=np.array([chan.upper() for chan in chans])
    #if isinstance(stats,str):
        #stats=[stats]
    #stats=np.array([stat.lower() for s\tat in stats])
    #if isinstance(nsts,str):
        #stats=[nsts]

    # directory of interest
    #fdir='/home/users/lina3509/Documents/Planets/S12_good'
    #if catalog in ['moonquake_trial']:
        #fdir2='moonquake_trial'
    #fdir=os.path.join(fdir,fdir2)
    
    # read the data    

    # This is fine - it reads everything into the same string. 

    


    #fnames=glob.glob(os.path.join(fdir,'A01','*','*.spike_cheader.SAC'))
    #print(fnames)
    
    #for fname in fnames:
        #sti=obspy.read(fname)
        #st=st+sti
     #select some
    #if minqual>0:
        #st=select_stacks(st,catalog=catalog)

    # note a pick time # this is just an estiamte from the catalog
    #d = datetime.timedelta(minutes = 30 ) #put this as the shift * delta as this is where the cross correlation peak is. 
    #tpick=obspy.UTCDateTime(1970,1,1,1,2,0)

    #for tr in st:
        #tr.stats.t3 = tr.stats.starttime + d 




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
        #print(line) 
        #print(line.split(',')[0]) 
        #print(line.split(',')[1]) 
        #dire='/home/users/lina3509/Documents/Planets/moonquake_trial/A01/S12/'
        #a = line.split(',')[0] 
        sti = obspy.read(line.split('\n')[0]) 
        #print(line.split(',')[0],line.split(',')[1])
        #sti[0].stats.picks = obspy.core.AttribDict() 
        #sti[0].stats.picks.cor1 = line.split(',')[1]
        #sti[0].stats.picks.cor1 = int(float(line.split(',')[1][:-1])) * sti[0].stats.delta
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
        iok[k]=np.diff(general.minmax(tr.data))>0

        # find maximum
        if iok[k]:
            find_maxenv(tr,pk='t4',flm=[0.75,1.75],f_low=1/60,twin=twin)
            tmax[k]=tr.stats.t4

            # also set a tentative pick time for now
            tspick[k]=tmax[k]-100

    # directory of interest
    #fdir='/home/users/lina3509/Documents/Planets'
    #if catalog in ['moonquake_trial']:
        #f#dir2='moonquake_trial'
    #fname=os.path.join(fdir,fdir2,'stack_quality')

    # write to file
    fname = '/home/users/lina3509/Documents/Planets/S16_peaked_pickes_MHN.txt'
    fl=open(fname,'w')
    for k in range(0,len(st)):
        tr=st[k]
        #print(lineSTORE[k])
        fl.write(str(lineSTORE[k])+","+str(tspick[k])+"\n")
    fl.close()

def read_stacks(catalog = 'Bulow' ,stats=['s12','s14','s15','s16'],chans=['*']):
                #minqual=1): #I have changed minqual
    """
    :param     catalog: which set of moonquakes stacks to grab
                         (default: 'Bulow')
    :param       stats: which stations to get
    :param       chans: which channels to get
    :param     minqual: only keep those stacks with some noted quality value
    :return         st: an obspy stream with one waveform per stack
    """

    # note the stations and channels
    if isinstance(chans,str):
        chans=[chans]
    chans=np.array([chan.upper() for chan in chans])
    if isinstance(stats,str):
        stats=[stats]
    stats=np.array([stat.lower() for stat in stats])

    # directory of interest
    fdir='/home/users/lina3509/Documents/Apollo'
    if catalog in ['bulow','Bulow']:
        fdir2='Bulow_et_al_2007_deep_moonquake_stacks'
    fdir=os.path.join(fdir,fdir2)
    
    # which nests are available
    nsts=glob.glob(os.path.join(fdir,'A*'))


    # read the data    
    st = obspy.Stream()
    for nst in nsts:
        for stn in stats:
            for chn in chans:
                fnames=glob.glob(os.path.join(fdir,nst,stn,'*.'+chn))
                for fname in fnames:
                    sti=obspy.read(fname)
                    st=st+sti

     #select some
    #if minqual>0:
        #st=select_stacks(st,minqual=minqual,catalog=catalog)

    # note a pick time
    tpick=obspy.UTCDateTime(1970,1,1,1,2,0)

    for tr in st:
        tr.stats.t3=tpick-tr.stats.starttime


    return st

def write_stack_info(catalog='moonquake_trial'):
    """
    :param     catalog: which set of moonquakes stacks to grab
                         (default: 'Bulow')
    """

    # read the data
    #st=read_stacks(catalog=catalog)
    st = read_single_moonquakes(lis='/home/users/lina3509/Documents/Apollo/A01_moonquakes/lists/S12_MHN_peaked.txt')  

    # time window that should include the maximum
    # this needs changing for the moonquake dates. 
    #twin=[obspy.UTCDateTime(1970,1,1,0,59),
          #obspy.UTCDateTime(1970,1,1,1,9)]
    
    # check that there's data for each
    iok=np.zeros(len(st),dtype=bool)
    #tmax=[obspy.UTCDateTime(1970,1,1,1)]*len(iok)
    #tspick=[obspy.UTCDateTime(1970,1,1,1)]*len(iok)
    for k in range(0,len(st)):
        tr=st[k]
        # check that they're not all zero
        iok[k]=np.diff(general.minmax(tr.data))>0

        # find maximum
        if iok[k]:
            #twin = [tr.stats.starttime + 10, tr.stats.starttime + 17000]
            find_maxenv(tr,pk='t4',flm=[0.75,1.75],f_low=1/60,twin=twin)
            tmax[k]=tr.stats.starttime+tr.stats.t4

            # also set a tentative pick time for now
            tspick[k]=tmax[k]-100

    # directory of interest
    fdir='/home/users/lina3509/Documents/Apollo/A01_moonquakes/'
    #if catalog in ['moonquake_trial']:
        #f#dir2='moonquake_trial'
    #fdir=os.path.join(fdir,fdir2)

    # write to file
    fname = '/home/users/lina3509/Documents/Apollo/A01_moonquakes/S12_peak_tent_pick_test.txt'
    fl=open(fname,'w')
    for k in range(0,len(st)):
        tr=st[k]
        tstr=tmax[k].strftime('%Y-%b-%d-%H-%M-%S')
        tstrp=tspick[k].strftime('%Y-%b-%d-%H-%M-%S')
        fl.write(tr.get_id()+',{:d},{:s},{:s}\n'.format(iok[k],tstr,tstrp))
    fl.close()


def read_corr_picks(lis):
	fil =open(lis,"r")  

	for line in fil: 
		#print(line) 
		#print(line.split(',')[0]) 
		dire = '/home/users/lina3509/Documents/Planets/s12_good_smaller/' 
		st = obspy.read(dire + line.split(',')[0]) 
		st[0].stats.picks = obspy.core.AttribDict() 
		st[0].stats.picks.cor1 = int(line.split(',')[1]) * st[0].stats.delta
		#print(st[0].stats)
		st[0].write(dire + line.split(',')[0])


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
    fdir='/home/users/lina3509/Documents/Planets'
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
    fdir='/home/users/lina3509/Documents/Planets'
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
    fdir='/home/users/lina3509/Documents/Planets'
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
    

def read_data(t1=None,t2=None,remresp=True):
    """
    :param         t1: start time (default: 2019-jun-10)ds.wr
    :param         t2: end time (default: t1+1 hour)
    :param    remresp: go ahead and remove the response
    :return        st: a Stream of waveforms
    """

    # time ranges
    if t1 is None:
        t1=obspy.UTCDateTime(2019,6,10)
    if t2 is None:
        t2=t1+3600.

    if remresp:
        # read buffered data
        bfr=300.
        st=read_data(t1=t1-bfr,t2=t2+bfr,remresp=False)
        
        # filter ranges
        fmin=1./np.minimum((t2-t1)/2.,500)
        fmax=st[0].stats.sampling_rate
        pre_filt=[fmin/2.,fmin,fmax*0.8,fmax*0.9]

        # taper and get rid of blanks
        st=st.trim(starttime=t1-bfr,endtime=t2+bfr,pad=True)
        msk=seisproc.prepfiltmask(st)
        st.taper(max_length=bfr,max_percentage=0.5)

        # remove response
        responsecorr.removeresponse(st,pre_filt=pre_filt)
        seisproc.addfiltmask(st,msk)

        # to desired time range
        st.trim(starttime=t1,endtime=t2)


    else:
        
        # read the data
        session=waveformdb.opendatabase('mars')
        q=session.query(Waveform)
        session.close()
        
        # just some time range
        qq=q.filter(Waveform.starttime<=t2,Waveform.endtime>=t1)    
        
        # read these values
        st=obspy.Stream()
        for wv in qq:
            st=st+wv.waveform(t1=t1,t2=t2)
            
        # merge if necessary
        st=st.merge()

    return st

