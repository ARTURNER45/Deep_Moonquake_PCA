import obspy
import waveformdb
from waveformtable import Waveform
import responsecorr
import seisproc
import numpy as np


def read_data(t1=None,t2=None,remresp=True):
    """
    :param         t1: start time (default: 2019-jun-10)
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

