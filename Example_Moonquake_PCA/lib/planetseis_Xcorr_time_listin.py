import obspy
import numpy as np
import apollo_single_data as adat
import matplotlib.pyplot as plt
from scipy import optimize
from scipy import signal
from scipy import stats
from xcorr import correlate_maxlag, correlate_template, get_lags
from matplotlib import gridspec


class dstack:
    def __init__(self,dsi=None):

        # assign the frequency bands to consider
        self.pickflms(flms=1)

        # copy from another object if necessary
        if dsi is not None:
            for nm in dsi.__dict__.keys():
                self.__setattr__(nm,dsi.__getattribute__(nm))

    def setup(self,lis='S16_MHN_flat_lowerthresh_new.txt',dtype='single_moonquake',time_window = [-300,300]):
        """
        :param        synth: make synthetic data
        """
        
        if dtype in ['synth']:
            self.makesynthnoise()
            self.pick_time_window([50,250])
            self.initpars()
            self.addsynthsignal()
        elif dtype in ['moonstacks']:
            self.read_moon_stacks()
            self.pick_time_window(time_window)
            self.initpars()
            self.pickflms(flms=[0.25,1.75])
        elif dtype in ['single_moonquake']:
            self.read_single_moonquakes_2(lis)
            self.pick_time_window(time_window)
            #self.pick_time_window([750,1750])
            self.initpars()
            self.pickflms(flms=[0.25,1.75])            

        # get filtering and initial extraction set up
        self.filterdata()
        self.initadata()

        # and initialize the search
        if dtype in ['synth']:
            self.initsearch(method = 'iterative_grid',tshfmax=60)
        elif dtype in ['moonstacks']:
            self.initsearch(tshfmax=20)
        elif dtype in ['single_moonquake']:
            self.initsearch(method = 'iterative_grid',tshfmax=60)

    def pickflms(self,flms=8):
        """
        pick some frequency bands, in Hz
        :param        flms: either 
                   an integer specifying the number of bands or
                   the desired frequency bands or
                   the frequency band boundaries
        """

        flms=np.atleast_1d(flms)
        if flms.size==1:
            self.flms=np.logspace(-1,0,int(flms)+1)
            self.flms=[self.flms[k:k+2] for k in range(0,self.flms.size-1)]
        elif flms.ndim==1:
            self.flms=[flms[k:k+2] for k in range(0,flms.size-1)]

        self.flms=np.atleast_2d(np.array(self.flms))
        self.Nfreq=self.flms.shape[0]

        
    def read_single_moonquakes_2(self,lis,normdata=True,tgrab=[-300,600]):
        st=adat.read_single_moonquakes(lis=lis)
        iev = np.arange(0,len(st))
        print(iev)
        stns,istn=np.unique([tr.stats.station for tr in st],return_inverse=True)
        chns,ichn=np.unique([tr.stats.channel for tr in st],return_inverse=True)
        self.stations=stns
        self.channels=chns
        self.Nev = len(st)

        # station info
        self.Nstat=len(self.stations)
        self.Nchan=len(self.channels)


        # length of window to read
        ndat=int(np.diff(tgrab)[0]/st[0].stats.delta)
        self.ndata=ndat
        self.data=np.ma.masked_array(np.ndarray([ndat,self.Nev,self.Nchan,self.Nstat]),
                                     mask=True)


        # grab out some portion
        for k in range(0,len(st)):
            tr=st[k].copy()
            #print(iev.shape)
            #tr.plot()
            tref=tr.stats.starttime + tr.stats.picks.cor1
            tr.trim(starttime=tref+tgrab[0], 
            endtime=tref+tgrab[1]+3*tr.stats.delta,pad=True)
            data=tr.data[0:ndat]
            #plt.plot(data)
            # normalise if desired
            if normdata:
                data=data/np.max(np.abs(data))
            self.data.mask[:,k,ichn[k],istn[k]]=False
            self.data[:,iev[k],ichn[k],istn[k]]=data


        # the timing
        self.delta= st[0].stats.delta
        self.tdata=np.arange(0,ndat)*self.delta
        self.dsamp=1./self.delta
        self.pick_time_window([-200,0])

        # only want data with all channels
        self.ident_okay_events(rembad=False)

    def makesynthnoise(self):
        """
        make some fake data to work with
        """

        # to mimic Mars data
        self.delta=1./20
        self.dsamp=1./self.delta
        self.Nstat=1
        self.Nchan=1
        self.Nev = 10
        self.ndata=int(600/self.delta)

        # data and timing
        self.data=np.random.randn(self.ndata,self.Nev,self.Nchan,self.Nstat)
        #self.data = np.ones((self.ndata,self.Nev,self.Nchan,self.Nstat))
        self.tdata=np.arange(0,self.ndata)*self.delta

        # set a central time
        self.pick_time_window([100,200])

        # create a set of waveforms
        tr = obspy.Trace()
        tr.data=self.data[:,0,0,0]
        tr.stats.delta=self.delta

        # low-pass filter
        fmax=self.dsamp/4.
        fmax=np.max(self.flms)
        for kev in range(0,self.Nev):
            for kchan in range(0,self.Nchan):
                for kstat in range(0,self.Nstat):
                    # bandpass the signal
                    tr.data=self.data[:,kev,kchan,kstat]
                    tr.filter('lowpass',freq=fmax,zerophase=True)
                    self.data[:,kev,kchan,kstat]=tr.data


    def savetshfs(self,fnameadd=''):
        """
        write the time shifts to a file
        :param        fnameadd: 
        """

        pass


    def addsynthsignal(self,flm=None,asig=5.):
        """
        add a synthetic signal to the data
        :param       flm: frequency band to filter to
        :param      asig: amplitude of the signal relative to the noise
        """

        if flm is None:
            flm=adat.minmax(self.flms)
        
        # make some random values and filter them
        Ncmp=self.Ncmp
        data=np.random.randn(self.ndata,Ncmp,self.Nchan,self.Nstat)

        # create a set of waveforms
        tr = obspy.Trace()
        tr.data=data[:,0,0,0]
        tr.stats.delta=self.delta

        # start and end time
        tst,tnd=tr.stats.starttime,tr.stats.endtime

        # tapered start and end time
        ttst=tst+self.tmean-self.twin/4.-self.tdata[0]
        ttnd=tst+self.tmean+self.twin/4.-self.tdata[0]

        # overall taper length
        tlen=np.minimum(self.twin,2./np.min(self.flms))

        # go through, filter and taper
        for kcmp in range(0,Ncmp):
            for kchan in range(0,self.Nchan):
                for kstat in range(0,self.Nstat):
                    # bandpass the signal
                    tr.data=data[:,kcmp,kchan,kstat]
                    tr.taper(type='hann',max_length=tlen,
                             max_percentage=0.5,side='both')
                    tr.filter('bandpass',freqmin=flm[0],
                              freqmax=flm[1],zerophase=True)

                    # grab portion of interest
                    tr.trim(starttime=ttst,endtime=ttnd)
                    tr.taper(type='hann',max_length=self.twin*2,
                             max_percentage=0.4,side='both')

                    # normalize to std of 1 in this window
                    tr.data=tr.data/(2 * np.std(tr.data))

                    # back to full length
                    tr.trim(starttime=tst,endtime=tnd,pad=True,
                            fill_value=0.)

                    # replace data
                    data[:,kcmp,kchan,kstat]=tr.data

        # how much to shift the signals at each point
        #self.stshfs = self.tmean+self.twin/5.*np.random.randn(self.Nev,self.Nstat,1,1)
        self.stshfs = self.tmean+15.*np.random.randn(self.Nev,self.Nstat,1,1)
        self.stshfs = self.stshfs - (np.median(self.stshfs)-self.tmean)
        self.stshfs = np.repeat(self.stshfs,self.Nfreq,2)
        self.stshfs = np.repeat(self.stshfs,self.Narr,3)
       


        # decide how much to put on each component
        self.champ=np.random.randn(self.Nev,self.Ncmp)
        self.champ=np.divide(self.champ,
                             np.power(np.sum(np.power(self.champ,2),axis=1,
                                             keepdims=True),0.5))

        # pick amplitudes
        amp=np.std(self.data,axis=0)
        amp=np.median(amp,axis=0)
        amp=amp*asig

        # will just roll the values
        nroll=((self.stshfs-self.tmean)/self.delta).astype(int)
        for kev in range(0,self.Nev):
            for kchan in range(0,self.Nchan):
                for kstat in range(0,self.Nstat):
                    for kcmp in range(0,self.Ncmp):
                        self.data[:,kev,kchan,kstat]=\
                             self.data[:,kev,kchan,kstat]+amp[kchan,kstat]*\
                             self.champ[kev,kcmp]*\
                             np.roll(data[:,kcmp,kchan,kstat],nroll[kev,kstat,0,0])

        # save the synthetic signal
        self.synthsig=data
                    
                    


    def filterdata(self):
        """
        filter to the desired frequency ranges
        """

        # initialize
        self.fdata=np.ndarray(list(self.data.shape)+[1],dtype=float)

        # create a set of waveforms
        tr = obspy.Trace()
        tr.data=self.data[:,0,0,0]
        tr.stats.delta=self.delta

        # taper length
        tlen=np.minimum(self.twin,2./np.min(self.flms))

        # and go through the parameter
        for kev in range(0,self.Nev):
            for kchan in range(0,self.Nchan):
                for kstat in range(0,self.Nstat):
                    # grab the data
                    tr.data=self.data[:,kev,kchan,kstat]
                    
                    # deal with masks and tapering
                    msk=adat.prepfiltmask(tr)
                    tr.taper(type='hann',max_length=tlen,
                             max_percentage=0.5,side='both')

                    for kfreq in range(0,self.Nfreq):
                        # filter 
                        trf=tr.copy()
                        trf.filter('bandpass',freqmin=self.flms[kfreq,0],
                                   freqmax=self.flms[kfreq,1],
                                   zerophase=True)

                        # add back mask?
                        adat.addfiltmask(trf,msk)
                        
                        # grab data
                        self.fdata[:,kev,kchan,kstat,kfreq]=trf.data

        # # mask some columns
        self.ident_okay_columns()



    def ident_okay_events(self,rembad=True):
        """
        identify the events that have all the necessary data
        :param     rembad: remove the events that don'thave them
        """

        if isinstance(self.data,np.ma.masked_array):
            eok=np.sum(self.data.mask,axis=0,keepdims=False)==0
        #else:
            shp=list(self.data.shape)[1:]
            eok=np.ones(shp,dtype=bool)

        # number of records available
        eok=np.sum(np.sum(eok,axis=2),axis=1)     
        eok=eok==self.Nstat*self.Nchan
        self.eok=eok

        # remove if desired
        if rembad:
            self.data=self.data[:,eok,:,:]
            self.Nev=np.sum(eok)
            self.eok=eok[eok]

            if 'nests' in self.__dict__.keys():
                self.nests=self.nests[eok]

	#for item in ds.data.mask: 
		#if item[0][0][1] == True:

    def ident_okay_columns(self):
        """
        identify columns in cdata that will have data
        """
        
        if isinstance(self.data,np.ma.masked_array):
            self.cokm=np.sum(self.data.mask,axis=0,keepdims=True)==0
            self.cokm=self.cokm.reshape(list(self.cokm.shape)+[1])
            self.cokm=np.repeat(self.cokm,self.Nfreq,axis=4)
        else:
            shp=list(self.fdata.shape)
            shp[0]=1
            self.cokm=np.ones(shp,dtype=bool)
        self.cok=self.cokm.flatten()

    def pick_time_window(self,trange=[-50,50]):
        """
        :param      trange: time window to examine and stack, relative to time zero 
                              in self.tdata
                            (will be trange + tshfs when appropriate)
        """

        self.trange=np.atleast_1d(trange).flatten()

        # window duration and center
        self.twin=np.diff(trange)[0]
        self.tmean=np.mean(self.trange)

        # number of points in the time window
        self.nwin = int(np.round(self.twin/self.delta))        

    def initpars(self):
        """
        set some initial values
        """
        
        # just one arrival for now
        self.Narr = 1

        # initialize window center times to tmean
        # not used with dual annealling anyway
        self.tshfs = self.tmean+\
            0.*self.twin*np.random.randn(self.Nev,1,1,1)
        self.tshfs = np.repeat(self.tshfs,self.Nstat,1)
        self.tshfs = np.repeat(self.tshfs,self.Nfreq,2)
        self.tshfs = np.repeat(self.tshfs,self.Narr,3)

        # and indices
        self.ishfs = adat.closest(self.tdata,self.tshfs.flatten()
                                     -self.twin/2.)
        self.ishfs = self.ishfs.reshape(self.tshfs.shape)

        # number of initial components (at least 6, for different 
        # moment tensors
        self.Ncmp = 1

        # to get everything
        self.iev=np.arange(0,self.Nev)
        self.ichan=np.arange(0,self.Nchan)
        self.istat=np.arange(0,self.Nstat)
        self.ifreq=np.arange(0,self.Nfreq)

        # how much to weight a mean value
        self.meanweight=1./self.Nev*0.1 #0.1

    def initadata(self):
        """
        initialize the aligned data from ishfs
        """
        

        # just initialize the matrix # this is where something is going wrong and I don't know why 
        self.adata=self.grabdata(self.ishfs[:,0,0,0],iev=self.iev,
                                 ichan=self.ichan,
                                 istat=self.istat,ifreq=self.ifreq)

        # then replace everything
        self.replaceadata(iev=self.iev,ichan=self.ichan,
                          istat=self.istat,ifreq=self.ifreq)



    def replaceadata(self,iev,ichan,istat,ifreq):

        # need to use different time shifts for each frequency
        for kstat in istat:
            ##print(kstat)
            jstat=np.arange(kstat,kstat+1)
            ##print(jstat)
            for kfreq in ifreq:
                jfreq=np.arange(kfreq,kfreq+1)
                vl=self.grabdata(self.ishfs[:,kstat,kfreq,0],iev=iev,
                                 ichan=ichan,istat=jstat,ifreq=jfreq)
                
                self.adata[:,:,:,kstat,kfreq]=vl[:,:,:,0,0]

    def tshfs2x(self):
        """
        collect the relevant time shifts and put them in an array x
        :return     x: the array
        """

        self.x = self.tshfs.flatten()
     
        

    def initsearch(self,method='DA',tshfmax=None):
        """
        set up the parameters for searching
        :param      method: which method to use for searching
        :param     tshfmax: maximum amount to allow to shift
                             from tmean (default: twin*1.5)
        """

        # max shift
        if tshfmax is None:
            tshfmax=self.twin*.5
        self.tshfmax=tshfmax
        self.bounds=[(self.tmean-self.tshfmax,
                      self.tmean+self.tshfmax)] * \
                  self.Nstat*self.Nev*self.Nfreq*self.Narr
        if method=='Nelder-Mead':
            self.searchfun=optimize.minimize
            self.method='Nelder-Mead'
            self.bounds=None
            self.opttype='local'
        elif method in ['SHGO','shgo']:
            self.searchfun=optimize.shgo
            self.method='SHGO'
            self.opttype='global'
        elif method in ['DA','dual_annealing']:
            self.searchfun=optimize.dual_annealing
            self.method='DA'
            self.opttype='global'
            self.maxiter = 5000
            self.initial_temp = 2325
            self.restart_temp_ratio = 2e-6
        elif method in ['basinhopping','BH']:
            self.searchfun=optimize.basinhopping
    
            self.method='BH'
            self.opttype='global_ub'
        elif method in ['differential_evolution','DE']:
            self.searchfun=optimize.differential_evolution
            self.method='DE'
            self.opttype='global'
        elif method in ['iterative_grid']:
            self.searchfun = optimize.brute
            self.method = 'iterative_grid'
            self.opttype = 'iterative_grid'
           
            self.ranges = [(self.tmean-self.tshfmax,
                      self.tmean + self.tshfmax)]
            self.stepsize = 0.1
            self.Ns = float(self.ranges[0][1] - self.ranges[0][0])/self.stepsize
            self.paramiter_search_range = np.arange(self.ranges[0][0],self.ranges[0][1]+self.stepsize,self.stepsize)
            self.Number_of_iterations = 150
            self.misfit_store = np.zeros(self.Number_of_iterations)
            self.TEMPRANGE = [1000]
            self.ref_misfit = -0.8
        elif method in ['DA']:
            self.opttype = 'DA'
            self.bounds=[(self.tmean-self.tshfmax,
                      self.tmean+self.tshfmax)] * \
                  self.Nstat*self.Nev*self.Nfreq*self.Narr          
       
    def findtshf(self):
        
        # starting point
        self.tshfs2x()

        # optimize
        if self.opttype in ['global']:
            self.res=self.searchfun(func=self.calcmisfit,
                                    bounds=self.bounds,maxiter = self.maxiter,restart_temp_ratio = self.restart_temp_ratio)
        elif self.opttype in ['global_ub']:
            self.res=self.searchfun(func=self.calcmisfit,
                                    x0=self.x)
        elif self.opttype in ['local']:
            self.res=self.searchfun(fun=self.calcmisfit,x0=self.x,
                                   method=self.method)
        elif self.opttype in ['DA']:
             num_dim = self.Nev
             srch = nbr.Searcher(
                       objective=self.calcmisfit,
                       limits=self.bounds,
                       num_samp=10,
                       num_resamp=5,
                       maximize=False,
                       verbose=True
                        )
             srch.update(500)
             srch.plot()
   
 
        self.calcperclow_saved()
    def iterative_grid_fine_tune(self,search_range =5, temp = 1000, Niter = 100):
        #must be run after cross correlation to fine tune the cross correaltion grids
        paramiter_search_range = np.zeros(shape=(self.Nev-1,int(2*search_range)+1))
        print('number of iterations:' + str(Niter))
        for k in range(0,Niter):
            print('iteration:'+str(k))
            for i in range(1,self.Nev):
                ev_search_range = np.arange(self.ishfs.flatten()[i] - search_range,self.ishfs.flatten()[i] + search_range +1,1)
                misfit = np.zeros(len(ev_search_range))
                for j in range(0,len(ev_search_range)):
                        #print(j)
                        x = self.ishfs.flatten()
                        x[i] = ev_search_range[j]
                        misfit[j]=self.calcmisfit_i(x)
                point = np.where(misfit == min(misfit))
                self.ishfs[i] = ev_search_range[point]	
	
        return min(misfit)

    def plot_svd_frc(self):
        plt.close()
        
        f = plt.figure()
        p = plt.axes()

        # cumulative numbers in SVD
        scum=np.cumsum(self.S)
        scum=np.append([0],scum/scum[-1])

        # plot sum, fractional
        nvl=scum.size-1
        p.plot(np.arange(0,nvl+1),scum,zorder=2,color='b',linewidth=2)
        
        p.plot([0,nvl],[0,1],color='k',linestyle='--',zorder=0,linewidth=1)

        p.plot([self.Ncmp]*2,[-1,2],color='brown',linewidth=1,zorder=1,
               linestyle=':')

        fs='medium'
        p.set_xlabel('principal component number',fontsize=fs)
        p.set_ylabel('cumulative energy accommodated',fontsize=fs)
        p.set_ylim([-0.02,1.02])
        p.set_xlim([0,nvl])
        plt.show()

    def xcorr(self,master):
        ishift = np.zeros(shape = (self.Nev -1 ))
        cc_store = np.zeros(shape = self.Nev)
        cc_min = np.zeros(shape = self.Nev)
        for i in range(0,self.Nev):
           cc = (correlate_template(self.adata[:,i,0,0,0],master))
           c2 = cc
           cc= abs(cc)
           cc_store[i] = max(cc)
           cc_min[i] = min(c2)
           ishift[i-1] = np.where(np.array(cc) == max(cc))[0][0]
           ishift[i-1] = self.ishfs[0,0,0,0] + int(ishift[i-1] - (len(self.adata[:,0,0,0,0])/2))
        self.ishfs[1:,0,0,0] = ishift
        return cc_store,cc_min

    def plotU(self,plotsynth=None,plotdata=False,Ncmp=None,ichan=None,istat=None,ifreq=None):

        if Ncmp is None:
            Ncmp=self.Ncmp

        if ichan is None:
            ichan=np.arange(0,self.Nchan)
        ichan=np.atleast_1d(ichan)
        Nchan=len(ichan)

        if istat is None:
            istat=np.arange(0,self.Nstat)
        istat=np.atleast_1d(istat)
        Nstat=len(istat)

        if ifreq is None:
            ifreq=np.arange(0,self.Nfreq)
        ifreq=np.atleast_1d(ifreq)
        Nfreq=len(ifreq)

        jchan,jstat,jfreq=np.meshgrid(ichan,istat,ifreq)
        jchan,jstat,jfreq=jchan.flatten(),jstat.flatten(),jfreq.flatten()
        N2=len(jchan)

        plt.close()
        f = plt.figure(figsize=(N2*5+.5,Ncmp*2))
        gs,p=gridspec.GridSpec(Ncmp,N2),[]
        gs.update(left=0.1,right=0.97)
        gs.update(bottom=0.07,top=0.97)
        gs.update(hspace=0.1,wspace=0.1)
        for gsi in gs:
            p.append(plt.subplot(gsi))
        p=np.array(p).reshape([Ncmp,N2])

        # whether to flip values
        tvl=np.arange(0,self.nwin)*self.delta
        tvl=tvl-np.mean(tvl)+self.tmean
        if 'synthsig' in self.__dict__.keys():
            shfmed=np.median(self.tshfs-self.stshfs,axis=0)
        else:
            shfmed=np.zeros([self.Nstat,self.Nchan,self.Nfreq],dtype=float)
        xcs=np.ones([Ncmp,N2],dtype=float)
        xcsd=np.ones([Ncmp,N2],dtype=float)

        if plotsynth is None:
            plotsynth='synthsig' in self.__dict__.keys()

        if plotsynth:
            for k in range(0,Ncmp):
                for kc in range(0,N2):
                    yvls=self.synthsig[:,k,jchan[kc],jstat[kc]]
                    p[k,kc].plot(self.tdata,yvls/np.max(np.abs(yvls)),color='k')

                    # cross-correlate to see if the PCs are flipped relative to U
                    xch=np.correlate(yvls,self.U[:,k,jchan[kc],jstat[kc],jfreq[kc]],
                                     mode='full')
                    if np.abs(np.min(xch))>np.max(xch):
                        xcs[k,kc]=-1.

                    # cross-correlate to see if the PCs are flipped relative to the data
                    xch=np.correlate(yvls,self.fdata[:,0,jchan[kc],jstat[kc],jfreq[kc]],
                                     mode='full')
                    if np.abs(np.min(xch))>np.max(xch):
                        xcsd[k,kc]=-1.


        if plotdata:
            for k in range(0,Ncmp):
                for kc in range(0,N2):
                    tshf=self.stshfs[0,jstat[kc],jfreq[kc],0]-self.tmean
                    yvl=self.fdata[:,0,jchan[kc],jstat[kc],jfreq[kc]]
                    p[k,kc].plot(self.tdata-tshf,
                                 xcsd[k,kc]*yvl/np.max(np.abs(yvl)),color='b')

                        

        i1=adat.closest(self.tdata,tvl)[0]

        for k in range(0,Ncmp):
            for kc in range(0,N2):
                yvl=self.U[:,k,jchan[kc],jstat[kc],jfreq[kc]]
                tmed=shfmed[jstat[kc],jfreq[kc],0]
                p[k,kc].plot(tvl+tmed,xcs[k,kc]*yvl/np.max(np.abs(yvl)),color='r')
                p[k,kc].set_xlim(adat.minmax(tvl))


    def calcmisfit(self,x):
        """
        calculate the misfit given a set of time shifts x
        :param       x: the concatenated time shifts
        :return   msft: a misfit accommodating the fraction in the 
                  first few components and the allowable time shifts
        """

        # set the time shifts
        self.x2tshfs(x)

        # grab shifted data
        self.replaceadata(iev=self.iev,ichan=self.ichan,
                          istat=self.istat,ifreq=self.ifreq)
  
        # calculate the fraction of the variance in the first components
        frc=self.calcperclow()
        #print(frc)
     

        # some regularization to make the shifts not all move
        #msft=self.meanweight*np.abs(np.mean(x)-self.tmean)/(self.twin/3)

        # for now, don't worry about the relative time shifts
        self.msft=-frc

        return self.msft


    def calcmisfit_i(self,x):
        self.isfhts_reshape(x)
        self.replaceadata(iev=self.iev,ichan=self.ichan,
                          istat=self.istat,ifreq=self.ifreq)
  
        # calculate the fraction of the variance in the first components
        frc=self.calcperclow()

        # for now, don't worry about the relative time shifts
        self.msft=-frc
        return self.msft
    
    
    def calcmisfit_loop(self,event_shift):
        """
        calculate the misfit given a set of time shifts x
        :param       x: the concatenated time shifts
                     event_shift: the shift for one event that is being looped over 
		     i: the index of the event that has been looped over
        :return   msft: a misfit accommodating the fraction in the 
                  first few components and the allowable time shifts
        """

        #replace the elements for the event 
        self.x[self.event_iterable] = event_shift 

        # set the time shifts
        self.x2tshfs(self.x)

        # grab shifted data
        self.replaceadata(iev=self.iev,ichan=self.ichan,
                          istat=self.istat,ifreq=self.ifreq)

        # calculate the fraction of the variance in the first components
        self.frc=self.calcperclow()
        ##print(frc)
     

        # some regularization to make the shifts not all move
        #msft=self.meanweight*np.abs(np.mean(self.x)-self.tmean)/(self.twin/3)

        # for now, don't worry about the relative time shifts
        self.msft=-self.frc
               
        return self.msft


    def reshape_final_misfit_vector(self):
       
       
       self.reshape_misfit = np.array(self.final_misfit_vector_for_plotting).reshape((int(self.Ns)+1,int(self.Nev)))

       return self.reshape_misfit.T


    def x2tshfs(self,x):
        """
        take an array of time shifts, and organize them  into 
        an [Nev x Nstat x Nfreq x Narr] matrix self.tshfs
        also set self.ishfs, the indices of the shifts to use
        :param      x: an array of time shifts, concatenated
        """
        
        # closest indices
        self.ishfs = adat.closest(self.tdata,x-self.twin/2.)
        ##print(self.ishfs)

        # desired shape
        self.ishfs = self.ishfs.reshape([self.Nev,self.Nstat,self.Nfreq,self.Narr])
        self.tshfs = x.reshape([self.Nev,self.Nstat,self.Nfreq,self.Narr])

    def isfhts_reshape(self,x):
        x = x.astype(int)
        self.ishfs = x.reshape([self.Nev,self.Nstat,self.Nfreq,self.Narr])

        ##print(self.ishfs.flatten())

    def grabdata(self,ishfs=None,iev=None,ichan=None,istat=None,ifreq=None):
        """
        :param        ishfs: shifts for each event
                                 starting indices of the windows
              ishfs is assumed to be station and channel independent
        :param          iev: which events are being extracted
        :param        ichan: indices of the channels to extract int 
        :param        istat: indices of the stations to extract
        :param        ifreq: indices of the frequency bands to extract
        :return       edata: extracted and aligned data
        """


        # grab the data
        shp1=list(self.fdata.shape)
        #ff#print(ichan.size,istat.size,ifreq.size)
        shp1=shp1[0:2]+[ichan.size,istat.size,ifreq.size]
        shp=[self.nwin,1,ichan.size,istat.size,ifreq.size]
        
        edata=self.fdata[:,:,ichan,:,:][:,:,:,istat,ifreq].reshape(shp1)

        
        edata=np.hstack([edata[ishfs[k]:ishfs[k]+self.nwin,
                               iev[k]:iev[k]+1,:,:,:].reshape(shp)
                         for k in range(0,iev.size)])
        for i in range(0,self.Nev):
             edata[:,i,:,:] = (edata[:,i,:,:]- np.mean(edata[:,i,:,:]))/np.std(edata[:,i,:,:])
             #edata[:,i,:,:] = edata[:,i,:,:]/np.max(abs(edata[:,i,:,:]))         


        return edata


    def datato2d(self,edata):
        """
        :param        edata: some extracted data set, with events 
                              indexed in the second dimension
        :return       cdata: data with third to last dimensions
                              concatenated along the column
        """

        # shape
        shp=list(edata.shape)
        if len(shp)==2:
            cdata=edata.copy()
        else:
            # change dimension ordering
            itrans=np.arange(0,len(shp))
            itrans=np.append(np.delete(itrans,1),1)
            cdata=edata.transpose(itrans)

            # and flatten first few dimensions
            cdata=cdata.reshape([shp[0]*np.prod(shp[2:]),shp[1]])

            # and flatten first few dimensions
        return cdata

    def datafrom2d(self,cdata):
        """
        :param       cdata: data with third to last dimensions
                              concatenated along the column
        :return      edata: extracted data
        """

        # the shape
        shp=[1,cdata.shape[1],self.Nchan,self.Nstat,self.Nfreq]
        shp[0]=int(np.round(cdata.shape[0]/np.prod(shp[2:])))
        shp2=shp[0:1]+shp[2:]+shp[1:2]

        # reshape
        edata=cdata.reshape(shp2)

        # and re-order dimensions
        itrans=np.arange(0,len(shp)-1)
        itrans=np.insert(itrans,1,len(shp)-1)             
        edata=edata.transpose(itrans)


        return edata
	

    def calcsvd(self,cdata):
        """
        calculate the svd of the collected data
        :param          cdata: the collected data
        :return             U: the potential time dependencies
        :return             S: the strengths of each one
        :return             V: the component of each time dependent on 
                                  each event
        """

        U,S,V=np.linalg.svd(cdata,full_matrices=False)
        
        return U,S,V


    def calcperclow(self):
        """
        calculate the percentage of the variance accommodated by the
        by the first few components of the SVD
        :return     frc: fraction in first self.Ncmp components
        """

        # re-organize input data
        cdata =self.datato2d(self.adata)

        # SVD
        U,S,V=self.calcsvd(cdata)

        # cumumulative S
        frc = np.sum(S[0:self.Ncmp])/np.sum(S)
        ##print(S[0],S[1],S[2],S[3],S[4],S[-1])
	

        ##print(self.Ncmp,S[0:self.Ncmp],np.sum(S))

        return frc

    def calcperclow_saved(self):
        """
        calculate the percentage of the variance accommodated by the
        by the first few components of the SVD
        :return     frc: fraction in first self.Ncmp components
        """

        # grab shifted data
        self.replaceadata(iev=self.iev,ichan=self.ichan,
                          istat=self.istat,ifreq=self.ifreq)

        # re-organize input data
        cdata =self.datato2d(self.adata)

        # SVD
        self.U,self.S,self.V=self.calcsvd(cdata)
        self.U=self.datafrom2d(self.U)

        # cumumulative S
        self.frc = np.sum(self.S[0:self.Ncmp])/np.sum(self.S)
        
    def calcperclow_auto(self):
        """
        calculate the percentage of the variance accommodated by the
        by the first few components of the SVD
        :return     frc: fraction in first self.Ncmp components
        """


        # re-organize input data
        cdata =self.datato2d(self.adata)

        # SVD
        self.U,self.S,self.V=self.calcsvd(cdata)
        self.U=self.datafrom2d(self.U)

        # cumumulative S
        self.frc = np.sum(self.S[0:self.Ncmp])/np.sum(self.S)

    def calcperclow_saved_cmp(self,flatten_tshfts):
        """
        calculate the percentage of the variance accommodated by the
        by the first few components of the SVD
        :return     frc: fraction in first self.Ncmp components
        """

        # grab shifted data
        self.replaceadata(iev=self.iev,ichan=self.ichan,
                          istat=self.istat,ifreq=self.ifreq)
 

        # re-organize input data
        cdata =self.datato2d(self.adata)

        # SVD
        self.U,self.S,self.V=self.calcsvd(cdata)
        self.U=self.datafrom2d(self.U)

        # cumumulative S
        self.frc = np.sum(self.S[0:self.Ncmp])/np.sum(self.S)
	



