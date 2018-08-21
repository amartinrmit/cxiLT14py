"""
cxilt14py
- Analysis scripts for performing angular correlation analysis on the cspad detector @ CXI @ LCLS

powder.py - calculate powder lineplots

September 2018

authors: Andrew Martin (andrew.martin@rmit.edu.au)
         

"""

import numpy as np


class powder:

#    def __init__(self, q2):

#        self.qgrids

    #
    # Precalculate the q indices for 1D powder plot
    #
    def set_up_q_indices( self, qarr, nqbins=500, qmin=-1, qmax=-1  ):

        # step size
        self.qstep = (qmax-qmin)/float(nqbins)
        self.nqbins = nqbins
        self.qsamps = np.zeros( nqbins )

        #
        # loop over qbins
        #
        self.iq = []
        for i in np.arange(nqbins):
            dq = qmin + i*self.qstep
            self.qsamps[i] = dq
#            iq = np.logical_and( (qarr>dq), qarr<dq+qstep )
            self.iq.append( np.where( (qarr>dq) * (qarr<dq+self.qstep) ) )


    #
    # otf - " one the fly".
    #
    # Uses precalcualted q indices for each q shell
    #
    def powderplot_otf( self, data, mask, stdout=False, ccout=False): 
        
        nqbins = self.nqbins
        # output arrays 
        pplot = np.zeros( nqbins )

        if stdout==True:
            splot = np.zeros( nqbins )

        if ccout==True:
            ccplot = np.zeros( nqbins )

        # step size
        qstep = self.qstep 
        
        # apply mask
        dtmp = data*mask

        #
        # loop over qbins
        #
        for i in np.arange(nqbins):
            iq = self.iq[i]

            pplot[i] = dtmp[ iq].sum()
            norm = mask[iq].sum()
            if norm>0:
                pplot[i] *= 1./norm

            if stdout==True:
                splot[i] = np.sum(dtmp[ iq]**2)
                if norm>0:
                    splot[i] *= 1./norm

            if ccout==True:
                shift =  iq[0].size/3
                iqr = ( np.roll( iq[0], shift ), np.roll( iq[1], shift ), np.roll( iq[2], shift ) )
                #                print len(iqr), len(iq), iq[0].shape

                ccnorm = np.sum(mask[ iq] * mask[iqr] ) 
#                ccplot[i] = np.sum((dtmp[ iq]-pplot[i]) * (dtmp[iqr]-pplot[i]) ) 
                ccplot[i] = np.sum(dtmp[ iq] * dtmp[iqr] ) 

                if ccnorm>0:
                    ccplot[i] *= 1./ccnorm
                    ccplot[i] += - pplot[i]*pplot[i]

        if stdout==True:
            splot = np.sqrt( splot - pplot*pplot) 

        output = [pplot]

        if stdout==True:
            output.append( splot )
        if ccout==True:
            output.append( ccplot )

        return output




    def powderplot( self, data, qarr, mask, nqbins, qmin=-1, qmax=-1,\
                    stdout=False, ccout=False): 
        
        # output arrays 
        pplot = np.zeros( nqbins )
        qsamps = np.zeros( nqbins )
        if stdout==True:
            splot = np.zeros( nqbins )

        if ccout==True:
            ccplot = np.zeros( nqbins )

        # step size
        qstep = (qmax-qmin)/float(nqbins)
        
        # apply mask
        dtmp = data*mask

        #
        # loop over qbins
        #
        for i in np.arange(nqbins):
            dq = qmin + i*qstep
            qsamps[i] = dq
#            iq = np.logical_and( (qarr>dq), qarr<dq+qstep )
            iq = np.where( (qarr>dq) * (qarr<dq+qstep) )
            pplot[i] = dtmp[ iq].sum()
            norm = mask[iq].sum()
            if norm>0:
                pplot[i] *= 1./norm

            if stdout==True:
                splot[i] = np.sum(dtmp[ iq]**2)
                if norm>0:
                    splot[i] *= 1./norm

            if ccout==True:
                shift =  iq[0].size/4
                iqr = ( np.roll( iq[0], shift ), np.roll( iq[1], shift ), np.roll( iq[2], shift ) )
#                print len(iqr), len(iq), iq[0].shape
                ccnorm = np.sum(mask[ iq] * mask[iqr] ) 
                ccplot[i] = np.sum(dtmp[ iq] * dtmp[iqr] ) 
                if ccnorm>0:
                    ccplot[i] *= 1./ccnorm
                ccplot[i] += - pplot[i]*pplot[i]


        if stdout==True:
            splot = np.sqrt( splot - pplot*pplot) 

        output = [pplot, qsamps]

        if stdout==True:
            output.append( splot )
        if ccout==True:
            output.append( ccplot )

        return output

    #
    # histogram of pixel values within a certain q range and angular range
    #
    def qhistogram( self, data, qarr, qmin, qmax, nbins=100,\
                    dmin=-1e10, dmax=-1e10):

        # get indices of pixels within the q range
        iq = np.logical_and( (qarr>qmin), qarr<qmax )
        if len(iq)==0:
            print "powder.py; qhistogram(); no pixel found within qrange."
            return [], []

        if dmin == -1e10:
            dmin = np.min( data[iq])
        if dmax < -1e10:
            dmax = np.max( data[iq])

        h, b = np.histogram( data[iq], range=(dmin,dmax), nbins=nbins )

        return h, b

    #
    # Precalculate the q indices for 1D powder plot
    #
    def set_up_qth_indices( self, qarr, tharr, nqbins=500, qmin=-1, qmax=-1,\
                          nthbins=500, thmin=0, thmax=2*3.14159):

        # step size
        self.qstep = (qmax-qmin)/float(nqbins)
        self.nqbins = nqbins
        self.qthsamps = np.zeros( (nqbins,nthbins,2) )

        self.thstep = (thmax - thmin) / float(nthbins)
        self.nthbins = nthbins 

        #
        # loop over qbins
        #
        self.iqth = []
        for i in np.arange(nqbins):
            print i
            dq = qmin + i*self.qstep
            for j in np.arange(nthbins):
                dth = thmin + j*self.thstep
                self.qthsamps[i,j,0] = dq
                self.qthsamps[i,j,1] = dth
                #            iq = np.logical_and( (qarr>dq), qarr<dq+qstep )
                self.iqth.append( np.where( (qarr>dq) * (qarr<dq+self.qstep) * (tharr>dth)*(tharr<dth+self.thstep ) ))

    #
    # otf - " one the fly".
    #
    # Uses precalcualted q indices for each q shell
    #
    def polarplot_otf( self, data, mask, stdout=False, ccout=False): 
        
        nqbins = self.nqbins
        nthbins = self.nthbins
        # output arrays 
        pplot = np.zeros( (nqbins, nthbins) )

        if stdout==True:
            splot = np.zeros( (nqbins, nthbins) )

        if ccout==True:
            ccplot = np.zeros( (nqbins, nthbins) )

        # step size
        qstep = self.qstep 
        thstep = self.thstep
        
        # apply mask
        dtmp = data*mask

        #
        # loop over qbins
        #
        for i in np.arange(nqbins):
            for j in np.arange(nthbins):
                iqth = self.iqth[i,j]

                pplot[i,j] = dtmp[ iqth].sum()
                norm = mask[iqth].sum()
                if norm>0:
                    pplot[i,j] *= 1./norm

            if stdout==True:
                splot[i,j] = np.sum(dtmp[ iqth]**2)
                if norm>0:
                    splot[i,j] *= 1./norm

#            if ccout==True:
#                shift =  iq[0].size/3
#                iqr = ( np.roll( iq[0], shift ), np.roll( iq[1], shift ), np.roll( iq[2], shift ) )
 

#                ccnorm = np.sum(mask[ iq] * mask[iqr] ) 
#                ccplot[i] = np.sum(dtmp[ iq] * dtmp[iqr] ) 

 #               if ccnorm>0:
#                    ccplot[i] *= 1./ccnorm
#                    ccplot[i] += - pplot[i]*pplot[i]

        if stdout==True:
            splot = np.sqrt( splot - pplot*pplot) 

        output = [pplot]

        if stdout==True:
            output.append( splot )
#        if ccout==True:
#            output.append( ccplot )

        return output
