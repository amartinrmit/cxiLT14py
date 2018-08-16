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
    # histogram of pixel values within a certain q range
    #
    def qhistogram( self, data, qarr, qmin, qmax, nbins=100,
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
