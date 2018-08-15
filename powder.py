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

    def powderplot( self, data, qarr, nqbins, qmin=-1, qmax=-1,\
                    stdout=False): 
        
        # output arrays 
        pplot = np.zeros( nqbins )
        if stdout==True:
            splot = np.zeros( nqbins )

        # step size
        qstep = (qmax-qmin)/float(nqbins)
        
        #
        # loop over qbins
        #
        for i in np.arange(qbins):
            dq = qmin + i*qstep
            iq = np.logical_and( (qarr>dq), qarr<dq+qstep )
            pplot[i] = data.mean( iq )
            if stdout==True:
                splot[i] = data.std( iq )

        if stdout==True:
            return pplot, splot
        else:
            return pplot

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
