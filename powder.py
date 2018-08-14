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

    def powderplot( self, qarr, data, nqbins, qmin=-1, qmax=-1): 
        
        # output arrays 
        pplot = np.zeros( nqbins )

        # step size
        qstep = (qmax-qmin)/float(nqbins)
        
        #
        # loop over qbins
        #
        for i in np.arange(qbins):
            dq = qmin + i*qstep
            iq = np.logical_and( (qarr>dq), qarr<dq+qstep )
            pplot[i] = data.mean( iq )

        return pplot
