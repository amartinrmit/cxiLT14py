"""
cxilt14py
- Analysis scripts for performing angular correlation analysis on the cspad detector @ CXI @ LCLS

correlation.py - calculate angular correlation function

September 2018

authors: Andrew Martin (andrew.martin@rmit.edu.au)
         

"""

import numpy as np
import scipy.ndimage as sdn


class angular_correlation:

    
    def polar_plot( x, y, data, nr, nth, rmin, rmax, thmin, thmax ):

        # r and theta arrays
        rarr = np.outer( np.arange(nr)*(rmax-rmin)/float(nr) + rmin, np.ones(nth) )
        tharr = np.outer( np.ones(nr), np.arange(nth)*(thmax-thmin)/float(nth) + thmin)
        
        newx = rarr*np.cos( tharr )
        newy = rarr*np.sin( tharr )
        
        newdata = sdn.map_coordinates( data, [newx.flatten(), newy.flatten()], order=3 )

        return newdata.reshape( nr, nth )

    def polarplot_angular_correlation( polar, polar2=None):

        fpolar = np.fft.fft( polar, axis=2 )

        if polar2 != None:
            fpolar2 = np.fft.fft( polar2, axis=2)
        else:
            fpolar2 = fpolar

        out = np.fft.ifft( fpolar2.conjugate() * fpolar )
        
        return out

        
    def apply_mask( func, mask ):
        return func*mask


    def correct_mask_correlation( corr, maskcorr ):
        imask = np.where( maskcorr != 0 )
        corr[imask] *= 1.0/maskcorr[imask]
        return corr

        
