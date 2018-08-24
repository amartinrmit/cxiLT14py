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

    
    def polar_plot( self,data, nr, nth, rmin, rmax, thmin, thmax, cenx, ceny, submean=False ):

        # r and theta arrays
        rarr = np.outer( np.arange(nr)*(rmax-rmin)/float(nr) + rmin, np.ones(nth) )
        tharr = np.outer( np.ones(nr), np.arange(nth)*(thmax-thmin)/float(nth) + thmin)
        
        newx = rarr*np.cos( tharr ) + cenx
        newy = rarr*np.sin( tharr ) + ceny
        
        newdata = sdn.map_coordinates( data, [newx.flatten(), newy.flatten()], order=3 )

        out = newdata.reshape( nr, nth )
        if submean == True:
            out = self.polar_plot_subtract_rmean( out  )

        return out

    def polar_plot_subtract_rmean( self, pplot ):

        av = np.average( pplot, 1 )
        out = pplot -np.outer( av, np.ones( pplot.shape[1] ) )
        return out

    def polarplot_angular_correlation( self, polar, polar2=None):

        fpolar = np.fft.fft( polar, axis=1 )

        if polar2 != None:
            fpolar2 = np.fft.fft( polar2, axis=1)
            out = np.fft.ifft( fpolar2.conjugate() * fpolar, axis=1 )
        else:
            out = np.fft.ifft( fpolar.conjugate() * fpolar, axis=1 )
       
        return out

        
    def apply_mask( self, func, mask ):
        return func*mask


    def mask_correction( self, corr, maskcorr ):
        imask = np.where( maskcorr != 0 )
        corr[imask] *= 1.0/maskcorr[imask]
        return corr

    #    
    # pairwise correlation of (flattened) arrays
    #
    # not for angular correlations; good for correlation of mean asic values
    #
    def allpixel_correlation( self, arr1, arr2 ):
        out = np.outer( arr1.flatten(), arr2.flatten() )
        return out

    # pearson correlation of a 2D area
    def pearsonCorrelation2D( self, arr1, arr2, lim=None):
        if lim == None:
            lim = [0, arr1.shape[0], 0, arr1.shape[1]]
      
        a1 = arr1[lim[0]:lim[1],lim[2]:lim[3]]
        a2 = arr2[lim[0]:lim[1],lim[2]:lim[3]]
        
        c1 = a1 - np.average(a1)
        c2 = a2 - np.average(a2)
        pc = np.sum( c1*c2 ) /np.sqrt( np.sum(c1*c1) * np.sum(c2*c2))
        return pc

# returns pearson correlation of each q ring
    def pearsonCorrelation2D_angular( self, arr1, arr2, lim=None):
        if lim == None:
            lim = [0, arr1.shape[0], 0, arr1.shape[1]]
      
        a1 = arr1[lim[0]:lim[1],lim[2]:lim[3]]
        a2 = arr2[lim[0]:lim[1],lim[2]:lim[3]]
        
        c1 = a1 - np.outer( np.average(a1, 1), np.ones( a1.shape[1]) )
        c2 = a2 - np.outer( np.average(a2, 1), np.ones( a2.shape[1]) )
        pc = np.sum( c1*c2, 1 ) /np.sqrt( np.sum(c1*c1, 1) * np.sum(c2*c2, 1))
        return pc

#    def gaussian_filter_correlation( self, corr, qsig, thsig ):
        
        # linspace for q and th values... 
        # then turn that into guassian function
        # then call convoluation function

    def process_asiccorrsum( self, corr, datasum, maskasic, n):

        # divide by number of processed frames
        c = corr / float(n)
        d = datasum / float(n)

        # correct asic correlation with mean
        mean_asic_corr = self.allpixel_correlation( d, d)
        c += - mean_asic_corr
        
        # pearson correlation of offdiagonal parts
        offdiag0 = c-np.diag(np.diag(c))

        return c, mean_asic_corr, offdiag0
