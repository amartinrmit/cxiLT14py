
import numpy as np
import matplotlib.pyplot as plt


class SVDthin:

    def __init__(self, rnkmax=10):

        self.ndata = 0
        self.slist = []
        self.rank_max = rnkmax
        self.rank_current = 0
        self.ulist = []

    def calc_coeffs( self, data ):
        
        c = []
        for i in np.arange(len(self.ulist) ):
            c.append( np.sum(self.ulist[i]*data) )
        return c

    def project_coefficients( self, c ):

        output = c[0]*self.ulist[0]
        for i in np.arange(len(self.ulist)-1)+1:
            output += c[i]*self.ulist[i]
        return output

    def update_eigenvalues( self, s ):

        for i in np.arange(self.rank_current):
            self.slist[i] = s[i]
            
        if self.rank_current < self.rank_max:
            self.slist.append( s[self.rank_current] )

    def update_eigenvectors( self, urot, p  ):

        uout = np.zeros( (self.rank_current+1, self.ulist[0].size) )
        u = np.array( self.ulist )
        uout[:self.rank_current,:] = u
        uout[self.rank_current,:] = p
        
#        print "u shapes", uout.shape, urot.shape
        uout = np.dot( np.transpose(urot), uout )
#        uout = np.dot( urot, uout )

 #       print "uout shape", uout.shape

        for i in np.arange(self.rank_current):
            self.ulist[i] = uout[i]
            
        if self.rank_current < self.rank_max:
            self.ulist.append( uout[self.rank_current] )
            
    def update_rank( self ):
        if self.rank_current < self.rank_max:
            self.rank_current += 1

    def add_column( self, data ):

        #
        # initialize variables if current rank is zero
        #
        if self.rank_current == 0:
            self.ulist.append( data / np.sqrt(np.sum(data**2)) )
            self.slist.append( np.sqrt(np.sum(data**2)) )
            self.update_rank()
        else:
            #
            # calculate projection of data onto u matrices
            # and project back again... (m vector)
            #
            c = self.calc_coeffs( data )
            proj = self.project_coefficients( c )

            #
            # construct p vector
            #
            p = data - proj
            pmod = np.sqrt( np.sum(p*p))

            #
            # make K matrix
            #
            kmat = np.zeros( (self.rank_current+1, self.rank_current+1))
            sdiag = np.diag(np.array(self.slist))
            kmat[:self.rank_current, :self.rank_current] = sdiag
            kmat[:self.rank_current, self.rank_current] = np.array(c)
            kmat[self.rank_current, self.rank_current] = pmod

            #print "kmat", kmat

            #
            # diagonalize K matrix
            #
            u, s, v = np.linalg.svd( kmat )
            
            #print "kmat u", u
            #print "kmat v", v
            #
            # update the eigenvalues
            #
            self.update_eigenvalues( s )

            #
            # Use u' and v' matrices to update the eigenvectors
            #
            self.update_eigenvectors( u, p/pmod )
            self.update_rank()
            
            self.ndata += 1
