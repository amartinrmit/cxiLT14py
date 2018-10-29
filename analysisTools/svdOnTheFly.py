
import numpy as np
import matplotlib.pyplot as plt
import h5py


class SVDthin:

    def __init__(self, rnkmax=10):

        self.ndata = 0
        self.slist = []
        self.rank_max = rnkmax
        self.rank_current = 0
        self.ulist = []

    def calc_coeffs( self, data, smax=-1, nstart=0 ):
        if smax==-1:
            smax = len(self.ulist)

        c = []
        for i in np.arange(smax)+nstart:
            c.append( np.sum(self.ulist[i]*data) )
        return c

    def project_coefficients( self, c, smax=-1, nstart=0 ):
        if smax == -1:
            smax = len(self.ulist)


        output = c[0]*self.ulist[nstart]
        for i in np.arange(smax-1-nstart)+1+nstart:
            output += c[i-nstart]*self.ulist[i]
        return output

    def project_data( self, data, smax=-1, nstart=0):
        c = self.calc_coeffs( data, smax, nstart)
        output = self.project_coefficients( c, smax, nstart )
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

    #
    # reads data from a h5 file
    #
    def h5read_svdmodes( self, filename, shape ):
        h5file = h5py.File(filename,"r")
        # print field
        
        # get rank
        field = "/parameters/"+"rank"
        self.current_rank = h5file[field][...]        

        # get eigenvalues
        #field = "/singular_values"
        #sv = h5file[field][...]        
        #self.update_eigenvalues( s )
        
        # get svd modes
        self.ulist = []
        for i in np.arange( self.current_rank):
            field = "/svdmodes/mode"+str(i)
            self.ulist.append( h5file[field][...].reshape(shape)  )

        h5file.close()

    #
    # writes an array into a h5 file
    #
    def h5write_svdmodes(self, filename):
        f = h5py.File(filename, 'w')    # overwrite any existing file

        # store rank
        field = "/parameters/"+"rank"
        f.create_dataset(field, data=self.rank_current)

        # store eigenvalues
        field = "/singular_values"
        f.create_dataset(field, data=np.array(self.slist) )

        # store svd modes
        for i, data in enumerate(self.ulist):
            field = "/svdmodes/mode"+str(i)
            dset = f.create_dataset(field, data=data)

        f.close()
