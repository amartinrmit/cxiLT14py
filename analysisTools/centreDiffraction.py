

import sys
import numpy as  np


class check_centre:

     def __init__(self, nx=0, ny=0, rmin=0, rmax=100):

          self.nx, self.ny = nx, ny
          self.rmin = rmin
          self.rmax = rmax
          x = np.outer( np.arange( nx) - nx/2, np.ones( ny ) )
          y = np.outer( np.ones( nx ), np.arange( ny ) - ny/2 )
          r = np.sqrt( x*x + y*y )
          self.ir = np.where( (r > rmin)*(r<rmax) )

     def centrosymmetry_metric( self, data, mask ):

          rotdata = np.rot90( data, 2 )
          rotmask = np.rot90( mask, 2 )
          tmask = mask*rotmask
          error = np.sum( np.abs( (rotdata[self.ir] - data[self.ir])*tmask[self.ir])  ) / np.sum( tmask[self.ir] )

          return error

     def line_search( self, data, mask, circ, dim='x', dim1mid=0, dim2shift=0, npoints=10, step=1 ):

          dim1start = dim1mid - (step*npoints)/2
          
          if dim=='x':
               dim1, dim2 = 0, 1
          elif dim=='y':
               dim1, dim2 = 1, 0
          else:
               print "line search error incorrect dimensions"
               sys.exit()

          data_dim2shift = np.roll( data, int(dim2shift), dim2 )
          mask_dim2shift = np.roll( mask, int(dim2shift), dim2 )
          minerr = 1e10
          err = np.zeros( npoints )
          for i in np.arange( npoints):
     
               tmp = np.roll( data_dim2shift, int(i+dim1start), dim1)
               tmp_mask = np.roll( mask_dim2shift, int(i+dim1start), dim1) * circ

               err[i] = self.centrosymmetry_metric( tmp, tmp_mask )
               if (err[i] < minerr):
                    minerr = err[i]

#          print "debug centreDiffraction.py linesearch", err[i], np.min(err), np.max(err), np.max(tmp), np.max(circ)
          pos = np.where( err == np.min(err) )[0][0] + dim1start
#          print err                                                                                                                          
          print pos, np.where( err == minerr )[0][0]                                                      
#          plt.plot( err )                                                                                                                   
#          plt.draw()                                                                                      
#          plt.show()                                                                                                                          
          return pos, minerr



     def twoD_line_search( self, data, mask, circ_filter, nsearch=1, xlim = [-10,10], ylim=[-10,10], xsteps=80, ysteps=80,\
                           verbosity=1):

          dx = (xlim[1]-xlim[0])/float(xsteps)
          dy = (ylim[1]-ylim[0])/float(ysteps)
          xpos = (xlim[1]+xlim[0])/2.0
          ypos = (ylim[1]+ylim[0])/2.0
         
          if verbosity>0:
               print "2D line search for pattern centre."
               print "line search parameter xrange, yrange, xsteps, ysteps, niterations", xlim, ylim, xsteps, ysteps, nsearch
          for i in np.arange(nsearch):
               xpos, xerr =  self.line_search( data, mask, circ_filter, dim='x', dim1mid=xpos, dim2shift=ypos, npoints=xsteps )
               ypos, yerr =  self.line_search( data, mask, circ_filter, dim='y', dim1mid=ypos, dim2shift=xpos, npoints=ysteps )
               if verbosity>0:
                    print "iteration ", i, "current centre", xpos, ypos, " x/y line search errors ", xerr, yerr

          if verbosity>0:
               print "line search final result. Pattern centre", xpos, ypos, " x/y line search errors ", xerr, yerr
 

          return xpos, ypos

     # shift - a 2D version of numpy's roll
     def array_shift(self, array,xshift=0,yshift=0):
          array = np.roll(array,xshift,0)
          array = np.roll(array,yshift,1)
          return array

     ## make an array with a circle set to one                 
     def circle(self, nx, ny, rad=None, cenx=None, ceny=None, invert=0 ):

          if rad is None: rad = np.min([nx,ny])/2
          if cenx is None: cenx = nx/2
          if ceny is None: ceny = ny/2

          print sys.version
          x = np.outer(np.arange(nx)-nx/2,np.ones(ny)) 
          y = np.outer(np.ones(nx),np.arange(ny)-ny/2) 

          dist = np.sqrt(x**2 + y**2)
          a = np.zeros([nx,ny])
          icirc = np.where(dist <= rad)
          a[icirc] = 1.
          
          if (invert==1):
               a = abs(1-a)

          out = self.array_shift(a, -nx/2, -ny/2)
          out = self.array_shift(out, cenx, ceny)
          return out

        
