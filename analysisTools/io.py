"""
cxilt14py
- Analysis scripts for performing angular correlation analysis on the cspad detector @ CXI @ LCLS

io.py - reading and writing files

September 2018

authors: Andrew Martin (andrew.martin@rmit.edu.au)
         

"""

import numpy as np
import os
import array
import struct
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

#
# reads a binary file of floating point numbers (no headers)
#      arguments:
#           fname : name of file to read
#           double: use double precision (True), or single precision (False)
#           shape:  shape of final array
#           dim:    dimension of final array if shape is not set, assumes sidelengths are equal
#

def read_dbin( fname, swapbyteorder=0, double=True, shape=(), dim=1 ):

    size = os.path.getsize(fname)
    # print "sizes:", size, size/8
    if double==True:
        code = 'd'
    else:
        code = 'f'
    b = array.array(code)
    f = open(fname, "r")
    b.fromfile( f, size/8 )
    f.close();
    l = b.tolist()
    
    if dim > 1:
        if len(shape)==0:
            n = len(l)
            nx = int(round(n**(1.0/float(dim))))
            shape = tuple( dim*[nx])
        output = np.array(l).reshape( shape )        
    
    if swapbyteorder == 1: output = output.newbyteorder()
    return output

#
# writes an array to a binary file as single or double precision
#
def write_dbin( fname, data, double=True ):
    if double==True:
        code = 'd'
    else:
        code = 'f'
    
    f = open( fname, "wb")
    fmt='<'+code*data.size
    bin = struct.pack(fmt, *data.flatten()[:] )
    f.write( bin )
    f.close()

#
# reads data from a h5 file
#
def h5read(filename,field="data/data1"):
     h5file = h5py.File(filename,"r")
     #print field
     h5data = h5file[field]
     image = h5data[...]
     h5file.close()
     return image

#
# writes an array into a h5 file
#
def h5write(filename,data,field="data/data",wa='w'):
     f = h5py.File(filename, wa)    # overwrite any existing file
     dset = f.create_dataset(field, data=data)
     f.close()

def saveImage( args, array, tail, ext=None, prog=None, h5field="/image", outname="None" ):

    if ext == None:
        ext = args.outputext

    if outname == "None":
        #set file name
        outname = args.outpath
    
        if prog != None:
            outname += prog[:-3]+"_"
    
        outname += args.exp+"_"+args.tag+"_"+args.run+"_nstart_"+str(args.nstart)+"_nframes_"+str(args.nframes)+"_"+tail+"."+ext

    if ext == "dbin":
        write_dbin( outname, array )
    elif ext == "h5":
        h5write( outname, array, field=h5field )
    else:
    #    try:
        plt.imsave( outname, array)
    #    except:
      #      print "error: cxiLT14py; analysisTools.io.saveImage"
    #        print "Either unknown extension:", args.outputext
    #        print "OR unknown location :", outname

    print "Saved image :", outname

def saveImshow( args, array, tail, ext=None, prog=None, h5field="/image", outname="None",
                cmin=0.0, cmax=0.0):

    if (cmin==0.0) and (cmax==0.0):
        cmin, cmax = np.min(array), np.max(array)

    if ext == None:
        ext = args.outputext

    if outname == "None":
        #set file name
        outname = args.outpath
    
        if prog != None:
            outname += prog[:-3]+"_"
    
        outname += args.exp+"_"+args.tag+"_"+args.run+"_nstart_"+str(args.nstart)+"_nframes_"+str(args.nframes)+"_"+tail+"."+ext

    if ext == "dbin":
        write_dbin( outname, array )
    elif ext == "h5":
        h5write( outname, array, field=h5field )
    else:
        plt.figure()
        plt.imshow( array )
        plt.clim( [cmin,cmax] )
        plt.colorbar()
        plt.draw()
        plt.savefig( outname )
    #    try:
        #plt.imsave( outname, array)
    #    except:
      #      print "error: cxiLT14py; analysisTools.io.saveImage"
    #        print "Either unknown extension:", args.outputext
    #        print "OR unknown location :", outname

    print "Saved image :", outname

def formatted_filename( args, tail, ext, prog=None ):

    outname = args.outpath
    if prog != None:
        outname += prog[:-3]+"_"
    outname += args.exp+"_"+args.tag+"_"+args.run+"_nstart_"+str(args.nstart)+"_nframes_"+str(args.nframes)+"_"+tail+"."+ext
    return outname
