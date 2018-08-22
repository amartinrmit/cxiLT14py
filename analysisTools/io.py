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
def h5write(filename,data,field="data/data"):
     f = h5py.File(filename, 'w')    # overwrite any existing file
     dset = f.create_dataset(field, data=data)
     f.close()
