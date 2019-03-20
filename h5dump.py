"""
cxilt14py
- Analysis scripts for performing angular correlation analysis on the cspad detector @ CXI @ LCLS

h5dump.py - writes cspad frames from xtc file to h5

September 2018

authors: Andrew Martin (andrew.martin@rmit.edu.au)
         

"""


import numpy as np
import analysisTools as at
import h5py

#
# set up parsing of command line arguments
#
atp = at.params.parameters()

# add another command line argument
atp.parser.add_argument( "--outputext", help="File extention for output image: .dbin or .h5", default="dbin")
atp.parser.add_argument( "--average", help="Choose to output data average or data sum", type=bool, default=True)
atp.parser.add_argument( "--assemble",  help="assemble data (true) or leave as array of asic data", type=bool, default=False)
atp.parser.add_argument( "--raw", help="output uncorrected data (true) or all corrections [dark/gain/common mode]", type=bool, default=False)
atp.parser.add_argument( "--applymask", help="apply a mask retrieved from psana", type=bool, default=True)
atp.parser.add_argument( "--frame_no", help="Number of event to output to image file", type=int, default=0)

# parse the command line arguments
atp.parse_arguments()

# write all the input arguments to a log file
atp.write_all_params_to_file( script=atp.parser.prog )


#
# Set up the psana variables from psana wrapper
#
print atp.args.exp, atp.args.run
psbb = at.psanaWrapper.psanaBlackBox( exp=atp.args.exp, run=atp.args.run )


#
# retrieve mask if required
#
if atp.args.applymask == True:
    mask = psbb.cspad.mask( psbb.run.event( psbb.times[0]), calib=True, status=True, edges=True, central=True, unbondnbrs8=True)

#
# sum nframes of data from a run
#
datasum = psbb.cspad.calib( psbb.evt ) * 0.0

outname = atp.args.outpath+atp.parser.prog[:-3]+"_"+atp.args.tag+"_"+atp.args.exp+"_"+atp.args.run+"_"+str(atp.args.nstart)+"_"+str(atp.args.nframes)+".h5"

f = h5py.File(outname, 'w')    # overwrite any existing file



for i, t in enumerate( psbb.times, atp.args.nstart ):
    if i>=(atp.args.nframes+atp.args.nstart):
        break


    if atp.args.verbose >0:
        print "Processing event ", i
    evt = psbb.run.event(t)

    # sum the data
    if atp.args.raw == True:
        data = psbb.cspad.raw( evt )
    else:
        data = psbb.cspad.calib( evt, cmpars=(5,0,0,0) )




        #if atp.args.average == True:
        #    datasum *= 1.0/float(atp.args.nframes)
        
    if atp.args.applymask == True:
        datasum *= mask


    if atp.args.assemble == True:
        output = psbb.cspad.image( evt, datasum)
    else:
        output = datasum

    field = "/frames/frame_"+str(i)
    dset = f.create_dataset(field, data=output)

if atp.args.assemble == True:
    output = psbb.cspad.image( evt, mask)
else:
    output = mask
    
field = "/mask"
dset = f.create_dataset(field, data=output)

f.close()


