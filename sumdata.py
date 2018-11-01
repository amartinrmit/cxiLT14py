"""
cxilt14py
- Analysis scripts for performing angular correlation analysis on the cspad detector @ CXI @ LCLS

sumdata.py - simple script to sum the data in a run

September 2018

authors: Andrew Martin (andrew.martin@rmit.edu.au)
         

"""


import numpy as np
import analysisTools as at


#
# set up parsing of command line arguments
#
atp = at.params.parameters()

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
# sum nframes of data from a run
#
datasum = psbb.cspad.calib( psbb.evt ) * 0.0

for i, t in enumerate( psbb.times, atp.args.nstart ):
    if i>=(atp.args.nframes+atp.args.nstart):
        break

    if atp.args.verbose >0:
        print "Processing event ", i
    evt = psbb.run.event(t)

    # sum the data
    data = psbb.cspad_calib( evt )

    datasum += data

output = psbb.cspad.image( evt, datasum)

outname = atp.args.outpath+'cspad_sum_'+atp.args.exp+"_"+atp.args.run+"_nstart_"+str(atp.args.nstart)+"_nframes_"+str(atp.args.nframes)+"_datasum.h5"

image = psbb.cspad.image(evt, datasum)


import h5py
with h5py.File(outname) as f:
    if 'data' in f :
        del f['data']
    if 'image' in f :
        del f['image']
    if 'mask' in f :
        del f['mask']
    if 'numevents' in f :
        del f['numevents']
    f['data'] = datasum
    f['image'] = image
    f['mask'] = psbb.mask
    f['numevents'] = i
#at.io.saveImage( atp.args, output, "datasum", prog=atp.parser.prog )   
