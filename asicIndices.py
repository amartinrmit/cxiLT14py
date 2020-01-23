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

# add another command line argument
atp.parser.add_argument( "--outputext", help="File extention for output image: .dbin or .h5", default="jpg")
atp.parser.add_argument( "--average", help="Choose to output data average or data sum", type=bool, default=True)
atp.parser.add_argument( "--assemble",  help="assemble data (true) or leave as array of asic data", type=bool, default=True)
atp.parser.add_argument( "--raw", help="output uncorrected data (true) or all corrections [dark/gain/common mode]", type=bool, default=False)
atp.parser.add_argument( "--applymask", help="apply a mask retrieved from psana", type=bool, default=True)
atp.parser.add_argument( "--frame_no", help="Number of event to output to image file", type=int, default=0)

# parse the command line arguments
atp.parse_arguments()

atp.args.outpath = "/reg/d/psdm/cxi/cxilt1417/scratch/amartin/results/run40/aduhist/"
atp.args.exp = 'cxilt1417'
atp.args.run = '40'

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

i = atp.args.frame_no
t = psbb.times[i]    
if atp.args.verbose >0:
    print "Processing event ", i
evt = psbb.run.event(t)

msk = np.zeros( (32,185,388) )
for i in np.arange(32):
    msk *= 0.0
    msk += 0.2
    msk[i,:,:194] = np.outer( np.arange(185)/float(185), np.ones(194) )
    msk[i,:,194:] = np.outer( np.ones(185), np.arange(194)/float(194) )
    
    output = psbb.cspad.image( evt, msk)
    outname = atp.args.outpath+atp.parser.prog[:-3]+"_"+str(i)+".jpg"
    at.io.saveImage( atp.args, output, "asic_index", prog=atp.parser.prog, outname=outname )   
