"""
cxilt14py
- Analysis scripts for performing angular correlation analysis on the cspad detector @ CXI @ LCLS

getPressureData.py - simple script to get pressure values

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
atp.parser.add_argument( "--outputext", help="File extention for output image: .dbin or .h5", default="dbin")
atp.parser.add_argument( "--reprate", help="rep rate of FEL", type=int, default=120)

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
# load pressure detector
#
psbb.loadPressureDet()
#help( psbb.press )

#
# sum nframes of data from a run
#
reprate = atp.args.reprate
pvals = []
atp.args.nframes = np.min([atp.args.nframes,len(psbb.times)])

#for i, t in enumerate( psbb.times, atp.args.nstart ):
for i in np.arange(atp.args.nstart,np.min([atp.args.nframes,len(psbb.times)]),reprate):
    if i>=(atp.args.nframes+atp.args.nstart):
        break
    if i > len(psbb.times):
        break

    t = psbb.times[i]

    if atp.args.verbose >0:
        print "Processing event ", i

    evt = psbb.run.event(t)
    pval = psbb.get_pressure( evt )

    evtId = evt.get(at.psanaWrapper.psana.EventId) #EventId contains useful info try "print evtId" if you like.

    seconds = evtId.time()[0] # get event time stamp (seconds)
    print "evtID, seconds, pressure :", evtId, seconds, pval
    pvals.append( pval )
pvals = np.array( pvals )

fname = at.io.formatted_filename( atp.args, "pressure", "txt", prog=atp.parser.prog )
np.savetxt( fname, pvals )
