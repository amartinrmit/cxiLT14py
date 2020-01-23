"""
cxilt14py
- Analysis scripts for performing angular correlation analysis on the cspad detector @ CXI @ LCLS

makePowder.py - make the average powder plot from many frames

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
atp.parser.add_argument( "--average", help="Choose to output data average or data sum", type=bool, default=True)
atp.parser.add_argument( "--raw", help="output uncorrected data (true) or all corrections [dark/gain/common mode]", type=bool, default=False)
atp.parser.add_argument( "--applymask", help="apply a mask retrieved from psana", type=bool, default=True)
atp.parser.add_argument( "--nqbins", help="no. of q bins in powder plot.", type=int, default=500)
atp.parser.add_argument( "--qmin", help="minimum q value (default is lowest value in q arrays)", type=float, default=-1)
atp.parser.add_argument( "--qmax", help="maximum q value (default is highest value in q arrays)", type=float, default=-1)
atp.parser.add_argument( "--std", help="Calculate the standard deviation of each q ring", type=bool, default=True)
atp.parser.add_argument( "--cc", help="Calculate cross-correlation of each ring shifted by a fixed no of pixels (hard coded in powder.py)", type=bool, default=False)

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
# get arrays with q values
#
qarrays = psbb.qarrays( psbb.evt )

# default q range
if atp.args.qmin == -1:
    atp.args.qmin = np.min(qarrays[3])

if atp.args.qmax == -1:
    atp.args.qmax = np.max(qarrays[3])

#
# set up powder parameters
#
pp = at.powder.powder()
pp.set_up_q_indices( qarrays[3], nqbins=atp.args.nqbins, qmin=atp.args.qmin, qmax=atp.args.qmax )

#
# sum nframes of data from a run
#
datasum = psbb.cspad.calib( psbb.evt ) * 0.0
pplotsum = np.zeros( atp.args.nqbins )
pvarsum = np.zeros( atp.args.nqbins )
pccsum = np.zeros( atp.args.nqbins )

for i in np.arange(len(psbb.times)-atp.args.nstart)+atp.args.nstart:
    t = psbb.times[i]
    if i>=(atp.args.nframes+atp.args.nstart):
        break


    if atp.args.verbose >0:
        print "Processing event ", i
    evt = psbb.run.event(t)

    # sum the data
    if atp.args.raw == True:
        data = psbb.cspad.raw( evt )
    else:
        data = psbb.cspad.calib( evt )

    ppout = pp.powderplot_otf( data, mask, stdout=atp.args.std, ccout=False )
    j=0
    pplotsum += ppout[j]
    if atp.args.std == True:
        j+=1
        pvarsum += ppout[j]**2

    if atp.args.cc == True:
        j+=1
        pccsum += ppout[j]


if atp.args.average == True:
    pplotsum *= 1.0/float(atp.args.nframes)
    pvarsum *= 1./float(atp.args.nframes)
    pccsum  *= 1./float(atp.args.nframes)

pstd = np.sqrt(pvarsum )



outname = atp.args.outpath+atp.args.exp+"_"+atp.args.run+"_nstart_"+str(atp.args.nstart)+"_nframes_"+str(atp.args.nframes)+"_powder."+atp.args.outputext
output = pplotsum
if atp.args.outputext == "dbin":
    at.io.write_dbin( outname, output )
elif atp.args.outputext == "h5":
    at.io.h5write( outname, field="/datasum" )

if atp.args.std == True:
    outname = atp.args.outpath+atp.args.exp+"_"+atp.args.run+"_nstart_"+str(atp.args.nstart)+"_nframes_"+str(atp.args.nframes)+"_powder."+atp.args.outputext
    output = pstd
    if atp.args.outputext == "dbin":
        at.io.write_dbin( outname, output )
    elif atp.args.outputext == "h5":
        at.io.h5write( outname, field="/datasum" )

if atp.args.cc == True:
    outname = atp.args.outpath+atp.args.exp+"_"+atp.args.run+"_nstart_"+str(atp.args.nstart)+"_nframes_"+str(atp.args.nframes)+"_powder."+atp.args.outputext
    output = pccsum
    if atp.args.outputext == "dbin":
        at.io.write_dbin( outname, output )
    elif atp.args.outputext == "h5":
        at.io.h5write( outname, field="/datasum" )
