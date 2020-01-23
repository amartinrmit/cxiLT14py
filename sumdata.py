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
atp.parser.add_argument( "--outputext", help="File extention for output image: .dbin or .h5", default="dbin")
atp.parser.add_argument( "--average", help="Choose to output data average or data sum", type=bool, default=True)
atp.parser.add_argument( "--assemble",  help="assemble data (true) or leave as array of asic data", type=bool, default=True)
atp.parser.add_argument( "--raw", help="output uncorrected data (true) or all corrections [dark/gain/common mode]", type=bool, default=False)
atp.parser.add_argument( "--applymask", help="apply a mask retrieved from psana", type=bool, default=True)
atp.parser.add_argument( "--minI", help="minimum intensity to add to sum", type=float, default=0.0)

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

fname = atp.args.outpath+"intensitydata.txt"
f = open( fname, "w" )
f.write( "#run = "+str(atp.args.run)+"\n")
f.close()

atp.args.nframes = np.min([atp.args.nframes,len(psbb.times)])

nframes = 0
for i in np.arange(atp.args.nframes) +atp.args.nstart:
#    if i>=(atp.args.nframes+atp.args.nstart):
#        break

    if atp.args.verbose >0:
        print "Processing event ", i
    evt = psbb.run.event(psbb.times[i])

    print "Pulse energy in mJ:", psbb.get_pulse_energy(evt)

    # sum the data
    if atp.args.raw == True:
        data = psbb.cspad.raw( evt )
    else:
        data = psbb.cspad.calib( evt )

    if data.sum() < atp.args.minI:
        print "skipping frame. intensity too low. (nframes="+str(nframes)+")"
        continue

    nframes += 1
    datasum += data
    print data.sum()
    f = open( fname, "a" )
    f.write( str(i)+","+str(data.sum())+"\n")
    f.close()

print "Total frames included in sum", nframes

if atp.args.average == True:
    datasum *= 1.0/float(nframes)

if atp.args.applymask == True:
    datasum *= mask


if atp.args.assemble == True:
    output = psbb.cspad.image( evt, datasum)
    maska = psbb.cspad.image( evt, mask)
else:
    output = datasum
    maska = mask



#outname = atp.args.outpath+atp.args.exp+"_"+atp.args.run+"_nstart_"+str(atp.args.nstart)+"_nframes_"+str(atp.args.nframes)+"_datasum."+atp.args.outputext
#if atp.args.outputext == "dbin":
#    at.io.write_dbin( outname, output )
#elif atp.args.outputext == "h5":
#    at.io.h5write( outname, output, field="/datasum" )
print output.shape
at.io.saveImage( atp.args, output, "datasum", prog=atp.parser.prog )   

if atp.args.applymask==True:
    at.io.saveImage( atp.args, maska, "mask", prog=atp.parser.prog )   
