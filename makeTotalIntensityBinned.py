"""
cxilt14py
- Analysis scripts for performing angular correlation analysis on the cspad detector @ CXI @ LCLS

makeTotalIntensityBinned.py - make a histogram of total intensity on the 2D detector

September 2020

authors: Andrew Martin (andrew.martin@rmit.edu.au)
         

"""


import numpy as np
import analysisTools as at
import h5py

def histindex( x, xmin, xmax, xbins):
    ib = int((x - xmin)*xbins/(xmax-xmin))
    if ib>xbins-1: ib = xbins-1
    if ib<0 : ib = 0
    return ib


#
# set up parsing of command line arguments
#
atp = at.params.parameters()

# add another command line argument
atp.parser.add_argument( "--outputext", help="File extention for output image: .dbin or .h5", default="dbin")
atp.parser.add_argument( "--average", help="Choose to output data average or data sum", type=bool, default=True)
atp.parser.add_argument( "--raw", help="output uncorrected data (true) or all corrections [dark/gain/common mode]", type=int, default=0)
#atp.parser.add_argument( "--raw", help="output uncorrected data (true) or all corrections [dark/gain/common mode]", type=bool, default=False)
atp.parser.add_argument( "--applymask", help="apply a mask retrieved from psana", type=bool, default=True)

atp.parser.add_argument( "--imax", help="max value in intensity bins", type=float, default=1e8)
atp.parser.add_argument( "--imin", help="min value in intensity bins", type=float, default=1e6)
atp.parser.add_argument( "--ibins", help="number of intensity bins", type=int, default=1)


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
# set intensity binning flag
#
if atp.args.ibins>0:
    ibinflag = True
    ibins = atp.args.ibins
else:
    ibinflag = False
    ibins = 1

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
totalIhist = np.zeros( ibins )

# list of all bin indices
indexlists = []
for i in range( ibins ):
    indexlists.append([])

for i in np.arange(len(psbb.times)-atp.args.nstart)+atp.args.nstart:
    t = psbb.times[i]
    if i>=(atp.args.nframes+atp.args.nstart):
        break


    if atp.args.verbose >0:
        print "Processing event ", i
    evt = psbb.run.event(t)

    # sum the data
    print "raw", atp.args.raw
    if atp.args.raw == True:
        data = psbb.cspad.raw( evt )
    else:
        data = psbb.cspad.calib( evt )
        dataraw = psbb.cspad.raw( evt )

    totalI = np.sum(data)
    totalIraw = np.sum(dataraw)
    ii = histindex( totalI, imin, imax, ibins )
    totalIhist[ii] += 1



#
# Output the histogram indices indices
#
outname = atp.args.outpath+atp.args.exp+"_"+atp.args.tag+"_"+atp.args.run+"_nstart_"+str(atp.args.nstart)+"_nframes_"+str(atp.args.nframes)+"_intensityHistogram.h5"

f = h5py.File(outname, 'w')      
field = "/ibinflag_imin_imax_ibins"
f.create_dataset(field, data=np.array([ibinflag, atp.args.imin,atp.args.imax,ibins]) )
field = "/totalIhist"
f.create_dataset(field, data=totalIhist )
for i in np.arange(ibins):
    field = "/indices/"+"bin"+str(i)
    f.create_dataset(field, data=np.array(indexlists[i]) )
f.close()
