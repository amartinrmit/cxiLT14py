"""
cxilt14py
- Analysis scripts for performing angular correlation analysis on the cspad detector @ CXI @ LCLS

asicSVD.py - calculate inter-asic SVD

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

atp.parser.add_argument( "--raw", help="output uncorrected data (true) or all corrections [dark/gain/common mode]", type=bool, default=False)

atp.parser.add_argument( "--normalize", help="normalize the data by the average pixels intensity", type=bool, default=False)
atp.parser.add_argument( "--diffCorr", help="Calculate the correlation of intenisty difference between current frame and a random frame.", type=bool, default=False)
atp.parser.add_argument( "--randomXcorr", help="Calculate the correlation between current frame and a random frame.", type=bool, default=False)

atp.parser.add_argument( "--pcrange", nargs=4, help="pixel range to evaluate the similarity (pearson correlation) of odd/even frame asic correlation functions ", type=int)

atp.parser.add_argument( "--useMinI", help="Reject frames below minimum intensity if true", type=bool )
atp.parser.add_argument( "--minI", help="minimum total integrated intensity ", type=float, default=100.)

#atp.parser.add_argument( "--cenx", help="beam centre in pixels from middle of assembled array : xcoord ", type=int)
#atp.parser.add_argument( "--ceny", help="beam centre in pixels from middle of assembled array : ycoord ", type=int)
atp.parser.add_argument( "--rebinx", help="bin asics to get sub asic intenisties : x ", type=int)
atp.parser.add_argument( "--rebiny", help="bin asics to get sub asic intensities : y", type=int)

# parse the command line arguments
atp.parse_arguments()
print "diffCorr", atp.args.diffCorr


#
# Set up the psana variables from psana wrapper
#
print atp.args.exp, atp.args.run
psbb = at.psanaWrapper.psanaBlackBox( exp=atp.args.exp, run=atp.args.run )
print "number of events in run :", len(psbb.times)

# check requested number of frames is not greater than total number of events
atp.args.nframes = np.min( [atp.args.nframes,len(psbb.times)] )

# write all the input arguments to a log file
atp.write_all_params_to_file( script=atp.parser.prog )



# angular correlation struct
ac = at.correlation.angular_correlation()



#
# retrieve mask and calculate its asic correlation
#
mask = psbb.cspad.mask( psbb.run.event( psbb.times[0]), calib=True, status=True, edges=True, central=True, unbondnbrs8=True)
# calculate mask intensity per asic
maskasic = psbb.asic_intensity( mask, atp.args.rebinx, atp.args.rebiny )
print mask.shape


#
# sum nframes of data from a run
#
s = psbb.cspad.calib( psbb.evt ).shape
datasum =  np.zeros( (s[0], s[1], s[2], 2 ))
nprocessed = np.zeros( 2 )

   
offdiagstats = []
cc_offdiagstats = []
totalIList = []
matrix = []
print "start loop"
for i, t in enumerate( psbb.times, atp.args.nstart ):
    if i>=(atp.args.nframes+atp.args.nstart):
        break

    m = i % 2

    if atp.args.verbose >0:
        print "Processing event ", i
    evt = psbb.run.event(t)

#    j = np.int( np.random.rand()*len(psbb.times) ) 
    j = np.int( np.random.rand()*atp.args.nframes ) + atp.args.nstart
    evtj = psbb.run.event(psbb.times[j])

    # sum the data
    if atp.args.raw == True:
        data = psbb.cspad.raw( evt )
        data2 = psbb.cspad.raw( evtj )
    else:
        data = psbb.cspad.calib( evt )
        data2 = psbb.cspad.calib( evtj )
#        data = psbb.cspad.raw( evt ) - psbb.cspad.pedestals( evt )
#        data2 = psbb.cspad.raw( evtj )- psbb.cspad.pedestals( evtj )
 

    totalI = data.sum()
    totalI2 = data2.sum()
    avIperpixel = (data*mask).mean()
    avIperpixel2 = (data2*mask).mean()
    if atp.args.verbose >0:
        print "total intenisty :", totalI, totalI2, atp.args.minI
        print "pixel average :", (data*mask).mean(), (data2*mask).mean()

    if atp.args.useMinI ==True:
        if (totalI < atp.args.minI) or (totalI2 < atp.args.minI):
            print "Skipping event. Integrated intensity too low"
            continue
        else:
            nprocessed[m] += 1
            totalIList.append( [totalI, totalI2] )
    else:
        nprocessed[m] += 1
        totalIList.append( [totalI, totalI2] )

    if atp.args.normalize==True:
        data *= 1.0/np.average(data*mask)
        data2 *= 1.0/np.average(data2*mask)

    datasum[:,:,:,m] += data*mask
    print "data size", data.size

    if atp.args.diffCorr == True:
        diff = (data - data2)*mask
        d = psbb.asic_intensity( diff, atp.args.rebinx, atp.args.rebiny ) /maskasic
    else:
        d = psbb.asic_intensity( data*mask, atp.args.rebinx, atp.args.rebiny ) /maskasic   #normalizes by number of pixels in an asic

    matrix.append( d.flatten() )

totalIarr = np.array( totalIList )

matrix = np.array( matrix )

nmat = np.dot( matrix.transpose(), matrix,  )

print nmat.shape, matrix.shape, totalIarr.shape
#exit()

u, s, vh = np.linalg.svd( nmat )

outname = at.io.formatted_filename( atp.args, "sing_vals", "txt", prog=atp.parser.prog )
np.savetxt( outname, s )
print "s", s[:20]
print "norm check", np.sum( u[:,0]*u[:,0] )
print "u ", u[:,:3]

outname = at.io.formatted_filename( atp.args, "eigenvectors", "txt", prog=atp.parser.prog )
np.savetxt( outname, u )


#
# make a eigenvector maps map
#
for asic in np.arange(mask.shape[0]):
    cmap = psbb.nda_from_asic_intensity( u[:,asic].reshape(32,atp.args.rebinx,atp.args.rebiny)\
                                         , mask.shape, rebinx=atp.args.rebinx,\
                                         rebiny=atp.args.rebiny )
    imgcmap = psbb.cspad.image(evt,cmap)
    at.io.saveImshow( atp.args , imgcmap, "_eigenvector"+str(asic), ext='jpg', prog=atp.parser.prog, cmin=-1.0, cmax=1.0  )

print np.dot( u[:,0], matrix.transpose() )

#
# project the data onto a basis vector and save the coefficients
#
asicmax = 10
for i in np.arange( asicmax ):
    outname = at.io.formatted_filename( atp.args, "e"+str(i)+"_coefficients", "txt", prog=atp.parser.prog )
    np.savetxt( outname, np.dot(u[:,i], matrix.transpose()) )

outname = at.io.formatted_filename( atp.args, "totalI", "txt", prog=atp.parser.prog )
np.savetxt( outname, totalIarr )
