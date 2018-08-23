"""
cxilt14py
- Analysis scripts for performing angular correlation analysis on the cspad detector @ CXI @ LCLS

asicCorrelation.py - simple script to calculate asic correlations and check convergence

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

atp.parser.add_argument( "--normalize", help="normalize the data by the average pixels intensity", type=bool, default=True)
atp.parser.add_argument( "--diffCorr", help="Calculate the correlation of intenisty difference between current frame and a random frame.", type=bool, default=False)
atp.parser.add_argument( "--randomXcorr", help="Calculate the correlation between current frame and a random frame.", type=bool, default=False)

atp.parser.add_argument( "--pcrange", nargs=4, help="pixel range to evaluate the similarity (pearson correlation) of odd/even frame asic correlation functions ", type=int)

atp.parser.add_argument( "--minI", help="minimum total integrated intensity ", type=float, default=100.)

#atp.parser.add_argument( "--cenx", help="beam centre in pixels from middle of assembled array : xcoord ", type=int)
#atp.parser.add_argument( "--ceny", help="beam centre in pixels from middle of assembled array : ycoord ", type=int)

# parse the command line arguments
atp.parse_arguments()
print "diffCorr", atp.args.diffCorr
# write all the input arguments to a log file
atp.write_all_params_to_file( script=atp.parser.prog )


#
# Set up the psana variables from psana wrapper
#
print atp.args.exp, atp.args.run
psbb = at.psanaWrapper.psanaBlackBox( exp=atp.args.exp, run=atp.args.run )


# angular correlation struct
ac = at.correlation.angular_correlation()



#
# retrieve mask and calculate its asic correlation
#
mask = psbb.cspad.mask( psbb.run.event( psbb.times[0]), calib=True, status=True, edges=True, central=True, unbondnbrs8=True)
# calculate mask intensity per asic
maskasic = psbb.asic_intensity( mask )



#
# sum nframes of data from a run
#
s = psbb.cspad.calib( psbb.evt ).shape
datasum =  np.zeros( (s[0], s[1], s[2], 2 ))
asiccorrsum = np.zeros( (32,32,2) )
nprocessed = np.zeros( 2 )
if atp.args.randomXcorr == True:
    asiccorrsumX = np.zeros( (32,32,2) )
   
totalIList = []
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

    totalI = data.sum()
    totalI2 = data2.sum()
    if atp.args.verbose >0:
        print "total intenisty :", totalI, totalI2

    if (totalI < atp.args.minI) or (totalI2 < atp.args.minI):
        print "Skipping event. Integrated intensity too low"
        continue
    else:
        nprocessed[m] += 1
        totalIList.append( [totalI, totalI2] )

    if atp.args.normalize==True:
        data *= 1.0/np.average(data*mask)
        data2 *= 1.0/np.average(data2*mask)

    datasum[:,:,:,m] += data*mask

    if atp.args.diffCorr == True:
        diff = (data - data2)*mask
        d = psbb.asic_intensity( diff ) /maskasic
    else:
        d = psbb.asic_intensity( data*mask ) /maskasic

    # correlate asic integrated values
    asiccorr =  ac.allpixel_correlation( d, d )

    asiccorrsum[:,:,m] += asiccorr


    if atp.args.randomXcorr == True:
        d = psbb.asic_intensity( data*mask ) /maskasic
        d2 = psbb.asic_intensity( data2*mask )/maskasic
        asiccorrsumX[:,:,m] += ac.allpixel_correlation( d, d2 )

# divide by number of processed frames
asiccorrsum[:,:,0] *= 1.0/float(nprocessed[0])
asiccorrsum[:,:,1] *= 1.0/float(nprocessed[1])
datasum[:,:,:,0] *= 1.0/float(nprocessed[0])
datasum[:,:,:,1] *= 1.0/float(nprocessed[1])

print "Number of processed frames even/odd", nprocessed
totalIarr = np.array(totalIList)
print "Min max integrated intensity :", np.min(totalIarr), np.max(totalIarr)


# correct asic correlation with mean
d = psbb.asic_intensity( datasum[:,:,:,0]*mask ) /maskasic
asiccorrsum[:,:,0] += - ac.allpixel_correlation( d, d)

d = psbb.asic_intensity( datasum[:,:,:,1]*mask ) /maskasic
asiccorrsum[:,:,1] += - ac.allpixel_correlation( d, d)


#
# pearson correlation to measure similarity of odd/even frame asic correlations
#
pc = ac.pearsonCorrelation2D( asiccorrsum[:,:,0], asiccorrsum[:,:,1], lim=atp.args.pcrange )


if atp.args.verbose > 0:
    print "Pearson value of odd/even frame asic correlation :", pc

outname = at.io.formatted_filename( atp.args, "pearson_values", "txt", prog=atp.parser.prog )
f = open( outname, 'w' )
f.write( "Pearson value of odd/even frame asic correlation :"+str(pc)+"\n" ) 
f.close()


#
# save some output
#
lim = atp.args.pcrange
at.io.saveImage( atp.args, asiccorrsum[:,:,0], "even_frame_asic_correlation", prog=atp.parser.prog )   
at.io.saveImage( atp.args, asiccorrsum[:,:,1], "odd_frame_asic_correlation", prog=atp.parser.prog )   
at.io.saveImage( atp.args, np.sum(asiccorrsum,2), "total_asic_correlation", prog=atp.parser.prog )   
img = psbb.cspad.image(evt,datasum[:,:,:,0])
at.io.saveImage( atp.args, img, "datasum_even", prog=atp.parser.prog )   

img = psbb.cspad.image(evt,datasum[:,:,:,1])
at.io.saveImage( atp.args, img, "datasum_odd", prog=atp.parser.prog )   

#
# make a correlation map
#
for asic in np.arange(mask.shape[0]):
    cmap = mask*0.0
    for i in np.arange(mask.shape[0]):
        cmap[i,:,:] = np.real(asiccorrsum[asic,i,0])
        
    imgcmap = psbb.cspad.image(evt,cmap)
    at.io.saveImage( atp.args , imgcmap, "_even_asic_correlation_map"+str(asic), ext='jpg', prog=atp.parser.prog  )

