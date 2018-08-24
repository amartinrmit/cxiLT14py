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
   
offdiagstats = []
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

    #
    #  get some statistics of the off diagonal asic correlation terms
    #
    d = psbb.asic_intensity( datasum[:,:,:,0]*mask ) /maskasic
    c0, mac0, od0 = ac.process_asiccorrsum( asiccorrsum[:,:,0], d, maskasic, nprocessed[0] )

    d = psbb.asic_intensity( datasum[:,:,:,1]*mask ) /maskasic
    c1, mac1, od1 = ac.process_asiccorrsum( asiccorrsum[:,:,1], d, maskasic, nprocessed[1] )

    offdiagstats.append( [od0.mean(), od0.std(), od1.mean(), od1.std()])



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
mean_asic_corr0 = ac.allpixel_correlation( d, d)
asiccorrsum[:,:,0] += - mean_asic_corr0

d = psbb.asic_intensity( datasum[:,:,:,1]*mask ) /maskasic
mean_asic_corr1 = ac.allpixel_correlation( d, d)
asiccorrsum[:,:,1] += - mean_asic_corr1


#
# pearson correlation to measure similarity of odd/even frame asic correlations
#
offdiag0 = asiccorrsum[:,:,0]-np.diag(np.diag(asiccorrsum[:,:,0]))
offdiag1 = asiccorrsum[:,:,1]-np.diag(np.diag(asiccorrsum[:,:,1]))

pc = ac.pearsonCorrelation2D( offdiag0, offdiag1, lim=atp.args.pcrange )

print "mean std offdiag even:", offdiag0.mean(), offdiag0.std()
print "mean std offdiag odd:", offdiag1.mean(), offdiag1.std()


if atp.args.verbose > 0:
    print "Pearson value of odd/even frame asic correlation :", pc

outname = at.io.formatted_filename( atp.args, "pearson_values", "txt", prog=atp.parser.prog )
f = open( outname, 'w' )
f.write( "Pearson value of odd/even frame asic correlation :"+str(pc)+"\n" ) 
f.close()

#
# normalize by diagonal to get cross-correlation coefficients
#
diag0 = np.sqrt(np.diag( np.abs(asiccorrsum[:,:,0]) ))
diag0corr = ac.allpixel_correlation( diag0, diag0)

diag1 = np.sqrt(np.diag( np.abs(asiccorrsum[:,:,1]) ))
diag1corr = ac.allpixel_correlation( diag1, diag1)

ccmatrix_even = asiccorrsum[:,:,0] / diag0corr
ccmatrix_odd = asiccorrsum[:,:,1] / diag1corr

print "diag0 (even)", diag0
print "diag1 (odd)", diag1
print "ccmatrix even diag", np.diag(ccmatrix_even)
print "ccmatrix odd  diag", np.diag(ccmatrix_odd)

#
# save some output
#
lim = atp.args.pcrange
at.io.saveImage( atp.args, asiccorrsum[:,:,0], "even_frame_asic_correlation", prog=atp.parser.prog )   
at.io.saveImage( atp.args, asiccorrsum[:,:,1], "odd_frame_asic_correlation", prog=atp.parser.prog )   
at.io.saveImage( atp.args, np.sum(asiccorrsum,2), "total_asic_correlation", prog=atp.parser.prog )   

at.io.saveImage( atp.args, mean_asic_corr0, "even_frame_mean_asic_correlation", prog=atp.parser.prog )   
at.io.saveImage( atp.args, mean_asic_corr1, "odd_frame_mean_asic_correlation", prog=atp.parser.prog )   

at.io.saveImage( atp.args, ccmatrix_even, "even_frame_asic_cross_correlation_coefficients", prog=atp.parser.prog )   
at.io.saveImage( atp.args, ccmatrix_odd, "odd_frame_mean_asic_cross_correlation_coefficients", prog=atp.parser.prog )   



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

#
# save the offdiagonal stats
#
offdiagstats = np.array(offdiagstats)
outname = at.io.formatted_filename( atp.args, "offdiagstats", "txt", prog=atp.parser.prog )
np.savetxt( outname, offdiagstats )
print "save file :", outname
