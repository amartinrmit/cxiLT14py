"""
cxilt14py
- Analysis scripts for performing angular correlation analysis on the cspad detector @ CXI @ LCLS

angularCorrelation.py - simple script to calculate angular correlations and check convergence

September 2018

authors: Andrew Martin (andrew.martin@rmit.edu.au)
         

"""


import numpy as np
import analysisTools as at
import time


script = "angularCorrelation3Dindexed.py"

#
# set up parsing of command line arguments
#
atp = at.params.parameters()

# add another command line argument
atp.parser.add_argument( "--outputext", help="File extention for output image: .dbin or .h5", default="dbin")

atp.parser.add_argument( "--raw", help="output uncorrected data (true) or all corrections [dark/gain/common mode]", type=bool, default=False)

atp.parser.add_argument( "--normalize", help="normalize the data by the average pixels intensity", type=bool, default=False)
atp.parser.add_argument( "--intensity_veto", help="skip frames with intensity lower than minI", type=bool, default=False)
atp.parser.add_argument( "--diffCorr", help="Calculate the correlation of intenisty difference between current frame and a random frame.", type=bool, default=False)
atp.parser.add_argument( "--randomXcorr", help="Calculate the correlation between current frame and a random frame.", type=bool, default=False)

atp.parser.add_argument( "--pcrange", nargs=4, help="pixel range to evaluate the similarity (pearson correlation) of odd/even frame angular correlation functions ", type=int)

atp.parser.add_argument( "--polarRange", nargs=4, help="pixel range to evaluate polar plots of assembled images ", type=float)
atp.parser.add_argument( "--nq", help="number of radial (q) polar bins ", type=int)
atp.parser.add_argument( "--nth", help="number of angular polar bins ", type=int)
atp.parser.add_argument( "--minI", help="minimum total integrated intensity ", type=float, default=1e5)

#atp.parser.add_argument( "--cenx", help="beam centre in pixels from middle of assembled array : xcoord ", type=int)
#atp.parser.add_argument( "--ceny", help="beam centre in pixels from middle of assembled array : ycoord ", type=int)


# parse the command line arguments
atp.parse_arguments()
print "raw", atp.args.raw
print "diffCorr", atp.args.diffCorr
# write all the input arguments to a log file
atp.write_all_params_to_file( script=atp.parser.prog )



params = {}
params = atp.params

nq = atp.args.nq
nth = atp.args.nth

#
# sum nframes of data from a run
#
s = psbb.cspad.calib( psbb.evt ).shape
datasum =  np.zeros( (s[0], s[1], s[2], 2 ))
pplotsum = np.zeros( (nq,nth,2) )
pplot_mean_sum = np.zeros( (nq,nth,2) )
corrqqsum = np.zeros( (nq,nq,nth,2) )
nprocessed = np.zeros( 2 )
if atp.args.randomXcorr == True:
    corrqqsumX = np.zeros( (nq,nq,nth,2) )

norm1sum, norm2sum = 0.0, 0.0
nprocessed1sum, nprocessed2sum = 0.0, 0.0



command = "python "+atp.args.script+" "
for d, e in params.items():
    command += d+" "+e+" "


#
# main batch loop
#
fbatch = open( batchlist, "r")
for i, fname in enumerate(fbatch):
    params["--indexfile"] = fname
    params["--tag"] = batchtag+"_"+str(i)
    sytem(command)

    # open the output files and add sum them

    outname = at.io.formatted_filename( atp.args, "intensity_norm", "txt", prog=atp.parser.prog )
    fnorm = open( outname, 'r' )
    bits = []
    for line in fnorm:
        bits.append( line.split(" ") )
    norm1 = bits[0][2]
    norm2 = bits[1][2]
    nprocessed1 = bits[2][2]
    nprocessed2 = bits[3][2]
    fnorm.close()

    norm1sum += norm1*nprocessed1
    norm2sum += norm2*nprocessed2
    nprocessed1sum += nprocessed1
    nprocessed2sum += nprocessed2

    ext = "dbin"
    tail ="even_frame_angular_correlation"
    fname = at.io.formatted_filename( atp.args, tail, prog=script )
    corrqqsum[:,:,:,0] += at.io.read_dbin( fname )*nprocessed1

    tail ="odd_frame_angular_correlation"
    fname = at.io.formatted_filename( atp.args, tail, prog=script )
    corrqqsum[:,:,:,1] += at.io.read_dbin( fname )*nprocessed2


    tail ="datasum_even"
    fname = at.io.formatted_filename( atp.args, tail, prog=script )
    datasum[:,:,0] = at.io.read_dbin( fname )

    tail ="datasum_odd"
    fname = at.io.formatted_filename( atp.args, tail, prog=script )
    sum[:,:,1] = at.io.read_dbin( fname )

    tail ="polarplotsum_even"
    fname = at.io.formatted_filename( atp.args, tail, prog=script )
    polarplotsum[:,:,0] += at.io.read_dbin( fname )

    tail ="polarplotsum_odd"
    fname = at.io.formatted_filename( atp.args, tail, prog=script )
    polarplotsum[:,:,1] += at.io.read_dbin( fname )

    
    tail ="polarplot_mean_sum_even"
    fname = at.io.formatted_filename( atp.args, tail, prog=script )
    pplot_mean_sum[:,:,0] = at.io.read_dbin( fname )

    tail ="polarplot_mean_sum_odd"
    fname = at.io.formatted_filename( atp.args, tail, prog=script )
    pplot_mean_sum[:,:,1] = at.io.read_dbin( fname )




# divide by number of processed frames
corrqqsum[:,:,:,0] *= 1.0/float(nprocessed1sum)
corrqqsum[:,:,:,1] *= 1.0/float(nprocessed2sum)
norm1sum *= 1.0/float(nprocessed1sum)
norm2sum *= 1.0/float(nprocessed2sum)


# output normalization data
#
outname = at.io.formatted_filename( atp.args, "intensity_norm", "txt", prog=atp.parser.prog )
fnorm = open( outname, 'w' )
fnorm.write( "norm1 :"+str(norm1sum)+"\n" ) 
fnorm.write("norm2 :"+str(norm2sum)+"\n" ) 
fnorm.write( "nprocessed even :"+str(nprocessed1sum)+"\n" ) 
fnorm.write("nprocessed odd :"+str(nprocessed2sum)+"\n" ) 
fnorm.close()




#
# save some output
#
lim = atp.args.pcrange
at.io.saveImage( atp.args, corrqqsum[:,:,:,0], "even_frame_angular_correlation", prog=atp.parser.prog )   
at.io.saveImage( atp.args, corrqqsum[:,:,:,1], "odd_frame_angular_correlation", prog=atp.parser.prog )   
at.io.saveImage( atp.args, np.sum(corrqqsum,3), "total_angular_correlation", prog=atp.parser.prog )   

img = psbb.cspad.image(evt,datasum[:,:,:,0])
at.io.saveImage( atp.args, img, "datasum_even", prog=atp.parser.prog )   
img = psbb.cspad.image(evt,datasum[:,:,:,1])
at.io.saveImage( atp.args, img, "datasum_odd", prog=atp.parser.prog )   

at.io.saveImage( atp.args, pplotsum[:,:,0], "polarplotsum_even", prog=atp.parser.prog )   
at.io.saveImage( atp.args, pplotsum[:,:,1], "polarplotsum_odd", prog=atp.parser.prog )   
at.io.saveImage( atp.args, pplot_mean_sum[:,:,0], "polarplot_mean_sum_even", prog=atp.parser.prog )   
at.io.saveImage( atp.args, pplot_mean_sum[:,:,1], "polarplot_mean_sum_odd", prog=atp.parser.prog )   

