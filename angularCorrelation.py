"""
cxilt14py
- Analysis scripts for performing angular correlation analysis on the cspad detector @ CXI @ LCLS

angularCorrelation.py - simple script to calculate angular correlations and check convergence

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

atp.parser.add_argument( "--normalize", help="normalize the data by the average pixels intensity", type=bool, default=True)
atp.parser.add_argument( "--diffCorr", help="Calculate the correlation of intenisty difference between current frame and a random frame.", type=bool, default=True)
atp.parser.add_argument( "--randomXcorr", help="Calculate the correlation between current frame and a random frame.", type=bool, default=False)
atp.parser.add_argument( "--normalize", help="normalize the data by the average pixels intensity", type=bool, default=True)

atp.parser.add_argument( "--pcrange", nargs=4, help="pixel range to evaluate the similarity (pearson correlation) of odd/even frame angular correlation functions ", type=int)

atp.parser.add_argument( "--polarRange", nargs=4, help="pixel range to evaluate polar plots of assembled images ", type=int)
atp.parser.add_argument( "--nq", help="number of radial (q) polar bins ", type=int)
atp.parser.add_argument( "--nth", help="number of angular polar bins ", type=int)

atp.parser.add_argument( "--cenx", help="beam centre in pixels from middle of assmebled array - xcoord ", type=int)
atp.parser.add_argument( "--ceny", help="beam centre in pixels from middle of assmebled array - ycoord ", type=int)


# parse the command line arguments
atp.parse_arguments()

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
# retrieve mask if required
#
if atp.args.applymask == True:
    mask = psbb.cspad.mask( psbb.run.event( psbb.times[0]), calib=True, status=True, edges=True, central=True, unbondnbrs8=True)

# copy some parameters for polar plot
nq, nth = atp.args.nq, atp.args.nth
rmin, rmax, thmin, thmax = atp.args.polarRange[0], atp.args.polarRange[1], atp.args.polarRange[2], atp.args.polarRange[3]
cenx, ceny = atp.args.cenx, atp.args.ceny

#
# sum nframes of data from a run
#
datasum = psbb.cspad.calib( psbb.evt ) * 0.0
corrqqsum = np.zeros( (nr,nth,m) )
if atp.args.randXcorr == True:
    corrqqsumX = np.zeros( (nr,nth,m) )
   

for i, t in enumerate( psbb.times, atp.args.nstart ):
    if i>=(atp.args.nframes+atp.args.nstart):
        break

    oddEven = i % 2

    if atp.args.verbose >0:
        print "Processing event ", i
    evt = psbb.run.event(t)

    j = np.int( np.random.rand()*len(psbb.times) )
    evtj = psbb.run.event(psbb.times[j])

    # sum the data
    if atp.args.raw == True:
        data = psbb.cspad.raw( evt )
        data2 = psbb.cspad.raw( evtj )
    else:
        data = psbb.cspad.calib( evt )
        data2 = psbb.cspad.calib( evtj )

    if atp.args.normalize==True:
        data *= 1.0/np.average(data*mask)
        data2 *= 1.0/np.average(data2*mask)

    if atp.args.diffCorr == True:
        diff = data - data2
        img = psbb.cspad.image(evt,diff*mask)
    else:
        img = psbb.cspad.image(evt,data*mask)

    pplot = ac.polar_plot( img, nq, nth, qmin, qmax, thmin, thmax, cenx+img.shape[0]/2, ceny+img.shape[1]/2, submean=True )
    corrqq = ac.polarplot_angular_correlation( pplot )      
    corrqqsum[:,:,m] += corrqq


    if atp.args.randomXcorr == True:
        img = psbb.cspad.image(evt,data*mask)
        img2 = psbb.cspad.image(evt,data2*mask)
        pplot = ac.polar_plot( img, nq, nth, qmin, qmax, thmin, thmax, cenx+img.shape[0]/2, ceny+img.shape[1]/2, submean=True )
        pplot2 = ac.polar_plot( img2, nq, nth, qmin, qmax, thmin, thmax, cenx+img.shape[0]/2, ceny+img.shape[1]/2, submean=True )
        corrqq = ac.polarplot_angular_correlation( pplot, pplot2 )      
        corrqqsumX[:,:,m] += corrqq
    
#
# correct correlations by mask correlation
#
corrqqsum[:,:,0] = ac.mask_correction( corrqqsum[:,:,0], corrqq_mask )
corrqqsum[:,:,1] = ac.mask_correction( corrqqsum[:,:,1], corrqq_mask )
if atp.args.randomXcorr == True:
    corrqqsum[:,:,0] = ac.mask_correction( corrqqsumX[:,:,0], corrqq_mask) 
    corrqqsum[:,:,1] = ac.mask_correction( corrqqsumX[:,:,1], corrqq_mask )


#
# pearson correlation to measure similarity of odd/even frame angular correlations
#
pc = ac.pearsoncorrelation2D( corrqqsum[:,:,0], corrqqsum[:,:,1], lim=atp.args.pcrange )

print "Pearson value of odd/even frame angular correlation :", pc


#
# save some output
#
at.io.saveImage( atp.args, corrqqsum[:,:,0], "even_frame_angular_correlation", prog=at.parser.prog )   
at.io.saveImage( atp.args, corrqqsum[:,:,1], "odd_frame_angular_correlation", prog=at.parser.prog )   
at.io.saveImage( atp.args, np.sum(corrqqsum,2), "total_angular_correlation", prog=at.parser.prog )   

