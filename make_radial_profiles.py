import numpy as np
import analysisTools as at


#
# set up parsing of command line arguments
#
atp = at.params.parameters()


# parse the command line arguments
atp.parser.add_argument( "--outputext", help="File extention for output image: .dbin or .h5", default="dbin")
atp.parse_arguments()

#
# Set up the psana variables from psana wrapper
#
print atp.args.exp, atp.args.run
psbb = at.psanaWrapper.psanaBlackBox( exp=atp.args.exp, run=atp.args.run )

# load the mask
mask = psbb.cspad.mask( psbb.run.event(psbb.times[0]), calib=True, status=True, edges=True, central=True, unbondnbrs8=True)

# load the x, y coordinates of each pixel
x, y, z = psbb.cspad.coords_xyz(psbb.run.event( psbb.times[0]))

# convert to pixel units
x, y, z = x/109.91974263, y/109.91974263, z/109.91974263


rad_avs = []

for i, t in enumerate( psbb.times, atp.args.nstart ):
    if i>=(atp.args.nframes+atp.args.nstart):
        break

    if atp.args.verbose >0:
        print "Processing event ", i
    evt = psbb.run.event(t)

    image = psbb.cspad.calib(evt)
    
    rad_avs.append(at.radial_profile.make_radial_profile(image, x, y, mask))
    

rad_avs = np.array(rad_avs)

at.io.saveImage( atp.args, rad_avs, "radial_profiles", prog=atp.parser.prog )   
