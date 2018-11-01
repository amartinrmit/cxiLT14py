import numpy as np
import analysisTools as at

pixsize = 109.91974263

def rad_prof(im, mask, psbb, rpix = None, cx = 0, cy = 0):
    # load the x, y coordinates of each pixel
    x, y, z = psbb.cspad.coords_xyz(psbb.run.event( psbb.times[0]))
    
    # convert to pixel units
    x, y, z = (x-cx)/pixsize, (y-cy)/pixsize, z/pixsize
    
    rad_av, rpix = at.radial_profile.make_radial_profile(im, x, y, mask, rpix)
    return rad_av


def fast_rad_prof(im, x, y, cx, cy):
    rs = np.round(
        np.sqrt((x-cx)**2 + (y-cy)**2)/pixsize, 0
        ).astype(np.uint16).ravel()
    
    r_count = np.bincount(rs)
    r_int   = np.bincount(rs, im.astype(np.float))
    return at.radial_profile.div_nonzero(r_int, r_count)



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
#x, y, z = psbb.cspad.coords_xyz(psbb.run.event( psbb.times[0]))

# convert to pixel units
#x, y, z = x/109.91974263, y/109.91974263, z/109.91974263


datasum = psbb.cspad.calib( psbb.evt ) * 0.0

#rad_avs = []
#rs = None
for i, t in enumerate( psbb.times, atp.args.nstart ):
    if i>=(atp.args.nframes+atp.args.nstart):
        break

    if atp.args.verbose >0:
        print "Processing event ", i
    evt = psbb.run.event(t)

    # sum the data
    cspad = psbb.cspad.calib(evt)
    
    #rad_av, rs = at.radial_profile.make_radial_profile(image, x, y, mask, rs)
    #rad_avs.append(rad_av)
    
    datasum += cspad

import h5py
with h5py.File('powder_3D_'+atp.args.exp+'_r'+atp.args.run.zfill(4)+'.h5') as f:
    f['data'] = datasum
    
#import pyqtgraph as pg
#pg.plot(rad_prof(datasum, mask, psbb))

"""

N = 40
cxs = np.linspace(-3*pixsize, -1*pixsize, N)
cys = np.linspace(-4*pixsize, -2*pixsize, N)
peak_heights = np.zeros((len(cxs), len(cys)), dtype= np.float)
rad_avs = []
min_rmax = np.inf

# load the x, y coordinates of each pixel
x, y, z = psbb.cspad.coords_xyz(psbb.run.event( psbb.times[0]))

# convert to pixel units
#x, y, z = x/pixsize, y/pixsize, z/pixsize

r_mask = np.sqrt(x**2 + y**2)/pixsize
r_mask = (r_mask>400) * (r_mask<600) * mask.astype(np.bool)
dm = datasum[r_mask]
xm = x[r_mask]
ym = y[r_mask]

for i, cx in enumerate(cxs):
    for j, cy in enumerate(cys):
            
        rav = fast_rad_prof(dm, xm, ym, cx, cy)
        
        #rad_avs.append(rad_prof(datasum, mask, psbb, cx = cx, cy = cy))
        rad_avs.append(rav)
        if len(rav) < min_rmax :
            min_rmax = len(rav)
        
        print i, j, cx, cy, rav.max()

peak_heights = np.max([r[:min_rmax] for r in rad_avs], axis=-1).reshape(peak_heights.shape)

i = np.argmax(peak_heights)
i, j = np.unravel_index(i, peak_heights.shape)
print('best cx, cy:', cxs[i], cys[j], 'um', cxs[i]/pixsize, cys[j]/pixsize, 'pix')

#rad_avs = np.array(rad_avs)
#at.io.saveImage( atp.args, rad_avs, "radial_profiles", prog=atp.parser.prog )   
"""
