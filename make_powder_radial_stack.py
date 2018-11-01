import numpy as np
import analysisTools as at
from analysisTools import correlation 

def dilate(mask):
    from scipy.ndimage.morphology import binary_dilation
    mask = ~binary_dilation(~mask)
    return mask

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

evt = psbb.run.event(psbb.times[0])
d0  = psbb.cspad_calib(evt)

ang = correlation.angular_correlation()
num_radial_steps, num_phi_steps = 512, 512

# assembled mask 
mask_im = psbb.cspad.image(evt, psbb.mask)

# polar plot of mask 
pplot_mask = ang.polar_plot(
        mask_im, 
        num_radial_steps, num_phi_steps, 
        0, mask_im.shape[0]//2+200, 
        0, 2*np.pi, 
        mask_im.shape[0]//2,
        mask_im.shape[1]//2,
        submean=False
)
pplot_mask  = dilate(pplot_mask >= 1.).astype(np.float)

# angular correlation of mask
ac_mask = ang.polarplot_angular_correlation(pplot_mask).real
m = ac_mask > 1.


data_sum = d0 * 0.0
ac_sum = np.zeros(ac_mask.shape, dtype=np.float)
imgs = []
#for i, t in enumerate( psbb.times, atp.args.nstart ):
for i in range(atp.args.nframes):
    if atp.args.verbose >0:
        print "Processing event ", i
    evt = psbb.run.event(psbb.times[i])

    d = psbb.cspad_calib(evt)
    im = psbb.cspad.image(evt, d)
    
    pplot = ang.polar_plot(
            im, 
            num_radial_steps, num_phi_steps, 
            0, im.shape[0]//2+200, 
            0, 2*np.pi, 
            im.shape[0]//2,
            im.shape[1]//2,
            submean=False
    ) * pplot_mask
        
    # radial average
    sum_mask   = np.sum(pplot_mask, axis=1)
    mm         = sum_mask >= 1.
    rad_av     = np.zeros((pplot.shape[0],))
    rad_av[mm] = np.sum(pplot, axis = 1)[mm] / sum_mask[mm] 
        
    # subtract mean
    pplot = pplot_mask * ang.polar_plot_subtract_rmean(
            pplot, pplot_mask = pplot_mask
            )
    
    # make the angular correlations
    ac_image = ang.polarplot_angular_correlation(pplot).real

    ac_sum += ac_image
    data_sum += d
    

# normalise by the mask
ac_sum_norm = ac_sum.copy()
ac_sum_norm[m] = ac_sum_norm[m] / ac_mask[m]
ac_sum_norm *= m

image_sum = psbb.cspad.image(evt, data_sum)

# take the autocorrelation of the total
pplot = ang.polar_plot(
        image_sum, 
        num_radial_steps, num_phi_steps, 
        0, im.shape[0]//2+200, 
        0, 2*np.pi, 
        im.shape[0]//2,
        im.shape[1]//2,
        submean=False
) * pplot_mask
    
# radial average
sum_mask   = np.sum(pplot_mask, axis=1)
mm         = sum_mask >= 1.
rad_av     = np.zeros((pplot.shape[0],))
rad_av[mm] = np.sum(pplot, axis = 1)[mm] / sum_mask[mm] 
    
# subtract mean
pplot = pplot_mask * ang.polar_plot_subtract_rmean(
        pplot, pplot_mask = pplot_mask
        )

# make the angular correlations
ac_of_tot = ang.polarplot_angular_correlation(pplot).real

# normalise by the mask
ac_of_tot_norm = ac_of_tot.copy()
ac_of_tot_norm[m] = ac_of_tot_norm[m] / ac_mask[m]
ac_of_tot_norm *= m

outname = 'acs_'+atp.args.exp+"_"+atp.args.run+"_nstart_"+str(atp.args.nstart)+"_nframes_"+str(atp.args.nframes)+"_datasum.h5"

import h5py
with h5py.File(outname) as f:
    if 'data' in f :
        del f['data']
    if 'image' in f :
        del f['image']
    if 'mask' in f :
        del f['mask']
    if 'numevents' in f :
        del f['numevents']
    f['data'] = data_sum
    f['image'] = image_sum
    f['mask'] = psbb.mask
    f['numevents'] = i
    f['ac_sum'] = ac_sum
    f['ac_of_tot'] = ac_of_tot
    f['ac_sum_norm'] = ac_sum
    f['ac_of_tot_norm'] = ac_of_tot
    f['ac_mask'] = ac_mask
#at.io.saveImage( atp.args, output, "datasum", prog=atp.parser.prog )   
