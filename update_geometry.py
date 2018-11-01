"""
I need a way to update the *.geom file
Then to convert this into a format that psana understands
"""
import h5py
import numpy as np
import shift_geom 
import analysisTools as at


pixsize = 109.91974263e-6


def fast_rad_prof(im, x, y, cx, cy):
    rs = np.round(
        np.sqrt((x-cx)**2 + (y-cy)**2), 0
        ).astype(np.uint16).ravel()
    
    r_count = np.bincount(rs)
    r_int   = np.bincount(rs, im.astype(np.float))
    return at.radial_profile.div_nonzero(r_int, r_count)

def find_offset(data, mask, x, y, rad_range = [0, 100], dx = 0.2, pix_range=10):
    N = pix_range
    cxs = np.linspace(-N/2 , N/2, N/dx)
    cys = np.linspace(-N/2 , N/2, N/dx)
    peak_heights = np.zeros((len(cxs), len(cys)), dtype= np.float)
    rad_avs = []
    min_rmax = np.inf
    
    r_mask = np.sqrt(x**2 + y**2)

    r_mask = (r_mask>rad_range[0]) * (r_mask<rad_range[1]) * mask.astype(np.bool)
    dm = data[r_mask]
    xm = x[r_mask]
    ym = y[r_mask]

    for i, cx in enumerate(cxs):
        for j, cy in enumerate(cys):
                
            rav = fast_rad_prof(dm, xm, ym, cx, cy)
            
            #rad_avs.append(rad_prof(datasum, mask, psbb, cx = cx, cy = cy))
            rad_avs.append(rav)
            if len(rav) < min_rmax :
                min_rmax = len(rav)
            
            #print i, j, cx, cy, rav.max()
    
    peak_heights = np.max([r[:min_rmax] for r in rad_avs], axis=-1).reshape(peak_heights.shape)

    i = np.argmax(peak_heights)
    i, j = np.unravel_index(i, peak_heights.shape)
    print('best cx, cy:', cxs[i]*pixsize, cys[j]*pixsize, 'um', 
                          cxs[i],         cys[j], 'pix')
    res = {'peak_heights': peak_heights, 'rad_avs': rad_avs}
    return cxs[i]*pixsize, cys[j]*pixsize, res

def ss_fs_to_ijkl(cspad_ij):
    """ 
    0: 388        388: 2 * 388  2*388: 3*388  3*388: 4*388
    (0, 0, :, :)  (1, 0, :, :)  (2, 0, :, :)  (3, 0, :, :)
    (0, 1, :, :)  (1, 1, :, :)  (2, 1, :, :)  (3, 1, :, :)
    (0, 2, :, :)  (1, 2, :, :)  (2, 2, :, :)  (3, 2, :, :)
    ...           ...           ...           ...
    (0, 7, :, :)  (1, 7, :, :)  (2, 7, :, :)  (3, 7, :, :)
    """
    cspad_psana_shape = (4, 8, 185, 388)
    cspad_geom_shape  = (1480, 1552)
    
    if cspad_ij.shape != cspad_geom_shape :
        raise ValueError('cspad input is not the required shape:' + str(cspad_geom_shape) )
    
    cspad_ijkl = np.zeros(cspad_psana_shape, dtype=cspad_ij.dtype)
    for i in range(4):
        cspad_ijkl[i] = cspad_ij[:, i * cspad_psana_shape[3]: (i+1) * cspad_psana_shape[3]].reshape((cspad_ijkl.shape[1:]))
    
    return cspad_ijkl

def get_pixel_map(fnam):
    import cfelpyutils.geometry_utils as gu
    import cfelpyutils.crystfel_utils as cu
    geom_dict = cu.load_crystfel_geometry(fnam)
    #pixel_map = gu.compute_visualization_pix_maps(geom_dict)
    pixel_map = gu.compute_pix_maps(geom_dict)
    x = ss_fs_to_ijkl(pixel_map.x).reshape((32, 185, 388))
    y = ss_fs_to_ijkl(pixel_map.y).reshape((32, 185, 388))
    return x, y

def get_pixel_map_psana(psbb):
    # load the x, y coordinates of each pixel
    x, y, z = psbb.cspad.coords_xyz(psbb.run.event( psbb.times[0]))

    # convert to pixel units
    x, y, z = x/(1e6*pixsize), y/(1e6*pixsize), z/(1e6*pixsize)
    return x, y

if __name__ == '__main__':
    fnam_powder = '/reg/neh/home/amorgan/2018/lt14/scratch/amorgan/powder/sum_1000_r0019.h5'
    fnam_powder = './nonecspad_sum_cxilt1417_40_nstart_0_nframes_2000_datasum.h5'
    #fnam_geom = '/reg/neh/home/amorgan/2018/lt14/scratch/amorgan/cxiLT14py/latest.geom'
    #fnam_geom = '/reg/neh/home/amorgan/2018/lt14/scratch/amorgan/cxiLT14py/cxilt14_r0019.geom'
    #fnam_geom_out = '/reg/neh/home/amorgan/2018/lt14/scratch/amorgan/cxiLT14py/cxilt14_r0019.geom'
    #fnam_geom_out = '/reg/neh/home/amorgan/2018/lt14/scratch/amorgan/cxiLT14py/cxilt14_r0019-test.geom'

    # load powder 
    powder = h5py.File(fnam_powder, 'r')['data'][()]
    mask = h5py.File(fnam_powder, 'r')['mask'][()]

    x, y = get_pixel_map_psana( 
            at.psanaWrapper.psanaBlackBox( exp='cxilt1417', run='40' )
            )

    cx, cy, res = find_offset(powder, mask, x, y, rad_range = [160, 190], pix_range=8)

    # now use shift geom
    #shift_geom.write_geom(fnam_geom, -cx, -cy, fnam_geom_out)
