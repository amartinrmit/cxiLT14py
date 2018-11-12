"""
Adapted from Derek's script by Andrew (Morgan)

From outside SLAC:
    ssh -X amorgan@pslogin.slac.stanford.edu
    ssh psana
    source /reg/g/psdm/etc/psconda.sh
    cd /reg/d/psdm/cxi/cxilt1417/scratch/amorgan/

Usage:
    - Single thread process 200 frames:
    $ python make_radial_profiles.py --exp=cxilt1417 --run=73 --outpath='/reg/neh/home/amorgan/2018/lt14/scratch/radial_profiles/test_waxs.73.h5' --verbose=1 --nframes=200

    - 4 cpu cores, process all frames:
    $ mpirun -np 4 python make_radial_profiles.py --exp=cxilt1417 --run=73 --outpath='/reg/neh/home/amorgan/2018/lt14/scratch/radial_profiles/test_waxs.73.h5' --verbose=1

    - 8 cpu cores, process all frames, submit to the batch que :
    $ bsub -J waxs.73 -o waxs.73.log -q psanaq -n 8 mpirun python make_radial_profiles.py --exp=cxilt1417 --run=73 --outpath='/reg/neh/home/amorgan/2018/lt14/scratch/radial_profiles/test_waxs.73.h5' --verbose=1

    - 8 cpu cores, process all frames, process runs 10-20, submit to the batch que :
    $ for r in {10..20};do bsub -J waxs.$r -o waxs.$r.log -q psanaq -n 8 mpirun python make_radial_profiles.py --exp=cxilt1417 --run=$r --outpath='/reg/neh/home/amorgan/2018/lt14/scratch/radial_profiles/test_waxs.$r.h5' --verbose=1;done
"""

import numpy as np
import analysisTools as at


#
# set up parsing of command line arguments
#
atp = at.params.parameters()

atp.parser.add_argument( "--gain", help="gain file in .npy format", default="/reg/data/ana04/cxi/cxilt1417/scratch/gain/gain_prelim.npy")
atp.parser.add_argument( "--mask", help="mask file in .npy format", default="/reg/data/ana04/cxi/cxilt1417/scratch/masks/better_mask.npy")

# parse the command line arguments
atp.parse_arguments()

#
# Set up the psana variables from psana wrapper
#
print atp.args.exp, atp.args.run
psbb = at.psanaWrapper.psanaBlackBox( exp=atp.args.exp, 
                                      run=atp.args.run, 
                                      gain_fnam = atp.args.gain,
                                      mask_fnam = atp.args.mask)

# convert to pixel units
x, y, z = psbb.x/psbb.pixel_size, psbb.y/psbb.pixel_size, psbb.z/psbb.pixel_size

events = psbb.ds_smd.events()

# use psana's outputter for MPI
smldata = psbb.ds_smd.small_data(atp.args.outpath, gather_interval=100)

rs = None
for i, evt in enumerate( events ) :
    if i>=atp.args.nframes:
        break

    if atp.args.verbose >0:
        print "Processing event ", i

    if evt is None :
        continue

    image = psbb.cspad_calib(evt)
    
    rad_av, rs = at.radial_profile.make_radial_profile(image, x, y, psbb.mask, rs)

    smldata.event( radialpro = rad_av) 

smldata.save()

