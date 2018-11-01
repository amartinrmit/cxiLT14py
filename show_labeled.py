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

import h5py
label = 3
with h5py.File('/reg/data/ana04/cxi/cxilt1417/scratch/radial_profiles/radial_profiles_cxilt1417_run6.h5', 'r') as f :
    labels = np.where(f['labels'][()] == label)[0]

imgs = []
#for i, t in enumerate( psbb.times, atp.args.nstart ):
for i, l in enumerate(labels):
    if i > (atp.args.nstart+atp.args.nframes):
        break
    t = psbb.times[l]
    if atp.args.verbose >0:
        print "Processing event ", i
    
    evt = psbb.run.event(t)
    imgs.append(psbb.cspad.image(evt))
    
import pyqtgraph as pg
imgs = np.array(imgs)
#pg.show(imgs[:100])


a = np.hstack(tuple([imgs[i, 400: -400, 400:-400] for i in range(10)]))
pg.show(a.T)

