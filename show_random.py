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


imgs = []
np.random.seed(0)
#for i, t in enumerate( psbb.times, atp.args.nstart ):
for i in range(atp.args.nframes):
    i = np.random.randint(0, len(psbb.times))

    if atp.args.verbose >0:
        print "Processing event ", i
    evt = psbb.run.event(psbb.times[i])

    d = mask * psbb.cspad.calib(evt, cmpars=(5, 0, 0, 0))
    #d = psbb.cspad.calib(evt)
    imgs.append(psbb.cspad.image(evt, d))
    
import pyqtgraph as pg
pg.show(np.array(imgs))
