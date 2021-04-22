"""
cxilt14py
- Analysis scripts for performing angular correlation analysis on the cspad detector @ CXI @ LCLS

This is just for convience if you don't want to type a lot of arguments intot the command line.
If this script is useful please copy it, rename it and modify it, but don't update it in the git repository 

September 2018

authors: Andrew Martin (andrew.martin@rmit.edu.au)
         

"""


from os import system


script = "angularCorrelation3Dindexed.py"
params = {}

# common parameters (all parameter values are given as strings)
params["--outpath"] = "/cds/home/a/amartin/2018/cxilt1417/results/run96/corr/"
params["--tag"] = 'p340_400_ewald'
params["--exp"] = 'cxilt1417'
params["--run"] = '96'
params["--nstart"] = '0'
params["--nframes"] = '300000'
params["--verbose"] = '1'

params["--indexfile"] = "/cds/home/a/amartin/2018/cxilt1417/results/sorting/waxs.96_340_400_nstart0_indexlist.npy"

# script specific parameters
params["--outputext"] = "dbin"
#params["--raw"] = '0'
#params["--weighted"] = 'True'
#params["--normalize"] = 'True'
#params["--intensity_veto"] = 'True'
params["--nq"] = "100"
params["--nth"] = "64"
params["--polarRange"] = "0 850 0 6.283185"
params["--pcrange"] = "8 45 10 50"
params["--cenx"] = "-6"     #-6
params["--ceny"] = "12"     #12
params["--diffCorr"] = "True"
params["--dz"] = "0.139"
params["--wavelength"] = "1.305e-10"
#params["--svdfile"] = "/reg/d/psdm/cxi/cxilt1417/scratch/amartin/results/run40/polarsvd_newcen/allpixelSVD_cxilt1417_40_2000_20000_evenframe_svdmodes.h5"
#params["--svdnsub"] = "30"
#params["--rankmax"] = "30" 

# generate the command line string from the dictionary of arguments
command = "python "+script+" "
for d, e in params.items():
    command += d+" "+e+" "

# run the script
print command
system( command )
