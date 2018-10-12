"""
cxilt14py
- Analysis scripts for performing angular correlation analysis on the cspad detector @ CXI @ LCLS

This is just for convience if you don't want to type a lot of arguments intot the command line.
If this script is useful please copy it, rename it and modify it, but don't update it in the git repository 

September 2018

authors: Andrew Martin (andrew.martin@rmit.edu.au)
         

"""


from os import system


script = "allpixelSVD.py"
params = {}

# common parameters (all parameter values are given as strings)
#params["--outpath"] = "/reg/neh/home/amartin/2018/cxi10016/results/run7/"
#params["--exp"] = 'cxip10016'
params["--outpath"] = "/reg/d/psdm/cxi/cxilt1417/scratch/amartin/results/run40/polarsvd/"
params["--tag"] = 'polarsvd'
params["--exp"] = 'cxilt1417'
params["--run"] = '40'
params["--nstart"] = '2000'
params["--nframes"] = '40000'
params["--verbose"] = '1'
#params["--rebinx"] = '1'
#params["--rebiny"] = '1'

# script specific parameters
params["--outputext"] = "dbin"
#params["--raw"] = '0'
params["--normalize"] = '1'
params["--diffCorr"] = '0'
params["--pcrange"] = "0 32 0 32"
params["--minI"] = "6.0e6"
params["--useMinI"] = "True"
params["--rankmax"] = "30"
params["--nq"] = "50"
params["--nth"] = "50"
params["--polarRange"] = "0 850 0 6.283185"
params["--cenx"] = "0"
params["--ceny"] = "0"


# generate the command line string from the dictionary of arguments
command = "python "+script+" "
for d, e in params.items():
    command += d+" "+e+" "

# run the script
print command
system( command )
