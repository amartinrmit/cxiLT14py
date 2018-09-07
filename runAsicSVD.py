"""
cxilt14py
- Analysis scripts for performing angular correlation analysis on the cspad detector @ CXI @ LCLS

This is just for convience if you don't want to type a lot of arguments intot the command line.
If this script is useful please copy it, rename it and modify it, but don't update it in the git repository 

September 2018

authors: Andrew Martin (andrew.martin@rmit.edu.au)
         

"""


from os import system


script = "subAsicSVD.py"
params = {}

# common parameters (all parameter values are given as strings)
#params["--outpath"] = "/reg/neh/home/amartin/2018/cxi10016/results/run7/"
#params["--exp"] = 'cxip10016'
params["--outpath"] = "/reg/neh/home/amartin/2018/cxic00217/results/run32/"
params["--tag"] = 'diffCorr_rebin12'
params["--exp"] = 'cxic00217'
params["--run"] = '32'
params["--nstart"] = '2000'
params["--nframes"] = '10000'
params["--verbose"] = '1'
params["--rebinx"] = '1'
params["--rebiny"] = '2'

# script specific parameters
params["--outputext"] = "dbin"
#params["--raw"] = '0'
params["--normalize"] = '1'
params["--diffCorr"] = '0'
params["--pcrange"] = "0 32 0 32"
params["--minI"] = "6.0e6"
params["--useMinI"] = "True"

# generate the command line string from the dictionary of arguments
command = "python "+script+" "
for d, e in params.items():
    command += d+" "+e+" "

# run the script
print command
system( command )
