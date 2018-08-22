"""
cxilt14py
- Analysis scripts for performing angular correlation analysis on the cspad detector @ CXI @ LCLS

This is just for convience if you don't want to type a lot of arguments intot the command line.
If this script is useful please copy it, rename it and modify it, but don't update it in the git repository 

September 2018

authors: Andrew Martin (andrew.martin@rmit.edu.au)
         

"""


from os import system


script = "sumdata.py"
params = {}

# common parameters (all parameter values are given as strings)
params["--datapath"] = "/reg/neh/home/amartin/2018/cxi10016/results/run7/"
params["--exp"] = 'cxip10016'
params["--run"] = '7'
params["--nstart"] = '0'
params["--nframes"] = '100'
params["--verbose"] = '1'

# script specific parameters
params["--outputext"] = "dbin"
params["--average"] = 'True'
params["--assemble"] = 'True'
params["--raw"] = 'False'
params["--applymask"] = 'True'


command = "python "+script+" "
for d, e in params.items():
    command += d+" "+e+" "

print command
system( command )
