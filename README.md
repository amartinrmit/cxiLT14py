# cxi-lt14

Scripts for angular correlation analysis of x-ray diffraction data collected at CXI beamline at LCLS

The analysisTools directory contains algorithms and functions that are generally useful
e.g.
- setting up psana variables (psanaWrapper.py)
- calculating powder plots (powder.py)
- angular correlation analysis (correlation.py)
- reading and writing images (io.py)
- finding the beam centre (centreDiffraction.py)
- reading command line parameters (params.py)

The main directory contains scripts for specific analyis tasks
e.g.
- summing data in a run (sumdata.py)

run_a_script.py is a template for using using python to call an analysis script. It saves having to write lots of command line arguments into the terminal. Please copy and modify run_a_script.py, but don't push changes back to the repository unless you really want to change the template for everyone.

Command line parameters are managed by argparse. Common parameters are defined in params.py, but each analysis script can append further command line parameters (see sumdata.py)
