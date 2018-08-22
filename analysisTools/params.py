
"""
cxilt14py
- Analysis scripts for performing angular correlation analysis on the cspad detector @ CXI @ LCLS

params.py - defines common command line arguments

September 2018

authors: Andrew Martin (andrew.martin@rmit.edu.au)
         

"""


import argparse


class parameters:


    #
    # Defines the par
    #
    def set_up_parser_with_common_args( self ):

        self.parser = argparse.ArgumentParser( description="cxi-lt14 analysis scripts for processing xtc files and performing angular correlation analysis." )


        self.parser.add_argument( "--datapath", "-d", nargs=1, help="Path to data directory", default="none" )
        self.parser.add_argument( "--exp", "-e", nargs=1, help="experiment code", default="cxiXXX")
        self.parser.add_argument( "--run", "-r", nargs=1, help="run number", default="1")
        self.parser.add_argument( "--cspadsrc", nargs=1, help="source code for the cspad detector", default='CXiDs1.0:Cspad.0')

        self.parser.add_argument( "--outpath", "-o", nargs=1, help="Path for output", default="none" )
        self.parser.add_argument( "--tag", "-t", nargs=1, help="A short string (tag) to name all output files.", default="test" )
        self.parser.add_argument( "--maskname", nargs=1, help="A filename (including path) of the mask image (= 0 for excluded pixels)", default="none" )

        self.parser.add_argument( "--nstart", nargs=1, help="Starting diffraction pattern number",type=float, default=0)
        self.parser.add_argument( "--npatterns", "-n", nargs=1, help="Number of diffraction patterns to process",type=int, default=1e8)
        self.parser.add_argument( "--nthreads", nargs=1, help="Number of threads to use",type=int, default=1)
	self.parser.add_argument( "--verbose", nargs=0, help="more output to screen", default=False)
	self.parser.add_argument( "--logname", nargs=1, help="Name of file to log input parameter", default="parameter_log.txt")

	self.parser.add_argument( "--cenx", nargs=1, help="Centre of the beam. x coordinate", type=float, default=0.)
	self.parser.add_argument( "--ceny", nargs=1, help="Centre of the beam. y coordinate", type=float, default=0.)
        self.parser.add_argument( "--outputfreq", nargs=1, help="Frequency that output is written.", type=int, default=100000)

       
    def parse_arguments( self ):
        self.args = self.parser.parse_args()

        
    #
    #  Writes all the input parameters to a log file
    #
    def write_all_params_to_file(self, name="None"):
        
        if name=="None":
            f = open( self.args.outpath+self.args.tag+self.args.logname, 'w')
        else:
            f = open( name, 'w' )

        a = self.args.__dict__
        for d, e in a.items():
            f.write( d, " = ", e )

        f.close()
