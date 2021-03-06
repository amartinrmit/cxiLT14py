
"""
cxilt14py
- Analysis scripts for performing angular correlation analysis on the cspad detector @ CXI @ LCLS

params.py - defines common command line arguments

September 2018

authors: Andrew Martin (andrew.martin@rmit.edu.au)
         

"""


import argparse


class parameters:

    def __init__( self ):
        self.set_up_parser_with_common_args()

    #
    # Defines the par
    #
    def set_up_parser_with_common_args( self ):

        self.parser = argparse.ArgumentParser( description="cxi-lt14 analysis scripts for processing xtc files and performing angular correlation analysis." )


        self.parser.add_argument( "--datapath", "-d", help="Path to data directory", default="none" )
        self.parser.add_argument( "--exp", "-e", help="experiment code", default="cxiXXX")
        self.parser.add_argument( "--run", "-r", help="run number", default="1")
        self.parser.add_argument( "--cspadsrc", help="source code for the cspad detector", default='CXiDs1.0:Cspad.0')

        self.parser.add_argument( "--outpath", "-o", help="Path for output", default="none" )
        self.parser.add_argument( "--tag", "-t", help="A short string (tag) to name all output files.", default="test" )
        self.parser.add_argument( "--maskname", help="A filename (including path) of the mask image (= 0 for excluded pixels)", default="none" )

        self.parser.add_argument( "--nstart", help="Starting diffraction pattern number",type=int, default=0)
        self.parser.add_argument( "--nframes", "-n", help="Number of diffraction patterns to process",type=int, default=1e8)
        self.parser.add_argument( "--nthreads", help="Number of threads to use",type=int, default=1)
        self.parser.add_argument( "--verbose", help="more output to screen; 0 - none; >0 some/all", type=int, default=1)
        self.parser.add_argument( "--logname", help="Name of file to log input parameter", default="parameter_log.txt")

        self.parser.add_argument( "--cenx", help="Centre of the beam. x coordinate", type=int)
        self.parser.add_argument( "--ceny", help="Centre of the beam. y coordinate", type=int)
        self.parser.add_argument( "--outputfreq", help="Frequency that output is written.", type=int, default=100000)

       
    def parse_arguments( self ):
        self.args = self.parser.parse_args()

        
    #
    #  Writes all the input parameters to a log file
    #
    def write_all_params_to_file(self, name="None", script="None"):
        
        if name=="None":
            f = open( self.args.outpath+script[:-3]+"_"+self.args.tag+"_"+self.args.logname, 'w')
        else:
            f = open( name, 'w' )
        f.write( "# log of input parameters (cxilt14py)\n")

        if script != "None":
            f.write( "# generated by "+script+"\n")

        a = self.args.__dict__
        for d, e in a.items():
            f.write( d+" = "+str(e)+"\n" )

        f.close()
