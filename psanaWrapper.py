"""
cxilt14py
- Analysis scripts for performing angular correlation analysis on the cspad detector @ CXI @ LCLS

psanaWrapper.py - main script for processing xtc files and performing angular correlation analysis

September 2018

authors: Andrew Martin (andrew.martin@rmit.edu.au)
         and others...

"""


import psana 


class psanaBlackBox:

    def __init__(self, exp='cxip10016', run='1' ):  

        #
        # Initialize variables 
        #     
        self.dsname_smd = 'exp='+exp+':run='+run+':smd'
        self.ds_smd = psana.DataSource( self.dsname_smd )
        self.env_smd = self.ds_smd.env()
        self.evt = self.ds_smd.events().next()

        self.dsname_idx = 'exp='+exp+':run='+run+':idx'
        self.ds_idx = psana.DataSource( self.dsname_idx )

        #
        # get event times
        #
        self.run = self.ds_idx.runs().next()
        self.times = self.run.times()
        self.nevents = len(self.times)

        # load information
        self.loadCspad()
        self.wavelength = self.get_wavelength( self.evt )
        self.energy = self.get_photon_beam_energy( self.evt )
        self.dz = self.get_detector_z( self.evt )

    def loadCspad( self, cspadsrc='CxiDs2.0:Cspad.0'):
        #
        # get a detector object
        #
        self.cspad = psana.Detector(cspadsrc, self.env_smd)
        
        self.shape = self.cspad.shape(par=0)
        self.size  = self.cspad.size(par=0)
        self.ndim  = self.cspad.ndim(par=0)
        self.instrument = self.cspad.instrument()

    def get_wavelength( self, evt, src='SIOC:SYS0:ML00:AO192'):
        wldet = psana.Detector( src, self.env_smd)
        return wldet(evt)


    def get_photon_beam_energy( self, evt, src='SIOC:SYS0:ML00:AO541'):
        phb = psana.Detector( src, self.env_smd )
        return phb(evt)

    def get_detector_z( self, evt, src='CXI:DS1:MMS:06.RBV'):
        dzp = psana.Detector( src, self.env_smd)
        return dzp(evt)

    def qarrays( self, evt, cx=0.0, cy=0.0 ):
        
        dz = self.get_detector_z( evt )
        wavelength = self.get_wavelength( evt )
        q0 = 1/wavelength

        qx = self.cspad.coords_x( evt ) * q0 / (self.cspad.coords_z + dz )
        qy = self.cspad.coords_y( evt ) * q0 / (self.cspad.coords_z + dz )
        qz = q0
        qz += q0

        qlen = sqrt( qx*qx + qy*qy + qz*qz )
        
        qx = qx / qlen
        qy = qy / qlen
        qz = (qz / qlen) - q0
        qabs = sqrt( qx*qx + qy*qy + qz*qz )

        return [qx, qy, qz, qabs]
