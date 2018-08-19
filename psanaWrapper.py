"""
cxilt14py
- Analysis scripts for performing angular correlation analysis on the cspad detector @ CXI @ LCLS

psanaWrapper.py - main script for processing xtc files and performing angular correlation analysis

September 2018

authors: Andrew Martin (andrew.martin@rmit.edu.au)
         and others...

"""


import psana 
import numpy as np

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
        self.env_idx = self.ds_idx.env()
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
        wldet = psana.Detector( src, self.env_idx)
        return wldet(evt)


    def get_photon_beam_energy( self, evt, src='SIOC:SYS0:ML00:AO541'):
        phb = psana.Detector( src, self.env_idx )
        return phb(evt)

    def get_detector_z( self, evt, src='CXI:DS1:MMS:06.RBV'):
        dzp = psana.Detector( src, self.env_idx)
        return dzp(evt)

    def get_pulse_length( self, evt, src='SIOC:SYS0:ML00:AO820' ):
        pl = psana.Detector( src, self.env_idx )
        return pl(evt)

    def get_pulse_energy( self, evt, src='SIOC:SYS0:ML00:AO569' ):
        pl = psana.Detector( src, self.env_idx )
        return pl(evt)

    def qarrays( self, evt, cx=0.0, cy=0.0 ):
        
        dz = self.get_detector_z( evt )
        wavelength = self.get_wavelength( evt )
        q0 = 1/wavelength

        qx = self.cspad.coords_x( evt ) * q0 / (self.cspad.coords_z(evt) + dz )
        qy = self.cspad.coords_y( evt ) * q0 / (self.cspad.coords_z(evt) + dz )
        qz = q0

        qlen = np.sqrt( qx*qx + qy*qy + qz*qz )
        print "debug qarrays() q0 qlen", q0, np.min(qlen), np.max(qlen)

        qx = qx * q0 / qlen
        qy = qy * q0 / qlen
        qz = (qz *q0 / qlen) - q0
        qabs = np.sqrt( qx*qx + qy*qy + qz*qz )

        return [qx, qy, qz, qabs]

    #
    # convolve an assembled image. Account for mask
    # assumes convoltion struct is real and positive
    #
    def convolve_assembled_image( self, image, mask, struct ):
        
        imc = np.convolve( image, struct )
        maskc = np.convolve( mask, struct )
        imask = np.where( maskc > 0.0 )
        imc[imask] *= 1.0/maskc[imask]
        return imc

    #
    # average intensity on 2x1 asics
    # rebinx, rebiny - optional divide 2x1 asics into tiles and calculate average of times
    # e.g. rebinx=1, rebiny=2 calculates intensity on 1x1 asics
    #
    def asic_intensity( self, data, rebinx=1, rebiny=1 ):

        if rebin==1:
            out = np.sum( data, 3)
            out = np.sum( out, 2 )
        else:
            out = np.zeros( (data.shape[0], data.shape[1], rebinx, rebiny )) 
            istep = data.shape[2]/rebinx
            jstep = data.shape[3]/rebiny
            for i in np.arange(rebinx):
                for j in np.arange(rebiny):
                    tmp = data[:,:,i*istep:(i+1)*istep,j*jstep:(j+1)*jstep]
                    tmp = np.sum(  tmp, 3 )
                    tmp2 = np.sum( tmp, 2 )
                    out[:,:,i,j] = tmp2
        return out
