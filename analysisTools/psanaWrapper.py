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

    def load_gain(self, fnam = '/reg/data/ana04/cxi/cxilt1417/scratch/gain/gain_prelim.npy'):
        self.gain = np.load(fnam)

    def load_mask(self, fnam = '/reg/data/ana04/cxi/cxilt1417/scratch/masks/better_mask.npy'):
        if fnam is None :
            self.mask = self.cspad.mask( evt, calib=True, status=True, edges=True, central=True, unbondnbrs8=True)
        else :
            self.mask = np.load(fnam)

    def loadCspad( self, cspadsrc='CxiDs2.0:Cspad.0'):
        #
        # get a detector object
        #
        self.cspad = psana.Detector(cspadsrc, self.env_smd)
        
        self.pixel_size = 109.91974263e-6
        self.shape = self.cspad.shape(par=0)
        self.size  = self.cspad.size(par=0)
        self.ndim  = self.cspad.ndim(par=0)
        self.instrument = self.cspad.instrument()
        self.x, self.y  = None, None
        self.mask = None
        self.gain = None

    def get_wavelength( self, evt, src='SIOC:SYS0:ML00:AO192'):
        wldet = psana.Detector( src, self.env_idx)
        return wldet(evt)


    def get_photon_beam_energy( self, evt, src='SIOC:SYS0:ML00:AO541'):
        phb = psana.Detector( src, self.env_idx )
        return phb(evt)

    def get_detector_z( self, evt, src='CXI:DS2:MMS:06.RBV', offset = 0.57538):
        try:
            dzp = psana.Detector( src, self.env_idx)
            dz = dzp(evt) * 1e-3 + offset
        except:
            dz = None
            print "(psanaWrapper.py) No detector distance found."
                
        return dz

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

        if (rebinx==1) and (rebiny==1):
            out = np.sum( data, 2)
            out = np.sum( out, 1 )
        else:
            out = np.zeros( (data.shape[0], rebinx, rebiny )) 
            istep = data.shape[1]/rebinx
            jstep = data.shape[2]/rebiny
            for i in np.arange(rebinx):
                for j in np.arange(rebiny):
                    tmp = data[:,i*istep:(i+1)*istep,j*jstep:(j+1)*jstep]
                    tmp = np.sum(  tmp, 2 )
                    tmp2 = np.sum( tmp, 1 )
                    out[:,i,j] = tmp2
        return out

    def nda_from_asic_intensity( self, asicI, shape, rebinx=1, rebiny=1 ):

        out  = np.zeros( shape )
        for asic in np.arange(shape[0]):
            if (rebinx==1) and (rebiny==1):
                out[asic,:,:] = asicI[asic]
            else:
                istep = shape[1]/rebinx
                jstep = shape[2]/rebiny
                for i in np.arange(rebinx):
                    for j in np.arange(rebiny):
                        out[asic,i*istep:(i+1)*istep,j*jstep:(j+1)*jstep] = asicI[asic,i,j]
        return out

    def cspad_calib(self, evt, polarisation_correction=True, gain_correction=True):
        """
        current recommended calibration pipleline. This is slow if done per event.
        """
        #
        # load the psana calibrated image with unbonded pixel based common mode correction
        #
        data = self.cspad.calib(evt, cmpars=(5, 0, 0, 0))

        #
        # Apply the mask
        #
        if self.mask is None :
            self.load_mask()

        data *= self.mask
        
        #
        # apply polarisation correction 
        #
        if polarisation_correction :
            #
            # get the pixel positions if not done already
            #
            if self.x is None or self.y is None :
                self.x, self.y, z = self.cspad.coords_xyz(evt)
                self.x *= 1e-6
                self.y *= 1e-6
             
            z = self.get_detector_z(evt)
            if z is not None :
                data = self.polarisation_correction(data, self.x, self.y, z)
        #
        # Apply the gain
        #
        if gain_correction :
            if self.gain is None :
                self.load_gain()
            data /= self.gain
        
        #
        # count photons
        # 
        data = self.cspad.photons(evt, data, self.mask, adu_per_photon=20.25)
        return data

    def polarisation_correction(self, cspad, x, y, z):
        pol = 1 - x**2 / (x**2 + y**2 + z**2)
        return cspad / pol

