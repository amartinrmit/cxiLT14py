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

    def __init__(self, exp='cxip10016', run='1', gain_fnam=None, mask_fnam=None ):  
        """
        By default gain = 1 and the mask is obtained from psana.
        For example:
           gain_fnam = '/reg/data/ana04/cxi/cxilt1417/scratch/gain/gain_prelim.npy' 
           mask_fnam = '/reg/data/ana04/cxi/cxilt1417/scratch/masks/better_mask.npy'
        """
        #
        # Initialize variables 
        #     
        print "Debug: loading psanabb"
        self.dsname_smd = 'exp='+exp+':run='+run+':smd'
        self.ds_smd = psana.MPIDataSource( self.dsname_smd )
        self.env_smd = self.ds_smd.env()
        self.evt = self.ds_smd.events().next()

        print "Debug: loading psanabb idx events"
        self.dsname_idx = 'exp='+exp+':run='+run+':idx'
        self.ds_idx = psana.DataSource( self.dsname_idx )
        self.env_idx = self.ds_idx.env()
        #
        # get event times
        #
        self.run = self.ds_idx.runs().next()
        self.times = self.run.times()
        self.nevents = len(self.times)

        print "Debug: loading cspad info"
        # load information
        self.loadCspad(gain_fnam=gain_fnam, mask_fnam=mask_fnam)
        self.wavelength = self.get_wavelength( self.evt )
        self.energy = self.get_photon_beam_energy( self.evt )
    
    def loadCspad( self, cspadsrc='CxiDs2.0:Cspad.0', gain_fnam=None, mask_fnam=None):
        #
        # get a detector object
        #
        self.cspad = psana.Detector(cspadsrc, self.env_smd)
        
        self.pixel_size = 109.91974263e-6
        self.shape = self.cspad.shape(par=0)
        self.size  = self.cspad.size(par=0)
        self.ndim  = self.cspad.ndim(par=0)
        self.instrument = self.cspad.instrument()
        
        # load x-y pixel coords and calculate polarisation
        self.x, self.y, z = self.cspad.coords_xyz(self.evt)
        self.x *= 1e-6
        self.y *= 1e-6
        self.z  = self.get_detector_z(self.evt)
        self.pol = 1 - self.x**2 / (self.x**2 + self.y**2 + self.z**2)

        print "DETECTOR Z DISTANCE :", self.z
        
        # load gain (such that gain correction: im = data / gain)
        if gain_fnam is not None :
            self.gain = np.load(gain_fnam)
        else :
            self.gain = 1.
        
        # load mask 
        if mask_fnam is not None :
            self.mask = np.load(mask_fnam)
        else :
            self.mask = self.cspad.mask( self.evt, calib=True, status=True, edges=True, central=True, unbondnbrs8=True)

    def cspad_calib(self, evt, common_mode=True, mask=True, gain=True, polarisation=True):
        """
        Return the corrected cspad data for event 'evt'.
        
        Applies calibration to the raw cspad data. This is on top of the 
        default psana calibration. 
        
        Parameters
        ----------
        evt : psana event object

        common_mode : bool
            if True then common mode correction using unbonded pixels.
        
        mask : bool
            if True then self.mask is multiplied.
        
        gain : bool
            if True then self.gain is divided.
        
        polarisation : bool
            if True then self.pol is divided.
        
        Returns
        -------
        cspad_data : numpy array
            Corrected cspad as a numpy array
        """
        if common_mode :
            cspad_data = self.cspad.calib(evt, cmpars=(5, 0, 0, 0))
        else :
            cspad_data = self.cspad.calib(evt)

        if mask :
            cspad_data *= self.mask

        if gain :
            cspad_data /= self.gain
            
        if polarisation :
            cspad_data /= self.pol
                
        return cspad_data

    def load_gain(self, fnam = '/reg/data/ana04/cxi/cxilt1417/scratch/gain/gain_prelim.npy'):
        self.gain = np.load(fnam)

    def loadPressureDet( self, src='CXI:LC20:SDS:Pressure'):
        self.press = psana.Detector( src,self.env_idx )

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

    def get_pressure( self, evt ):
        pval = self.press(evt)
        if pval == None:
            output = 0.0
        else:
            output = pval
        return output

    def qarrays( self, evt, cx=0.0, cy=0.0, dz=-1, wavelength=-1 ):
        
        if dz == -1:
            dz = self.get_detector_z( evt ) + (0.139-0.567)

        if wavelength == -1:
            wavelength = self.get_wavelength( evt )

        q0 = 1/wavelength

        qx = self.cspad.coords_x( evt ) * q0 / dz #(self.cspad.coords_z(evt) + dz )
        qy = self.cspad.coords_y( evt ) * q0 / dz #(self.cspad.coords_z(evt) + dz )
        qz = q0
        print "debug qx, qy", np.max(self.cspad.coords_x( evt )) , dz, np.max(self.cspad.coords_z(evt)), np.max(qx), q0

        qlen = np.sqrt( qx*qx + qy*qy + qz*qz )
        print "debug qarrays() q0 qlen", q0, np.min(qlen), np.max(qlen)

        qx = qx * q0 / qlen
        qy = qy * q0 / qlen
        qz = (qz *q0 / qlen) - q0
        qabs = np.sqrt( qx*qx + qy*qy + qz*qz )
        qangle = np.arctan2( qx, qy )
        print "debug qarrays() qabs max", np.max(qabs), np.max(qx ), np.max(qy), np.max(qz)
        
        return [qx, qy, qz, qabs, qangle]

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
