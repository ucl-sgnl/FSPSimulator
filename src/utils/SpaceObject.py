### class for an individual space object
import uuid
import math 
import datetime
import numpy as np
import sgp4
from src.utils.coords import kep2car, trueanom2meananom, calculate_kozai_mean_motion

class Satellite:
    def __init__(self, cospar_id=None, rso_name=None, rso_type=None, payload_operational_status=None, orbit_type=None, application=None, source=None, 
                 orbital_status_code=None, launch_site=None, decay_date=None, mass=None, maneuverable=False, spin_stabilized=False, 
                 orbital_period=None, inclination=None, apogee_altitude=None, perigee_altitude=None, radar_cross_section=None, 
                 characteristic_area=None, characteristic_length=None, propulsion_type=None, epoch=None, sma=None, inc=None, 
                 argp=None, raan=None, tran=None, ecc=None):
        
        self.CATID = None # This is a computed property based on the GUID and isn't included as a parameter
        self.GUID = self._compute_catid() # This will be set by the system and isn't included as a parameter
        self.cospar_id = cospar_id
        self.rso_name = rso_name
        self.rso_type = rso_type
        self.payload_operational_status = payload_operational_status
        self.orbit_type = orbit_type
        self.application = application
        self.source = source # http://celestrak.org/satcat/sources.php
        self.orbital_status_code = orbital_status_code
        self.launch_site = launch_site
        self.decay_date = decay_date
        self.mass = mass #in Kg
        self.maneuverable = maneuverable
        self.spin_stabilized = spin_stabilized
        self.orbital_period = orbital_period #in minutes
        self.apogee_altitude = apogee_altitude #in Km
        self.perigee_altitude = perigee_altitude #in Km
        self.radar_cross_section = radar_cross_section #in meters^2
        self.characteristic_area = characteristic_area #in meters^2
        self.characteristic_length = characteristic_length #in meters
        self.propulsion_type = propulsion_type
        self.epoch = epoch
        self.sma = sma #in km
        self.inc = inc
        self.argp = argp
        self.raan = raan
        self.tran = tran
        self.eccentricity = ecc
        self.meananomaly = None
        self.cart_state = None #cartesian state vector [x,y,z,u,v,w] to be computed using generate_cart (from keplerian elements)
        self.C_d = 2.2 #Drag coefficient


    # This is a private method because its name starts with an underscore
    def _compute_catid(self):
        return uuid.uuid4()
    
    def generate_cart(self):
        x,y,z,u,v,w = kep2car(a = self.sma, e=self.eccentricity, i = self.inc, w = self.argp, W=self.raan, V=self.tran)
        self.cart_state = np.array([x,y,z,u,v,w])

    def sgp4_prop(self, jd):
        #jd is julian date of the epoch we want to propagate to
        #NOTE: from sgp4 docs
        # Note that ndot and nddot are ignored by the SGP4 propagator, 
        # so you can leave them 0.0 without any effect on the resulting satellite positions. 
        # But they do at least get saved to the satellite object, and written out if you write 
        # the parameters to a TLE or OMM file (see the “Export” section, above).
        sgp4_object = sgp4.api.Satrec()
        density = self.get_density()
        self.bstar = (self.C_d * self.characteristic_area * density)/2*self.mass #BStar = Cd * A * rho / 2m. Where Cd is the drag coefficient, A is the cross-sectional area of the satellite, rho is the density of the atmosphere, and m is the mass of the satellite.
        self.ndot = 0.0
        self.nddot = 0.0
        self.meananomaly = trueanom2meananom(self.tran, self.eccentricity)
        self.no_kozai = calculate_kozai_mean_motion(a = self.sma, mu = 398600.4418)
        sgp4_object.sgp4init('WGS72', 'i', self.satnum, self.epoch, self.bstar, self.ndot, self.nddot, self.eccentricity, self.argp, self.inc, self.meananomaly, self.no_kozai, self.raan)
        fr = 0.0
        error, position, velocity = sgp4_object.sgp4(jd, fr)
        if error != 0:
            print("Error in SGP4 propagation")
            return None
        else :
            self.cart_state = np.array(position + velocity)

    def get_density(self):
        print("fetching density from X Model")
        return 1e-12 #Placeholder value #in kg/m^3