### class for an individual space object
import uuid
import math 
import datetime
import numpy as np
import sgp4
import matplotlib.pyplot as plt
from sgp4.api import Satrec, WGS72
from src.utils.coords import kep2car, trueanom2meananom, calculate_kozai_mean_motion, expo_simplified, utc_to_jd
import matplotlib.cm as cm

class Satellite:
    def __init__(self, cospar_id=None, rso_name=None, rso_type=None, payload_operational_status=None, orbit_type=None, application=None, source=None, 
                 orbital_status_code=None, launch_site=None, decay_date=None, mass=None, maneuverable=False, spin_stabilized=False, 
                 orbital_period=None, inclination=None, apogee_altitude=None, perigee_altitude=None, radar_cross_section=None, 
                 characteristic_area=None, characteristic_length=None, propulsion_type=None, epoch=None, sma=None, inc=None, 
                 argp=None, raan=None, tran=None, ecc=None):
        
        self.CATID = None # This is a computed property based on the GUID and isn't included as a parameter
        self.GUID = self._compute_catid() # This will be set by the system and isn't included as a parameter
        self.cospar_id = str(cospar_id)
        self.rso_name = str(rso_name)
        self.rso_type = str(rso_type)
        self.payload_operational_status = payload_operational_status

        self.orbit_type = orbit_type
        self.application = application
        self.source = source # http://celestrak.org/satcat/sources.php
        self.orbital_status_code = str(orbital_status_code)
        self.launch_site = str(launch_site)
        self.decay_date = datetime.datetime.strptime(decay_date, '%Y-%m-%d %H:%M:%S') #in UTC
        self.mass = float(mass) #in Kg
        self.maneuverable = str(maneuverable)
        self.spin_stabilized = str(spin_stabilized)# I have left spin_stabilized and maneuverable as strings in case we later want to add more options than just True/False (for example different thruster types resulting in different kinds of maneuverability)
        self.orbital_period = float(orbital_period) #in minutes
        self.apogee_altitude = float(apogee_altitude) #in Km
        self.perigee_altitude = float(perigee_altitude) #in Km
        self.radar_cross_section = float(radar_cross_section) #in meters^2
        self.characteristic_area = float(characteristic_area) #in meters^2
        self.characteristic_length = float(characteristic_length) #in meters
        self.propulsion_type = propulsion_type
        #epoch must be cast to datetime object and be specified in UTC time in the format: datetime(year-month-day hour:minute:second)
        self.epoch = datetime.datetime.strptime(epoch, '%Y-%m-%d %H:%M:%S') #in UTC
        self.sma = sma #in km
        self.inc = inc
        self.argp = argp
        self.raan = raan
        self.tran = tran
        self.eccentricity = ecc
        self.meananomaly = None
        self.cart_state = None #cartesian state vector [x,y,z,u,v,w] to be computed using generate_cart (from keplerian elements)
        self.C_d = 2.2 #Drag coefficient

        #Value Checks
        possible_operational_status = ['+', '-', 'P', 'B', 'S', 'X', 'D', '?'] #https://celestrak.org/satcat/status.php
        if self.payload_operational_status not in possible_operational_status:
            raise ValueError('payload_operational_status must be one of the following: {}'.format(possible_operational_status))
        
        possible_spin_stabilized = ['y', 'n']
        if self.spin_stabilized not in possible_spin_stabilized:
            raise ValueError('spin_stabilized must be one of the following: {}'.format(possible_spin_stabilized))
        
        possible_maneuverable = ['y', 'n']
        if self.spin_stabilized not in possible_spin_stabilized:
            raise ValueError('spin_stabilized must be one of the following: {}'.format(possible_spin_stabilized))
        
        possible_launch_sites = ['CCSFS', 'KSC', 'VAFB', 'Kwajalein', 'Other'] #TODO: update with actual values
        if self.launch_site not in possible_launch_sites:
            raise ValueError('launch_site must be one of the following: {}'.format(possible_launch_sites))
        
        possible_orbital_status_codes = ['EA0']
        
        #check that inc, argp, raan and tran are all between 0 and 2pi
        if self.inc < 0 or self.inc > 2*math.pi:
            raise ValueError('inc must be in Radians. Current value not between 0 and 2pi')
        if self.argp < 0 or self.argp > 2*math.pi:
            raise ValueError('argp must be in Radians. Current value not between 0 and 2pi')
        if self.raan < 0 or self.raan > 2*math.pi:
            raise ValueError('raan must be in Radians. Current value not between 0 and 2pi')
        if self.tran < 0 or self.tran > 2*math.pi:
            raise ValueError('tran must be in Radians. Current value not between 0 and 2pi')
        
        if self.ecc<0 or self.ecc>1:
            raise ValueError('eccentricity must be between 0 and 1')

    # This is a private method because its name starts with an underscore
    def _compute_catid(self):
        return uuid.uuid4()
    
    def generate_cart(self):
        # generate cartesian state vector from keplerian elements
        # Keplerian elements are in radians
        x,y,z,u,v,w = kep2car(a = self.sma, e=self.eccentricity, i = self.inc, w = self.argp, W=self.raan, V=self.tran)
        self.cart_state = np.array([x,y,z,u,v,w])
        return

    def sgp4_epoch(self):
        # calculate the number of days since 1949 December 31 00:00 UT
        sgp4epoch = self.epoch - datetime.datetime(1949, 12, 31, 0, 0, 0) 
        # convert to number of days
        sgp4epoch = sgp4epoch.days
        return sgp4epoch
    
    def get_density(self, altitude, model = "exponential"):
        if model == "exponential":
            print("fetching density from exponential density model")
            return expo_simplified(altitude)
        #TODO: other density models here when ready (USSA 76 probably only one we need)
        else:
            print("fetching density from placeholder value")
        return 1e-12 #Placeholder value #in kg/m^3

    def sgp4_prop(self, jd):
        #jd is julian date of the epoch we want to propagate to
        #NOTE: from sgp4 docs
        # Note that ndot and nddot are ignored by the SGP4 propagator, 
        # so you can leave them 0.0 without any effect on the resulting satellite positions. 
        # But they do at least get saved to the satellite object, and written out if you write 
        # the parameters to a TLE or OMM file (see the “Export” section, above).
        sgp4_object = Satrec()
        altitude = (self.perigee_altitude + self.apogee_altitude)/2 #TODO: make this not allowed for non-cirular orbits
        density = self.get_density(altitude, model = "exponential")
        self.bstar = (self.C_d * self.characteristic_area * density)/2*self.mass #BStar = Cd * A * rho / 2m. Where Cd is the drag coefficient, A is the cross-sectional area of the satellite, rho is the density of the atmosphere, and m is the mass of the satellite.
        self.ndot = 0.0
        self.nddot = 0.0
        self.meananomaly = trueanom2meananom(self.tran, self.eccentricity)
        self.no_kozai = calculate_kozai_mean_motion(a = self.sma, mu = 398600.4418)
        self.satnum = 0 #Placeholder value. Believe this is NORAD ID? But doesnt really matter as not used in SGP4
        sgp4epoch = self.sgp4_epoch()
        sgp4_object.sgp4init(WGS72, 'i', self.satnum, sgp4epoch, self.bstar, self.ndot, self.nddot, self.eccentricity, self.argp, self.inc, self.meananomaly, self.no_kozai, self.raan)
        fr = 0.0
        error, position, velocity = sgp4_object.sgp4(jd, fr)
        if error != 0:
            print("Error in SGP4 propagation")
            return None
        else :
            self.cart_state = np.array(position + velocity)

def test_sgp4_prop():
    # quick test to make sure sgp4 propagator is working
    sat = Satellite(sma = 6378.137+500, perigee_altitude=250, apogee_altitude=250, ecc=0.0004, inc = np.deg2rad(85.9), argp = np.deg2rad(58), raan=np.deg2rad(88), tran=np.deg2rad(60), characteristic_area=20.0, mass = 100, epoch = "2022-04-16 00:00:00")
    print(sat.GUID)
    ephemeris = []
    dt_startday = datetime.datetime.strptime('2022-04-16 00:00:00', '%Y-%m-%d %H:%M:%S')
    dt_enday = datetime.datetime.strptime('2022-05-16 00:00:00', '%Y-%m-%d %H:%M:%S')
    start_day = utc_to_jd([dt_startday])
    end_day = utc_to_jd([dt_enday])
    step_size = 0.125/2 #in days
    date_range = np.arange(start_day[0], end_day[0], step_size)
    for date in date_range:
        sat.sgp4_prop(date)
        ephemeris.append(sat.cart_state)
        print("altitude:",np.linalg.norm(sat.cart_state[0:3])-6378.137)
    ephemeris = np.array(ephemeris)

    #plot ephemeris
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(ephemeris[:,0], ephemeris[:,1], ephemeris[:,2], c=date_range, s=1)
    #force aspect ratio to be 1:1:1
    ax.set_xlim(-7000, 7000)
    ax.set_ylim(-7000, 7000)
    ax.set_zlim(-7000, 7000)
    ax.legend()
    #add colorbar
    m = cm.ScalarMappable(cmap=cm.jet)
    m.set_array(date_range)
    fig.colorbar(m)
    plt.show()

if __name__ == "__main__":
    test_sgp4_prop()
    

   