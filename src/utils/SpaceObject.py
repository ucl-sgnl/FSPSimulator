### class for an individual space object
import uuid
class Satellite:
    def __init__(self, cospar_id=None, rso_name=None, rso_type=None, payload_operational_status=None, orbit_type=None, application=None, source=None, 
                 orbital_status_code=None, launch_site=None, decay_date=None, mass=None, maneuverable=False, spin_stabilized=False, 
                 orbital_period=None, inclination=None, apogee_altitude=None, perigee_altitude=None, radar_cross_section=None, 
                 characteristic_area=None, characteristic_length=None, propulsion_type=None, epoch=None, sma=None, inc=None, 
                 argp=None, raan=None, tran=None):
        
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
        self.mass = mass
        self.maneuverable = maneuverable
        self.spin_stabilized = spin_stabilized
        self.orbital_period = orbital_period
        self.inclination = inclination
        self.apogee_altitude = apogee_altitude
        self.perigee_altitude = perigee_altitude
        self.radar_cross_section = radar_cross_section
        self.characteristic_area = characteristic_area
        self.characteristic_length = characteristic_length
        self.propulsion_type = propulsion_type
        self.epoch = epoch
        self.sma = sma
        self.inc = inc
        self.argp = argp
        self.rann = raan
        self.train = tran
    
        
    # This is a private method because its name starts with an underscore
    def _compute_catid(self):
        return uuid.uuid4()