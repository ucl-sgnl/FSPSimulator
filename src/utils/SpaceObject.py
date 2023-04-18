### class for an individual space object
class Satellite:
    def __init__(self, cospar_id, rso_name, rso_type, payload_operational_status, orbit_type, application, source, 
                 orbital_status_code, launch_site, decay_date, mass=None, maneuverable=False, spin_stabilized=False, 
                 orbital_period=None, inclination=None, apogee_altitude=None, perigee_altitude=None, radar_cross_section=None, 
                 characteristic_area=None, characteristic_length=None, propulsion_type=None, epoch=None, sma=None, inc=None, 
                 argp=None, raan=None, tran=None):
        
        self.CATID = None # This is a computed property based on the GUID and isn't included as a parameter
        self.GUID = None # This will be set by the system and isn't included as a parameter
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
        return f"{self.GUID:08d}"

    def get_cospar_id(self):
        return self.cospar_id
    
    def get_rso_name(self):
        return self.rso_name
    
    def get_rso_type(self):
        return self.rso_type
    
    def get_payload_operational_status(self):
        return self.payload_operational_status
    
    def get_orbit_type(self):
        return self.orbit_type
    
    def get_application(self):
        return self.application
    
    def get_source(self):
        return self.source
    
    def get_orbital_status_code(self):
        return self.orbital_status_code
    
    def get_launch_site(self):
        return self.launch_site
    
    def get_decay_date(self):
        return self.decay_date
    
    def get_mass(self):
        return self.mass
    
    def is_maneuverable(self):
        return self.maneuverable
    
    def is_spin_stabilized(self):
        return self.spin_stabilized
    
    def get_orbital_period(self):
        return self.orbital_period
    
    def get_inclination(self):
        return self.inclination
    
    def get_apogee_altitude(self):
        return self.apogee_altitude
    
    def get_perigee_altitude(self):
        return self.perigee_altitude
    
    def get_radar_cross_section(self):
        return self.radar_cross_section
    
    def get_characteristic_area(self):
        return self.characteristic_area
    
    def get_characteristic_length(self):
        return self.characteristic_length
    
    def get_propulsion_type(self):
        return self.propulsion_type
    
    def get_epoch(self):
        return self.epoch
    
    def get_sma(self):
        return self.sma
    
    def get_inc(self):
        return self.inc
    
    def get_argp(self):
        return self.argp
    
    def get_raan(self):
        return self.raan
    
    def get_tran(self):
        return self.tran