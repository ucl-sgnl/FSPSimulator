import math
import datetime
import warnings
from enum import Enum
import numpy as np
from astropy.time import Time
from fspsim.utils.Conversions import kep2car, true_to_mean_anomaly, orbital_period, get_day_of_year_and_fractional_day ,write_tle, utc_to_jd
from fspsim.utils.Propagators import kepler_prop

class OperationalStatus(Enum):
    POSITIVE = '+'
    NEGATIVE = '-'
    PARTIAL = 'P'
    B = 'B'
    S = 'S'
    X = 'X'
    D = 'D'
    UNKNOWN = '?'

class ObjectType(Enum):
    DEB = 'DEB'
    PAYLOAD = 'PAY'
    ROCKET_BODY = 'R/B'
    UNKNOWN = 'UNK'
    OTHER = '?'

def verify_value(value, impute_function):
    """
    Verifies that the value is a float and is not None. If it is None or not a float or 
    too small, then impute the value using the impute_function.

    :param value: The value to verify
    :type value: float or None or str
    :param impute_function: Function to impute value if necessary
    :type impute_function: function
    :return: The verified or imputed value
    :rtype: float
    """   
    try:
        value = float(value)
        if value >= 0.1:  # if value is acceptable
            return value
    except (TypeError, ValueError):
        pass
    return impute_function(value)  # pass the value to the impute_function

def verify_angle(value, name, random=False):
    """
    Verifies that the value is a float and is between 0 and 360. If it is None or not 
    a float or out of range, raise a ValueError.

    :param value: The value to verify
    :type value: float or None or str
    :param name: Name of the object
    :type name: str
    :param random: If `True`, return a random angle between 0 and 360 for invalid values, defaults to `False`
    :type random: bool, optional
    :return: The verified angle or a random angle if `random` is `True`
    :rtype: float
    :raises ValueError: If the value is not a float or out of range and `random` is `False`
    """
    try:
        value = float(value)
        if 0 <= value <= 360:  # if value is acceptable
            return value
        raise ValueError(f'Object {name} is negative or more than 360 degrees. {value} is not between 0 and 360')
    except (TypeError, ValueError):
        if random:
            return np.random.uniform(0, 360)
        raise ValueError(f'{value} is None or not a float. Please check the input data')

def verify_eccentricity(value):
    """
    Verifies that the value is a float and is between 0 and 1. If it is None or not a 
    float or out of range, raise a ValueError.

    :param value: The value to verify
    :type value: float or None or str
    :return: The verified value
    :rtype: float
    :raises ValueError: If the value is not a float or out of range
    """
    try:
        value = float(value)
        if 0 <= value <= 1:  # if value is acceptable
            return value
        raise ValueError('Object eccentricity is negative or more than 1. Please check the input data')
    except (TypeError, ValueError):
        raise ValueError('Object eccentricity is None or not a float. Please check the input data')

class SpaceObject:
    def __init__(self, rso_name=None, rso_type=None, payload_operational_status=None, application=None, source=None, 
                 launch_site=None, mass=None, maneuverable=False, spin_stabilized=False,
                 object_type = None, apogee=None, perigee=None, characteristic_area=None, characteristic_length=None, propulsion_type=None, epoch=None, sma=None, inc=None, 
                 argp=None, raan=None, tran=None, eccentricity=None, operator=None, launch_date=None,  decay_date=None, tle=None, station_keeping=None, orbit_source=None):
        
        # launch_date must be datetime.datetime.strptime('%Y-%m-%d') format between 1900-00-00 and 2999-12-31
        self.launch_date = datetime.datetime.strptime(launch_date, '%Y-%m-%d')
        self.decay_date = decay_date
        self.rso_name = str(rso_name)
        self.rso_type = str(rso_type) if rso_type is not None else None

        #Operational status
        self.payload_operational_status = str(payload_operational_status) # read in the value that is passed and cast to string
        possible_operational_status = ['+', '-', 'P', 'B', 'S', 'X', 'D', '?'] #https://celestrak.org/satcat/status.php -> this is the FSP operational status code
        if self.payload_operational_status not in possible_operational_status: 
            #raise a warning that still lets the code run and forces the operational status to be '?'
            # warnings.warn('WARNING: payload_operational_status should be one of the following: {}'.format(possible_operational_status))
            self.payload_operational_status = '?'

        #Object type
        self.object_type = str(object_type) # read in the value that is passed and cast to string
        possible_object_types = ['DEB', 'PAY', 'R/B', 'UNK', '?'] # i have added ? as a possible object type so that the code can still run if this is not specified
        if self.object_type not in possible_object_types:
            #raise a warning that still lets the code run and forces the object type to be '?'
            # print('WARNING: object_type should be one of the following: {}'.format(possible_object_types))
            self.object_type = '?'

        self.application = str(application) #TODO: get the code from FSP documentation for this and implement the checks

        self.operator = str(operator) # anything is fine here, if it is None or empty then we just set it to unknown
        if self.operator == None or self.operator == '':
            self.operator = 'Unknown'

        self.characteristic_length = verify_value(characteristic_length, self.impute_char_length)
        self.characteristic_area = verify_value(characteristic_area, self.impute_char_area)
        self.mass = verify_value(mass, self.impute_mass)

        self.source = source # http://celestrak.org/satcat/sources.php #TODO: do we really want source? This is just the country of origin, but I think maybe owner is more relevant nowadays
        self.launch_site = str(launch_site)
        self.decay_date = datetime.datetime.strptime(decay_date, '%Y-%m-%d') if (self.decay_date is not None and self.decay_date != '-') else datetime.datetime(2999, 1, 1) #set decay date to 01/01/2999 if is it None or '-' or empty
        #checl the difference between the two datetime objects is positive (i.e. decay date is after launch date)
        if (self.decay_date - self.launch_date).days < 0:
            raise ValueError('Object decay date cannot be earlier than launch date.')
        
        self.maneuverable = str(maneuverable)
        self.spin_stabilized = str(spin_stabilized)# I have left spin_stabilized and maneuverable as strings in case we later want to add more options than just True/False (for example different thruster types resulting in different kinds of maneuverability)
        
        # Apogee and Perigee
        #raise an error if these are None, 0, negative or more than 100000 
        self.apogee = float(apogee) #in Km
        if self.apogee is None or self.apogee == 0 or self.apogee < 0 or self.apogee > 200000:
            raise ValueError('Object apogee altitude is None, 0, negative or more than 200000km. Please check the input data')
        
        self.perigee = float(perigee) #in Km
        if self.perigee is None or self.perigee == 0 or self.perigee < 0 or self.perigee > 200000:
            raise ValueError('Object perigee altitude is None, 0, negative or more than 200000km. Please check the input data')
    
        self.propulsion_type = str(propulsion_type)
        self.epoch = epoch
        if self.epoch is None:
            self.epoch = datetime.datetime.utcnow()
        else:
            try:
                self.epoch = datetime.datetime.strptime(self.epoch, '%Y-%m-%d %H:%M:%S') #in UTC
            except:
                self.epoch = datetime.datetime.strptime(self.epoch, '%Y-%m-%dT%H:%M:%S.%f') # this is space-track

        self.day_of_year = get_day_of_year_and_fractional_day(self.epoch)
        self.station_keeping = station_keeping

        if self.station_keeping == True: # if station_keeping is just set to True, then the object will station keep from launch to decay
            self.station_keeping = [self.launch_date, self.decay_date] 
        elif self.station_keeping == False or self.station_keeping == None:
            self.station_keeping = None
        else:
            #if only one is specified then we assume that station keeping is from launch to that date
            if len(self.station_keeping) == 1:
                self.station_keeping = [datetime.datetime.strptime(self.station_keeping[0], '%Y-%m-%d %H:%M:%S')]
            #if two dates are specified, we assume these are the start and end dates of station keeping
            elif len(self.station_keeping) == 2:
                self.station_keeping = [datetime.datetime.strptime(self.station_keeping[0], '%Y-%m-%d %H:%M:%S'), datetime.datetime.strptime(self.station_keeping[1], '%Y-%m-%d')]
            else:
                raise ValueError('station_keeping must be either None, True, False or a list of length 1 or 2')
        
        self.ephemeris = None # This is where the ephemeris of the propagations will be stored

        # Semi-major axis
        self.sma = (self.apogee + self.perigee)/2 + 6378.137 #in km #TODO: this is not correct for non-circular orbits.
        #raise an error if sma is None, 0, negative or more than 100000
        if self.sma is None or self.sma == 0 or self.sma < 0 or self.sma > 200000:
            raise ValueError('Object sma is None, 0, negative or more than 100000km. Please check the input data')
        
        # Orbital period
        self.orbital_period = orbital_period(self.sma) #in minutes
        #raise an error if orbital_period is None, 0, or negative
        if self.orbital_period is None or self.orbital_period == 0 or self.orbital_period < 0:
            raise ValueError('Object orbital_period is None, 0, or negative. Please check the input data')
        
        self.inc = verify_angle(inc, 'inc')
        self.argp = verify_angle(argp, 'argp')
        self.raan = verify_angle(raan, 'raan')
        self.tran = verify_angle(tran, 'tran', random=True) # random=True means that if tran is None, then a random value between 0 and 360 is imputed
        self.eccentricity = verify_eccentricity(eccentricity)
        self.meananomaly = true_to_mean_anomaly(self.tran, self.eccentricity)

        # This is to store the result of the conversion of the Keplerian elements to cartesian ECI coordinates
        self.cart_state = self.generate_cart() #cartesian state vector [x,y,z,u,v,w] to be computed using generate_cart (from keplerian elements)

        self.C_d = 2.2 #Drag coefficient

        mean_motion = 2*math.pi/self.orbital_period # in radians per minute
        mean_motion_revs_per_day = mean_motion * 24 * 60 / (2*math.pi) # in revolutions per day

        # if tle is not a string of the format: "{69 characters}\n{69 characters}""  i.e. 2 lines of 69 characters, then raise a warning and set tle to None
        self.tle=write_tle(catalog_number=00000, classification="U", launch_year=2020, launch_number= 67, launch_piece="A",
              epoch_year=20, epoch_day=264.51782528, first_derivative='-.00002182', second_derivative='00000-0', drag_term='-11606-4',
              ephemeris_type=0, element_set_number=292,inclination= self.inc, raan=self.raan, eccentricity=self.eccentricity, arg_perigee=self.argp,
              mean_anomaly=self.meananomaly, mean_motion=mean_motion_revs_per_day,revolution_number=56353)
        
        if len(self.tle) != 2*69+1 or self.tle[69] != '\n':
            warnings.warn('WARNING: tle must be a string of the format: "{69 characters}\\n{69 characters}" i.e. 2 lines of 69 characters. Setting tle to None')
            # Example TLE that works: "1 53544U 22101T   23122.20221856  .00001510  00000-0  11293-3 0  9999\n2 53544  53.2176  64.0292 0001100  79.8127 280.2989 15.08842383 38928"

    def impute_char_length(self, char_length):
        """
        Imputes a characteristic length based on the object type.

        :param char_length: Characteristic length to verify or impute.
        :type char_length: float or None
        :return: Imputed or verified characteristic length.
        :rtype: float
        :raises ValueError: If the imputed or given char_length is None or 0.
        """
        if char_length is None or char_length == 0:
            if self.object_type == 'PAYLOAD':
                char_length = 1.3
            elif self.object_type == "ROCKET BODY":
                char_length = 8
            else:
                char_length = 1.3
        # raise an error if this is None or 0
        if char_length is None or char_length == 0:
            raise ValueError('Object characteristic length is None or 0. Please check the input data')
        
        return char_length

    def impute_char_area(self, char_area):
        """
        Imputes a characteristic area based on the object type.

        :param char_area: Characteristic area to verify or impute.
        :type char_area: float or None
        :return: Imputed or verified characteristic area.
        :rtype: float
        :raises ValueError: If the imputed or given char_area is None or 0.
        """
        if char_area is None or char_area == 0:
            if self.object_type == 'PAYLOAD':
                char_area = 1.3 * 3
            elif self.object_type == "ROCKET BODY":
                char_area = 8 * 3
            else:
                char_area = 1.3 * 3
        # raise an error if this is None or 0
        if char_area is None or char_area == 0:
            raise ValueError('Object characteristic area is None or 0. Please check the input data')

        return char_area
        
    def impute_mass(self, mass):
        """
        Imputes a mass based on the object type.

        :param mass: Mass to verify or impute.
        :type mass: float or None or str
        :return: Imputed mass value.
        :rtype: float
        """
        # From Alfano, Oltrogge and Sheppard, 2020
        if mass is None or mass == 0: # if mass is None or 0, then impute it based on the object length
            imputed_mass = 0.1646156590 * self.characteristic_length ** (-1.0554951607) # TODO:this basically gives ~0.12 Kgs for a 1.3m length object. Means cubesats are very light. Maybe we should change this?
        else:
            imputed_mass = float(mass) # cast to float in case it is a string
        return imputed_mass
    
    def generate_cart(self):
        """
        Generates a cartesian state vector from keplerian elements.

        :return: Cartesian state.
        :rtype: np.array
        """
        x,y,z,u,v,w = kep2car(a = self.sma, e=self.eccentricity, i = math.radians(self.inc), w = math.radians(self.argp), W=math.radians(self.raan), V=math.radians(self.tran), epoch=Time(self.epoch, format='datetime'))
        self.cart_state = np.array([[x, y, z], [u, v, w]])
        return self.cart_state

    def prop_catobject(self, jd_start, jd_stop, step_size, output_freq):
        """
        Propagates a celestial object based on initial conditions, propagator type, and station keeping preferences.

        :param jd_start: Julian start date for the simulation.
        :type jd_start: float
        :param jd_stop: Julian stop date for simulation.
        :type jd_stop: float
        :param step_size: Step size for propagation.
        :type step_size: float
        :param output_freq: Frequency at which to output the ephemeris (in seconds).
        :type output_freq: float
        :param integrator_type: Numerical integrator to use.
        :type integrator_type: str
        :param use_sgp4_propagation: Propagate using SGP4 for 100-minute segments.
        :type use_sgp4_propagation: bool
        :return: None. Updates the `ephemeris` attribute of the object.
        """
        
        #check launch date - this is for objects that are launched after the start of the simulation
        if utc_to_jd(self.launch_date) > jd_start:
            jd_start = utc_to_jd(self.launch_date)

        # If output_freq is not specified, set it to equal the step_size
        if output_freq is None:
            output_freq = step_size
        else:
            output_freq_steps = max(1, round(output_freq / step_size))  # Calculate the output frequency in steps

        if self.station_keeping == True:
        # Object will station keep from launch to decay
            print("kepler for: ", self.rso_name)
            self.ephemeris = kepler_prop(jd_start, jd_stop, step_size, a=self.sma, e=self.eccentricity, i=self.inc, w=self.argp, W=self.raan, V=self.tran)
            self.ephemeris = self.ephemeris[::output_freq_steps]

        # Check if station keeping is specified as a list. If it is propagate using Kepler from the start date to the start of the station keeping date
        elif isinstance(self.station_keeping, list):
            # Propagate using Kepler from the start date to the end of the station keeping date
            # This assumes that if station keeping is going to occurr, it always occurs from the beggining of the orbit
            station_keeping_start_jd = utc_to_jd(self.station_keeping[0])
            station_keeping_end_jd = utc_to_jd(self.station_keeping[1])
            ephemeris_station_keep = kepler_prop(station_keeping_start_jd, station_keeping_end_jd, step_size, a=self.sma, e=self.eccentricity, i=self.inc, w=self.argp, W=self.raan, V=self.tran, drag_decay=False)

            # Propagate using dragdecay for the rest of the orbit
            dragdecay_ephemeris = kepler_prop(jd_start=station_keeping_end_jd, jd_stop=jd_stop, step_size=step_size, a=self.sma, e=self.eccentricity, i=self.inc, w=self.argp, W=self.raan, V=self.tran, area=self.characteristic_area, mass=self.mass, cd=self.C_d, drag_decay=True)

            ephemeris_station_keep = np.array(ephemeris_station_keep)
            dragdecay_ephemeris = np.array(dragdecay_ephemeris)

            # Combine arrays
            self.ephemeris = np.concatenate((ephemeris_station_keep, dragdecay_ephemeris), axis=0)

        else: # if no station keeping involved at all, then just propagate using dragdecay the whole way through
            # Object will not station keep, propagate using only decay
            self.ephemeris = kepler_prop(jd_start=jd_start, jd_stop=jd_stop, step_size=step_size, a=self.sma, e=self.eccentricity, i=self.inc, w=self.argp, W=self.raan, V=self.tran, area=self.characteristic_area, mass=self.mass, cd=self.C_d, drag_decay=True)
            self.ephemeris = self.ephemeris[::output_freq_steps]

if __name__ == "__main__":
    pass