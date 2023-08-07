import math
import datetime
import warnings
from enum import Enum
import numpy as np
from astropy.time import Time
from utils.Conversions import kep2car, true_to_mean_anomaly, orbital_period, get_day_of_year_and_fractional_day,prep_ephemeris_for_tle_fitting, fit_TLE_to_ephemeris
from utils.Propagators import numerical_prop, kepler_prop, sgp4_prop_TLE

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
    """Verifies that the value is a float and is not None. If it is None or not a float or too small, then impute the value using the impute_function."""
    try:
        value = float(value)
        if value >= 0.1:  # if value is acceptable
            return value
    except (TypeError, ValueError):
        pass
    return impute_function(value)  # pass the value to the impute_function

def verify_angle(value, name, random=False):
    """Verifies that the value is a float and is between 0 and 360. If it is None or not a float or out of range, raise a ValueError."""
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
    """Verifies that the value is a float and is between 0 and 1. If it is None or not a float or out of range, raise a ValueError."""
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
                 argp=None, raan=None, tran=None, eccentricity=None, operator=None, launch_date=None,  decay_date=None, tle=None, station_keeping=None):
        
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
        
        # if tle is not a string of the format: "{69 characters}\n{69 characters}""  i.e. 2 lines of 69 characters, then raise a warning and set tle to None
        self.tle=tle
        if tle is not None:
            if len(tle) != 2*69+1 or tle[69] != '\n':
                warnings.warn('WARNING: tle must be a string of the format: "{69 characters}\\n{69 characters}" i.e. 2 lines of 69 characters. Setting tle to None')
                # Example TLE that works: "1 53544U 22101T   23122.20221856  .00001510  00000-0  11293-3 0  9999\n2 53544  53.2176  64.0292 0001100  79.8127 280.2989 15.08842383 38928"
                self.tle = None
            else:
                self.tle = tle
        elif tle is None:
            self.tle = None

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
    
    def impute_char_length(self, char_length):
        """
        This function will impute a characteristic length based on the object type
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
        This function will impute a characteristic area based on the object type
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
        This function will impute a mass based on the object type
        """
        # From Alfano, Oltrogge and Sheppard, 2020
        if mass is None or mass == 0: # if mass is None or 0, then impute it based on the object length
            imputed_mass = 0.1646156590 * self.characteristic_length ** (-1.0554951607) # TODO:this basically gives ~0.12 Kgs for a 1.3m length object. Means cubesats are very light. Maybe we should change this?
        else:
            imputed_mass = float(mass) # cast to float in case it is a string
        return imputed_mass
    
    def generate_cart(self):
        # generate cartesian state vector from keplerian elements
        # Keplerian elements are in radians
        
        #This is the format of the epoch: datetime.datetime.strptime(self.epoch, '%Y-%m-%d %H:%M:%S') #in UTC
        # COnvert it to astropy time
        x,y,z,u,v,w = kep2car(a = self.sma, e=self.eccentricity, i = math.radians(self.inc), w = math.radians(self.argp), W=math.radians(self.raan), V=math.radians(self.tran), epoch=Time(self.epoch, format='datetime'))
        self.cart_state = np.array([[x, y, z], [u, v, w]])
        return self.cart_state

def prop_catobject(self, jd_start, jd_stop, step_size, output_freq, integrator_type, force_model, long_term_sgp4):
    """
    Function to propagate a celestial object based on initial conditions, propagator type, and station keeping preferences.

    Parameters:
    jd_start (float): Julian start date for propagation
    jd_stop (float): Julian stop date for propagation
    step_size (float): Step size for propagation
    output_freq (float): Frequency at which to output the ephemeris (in seconds)
    integrator_type (str): String indicating which numerical integrator to use, default is "RK45"
    use_sgp4_propagation (bool): Propagate using SGP4 for 100-minute segments. Default is False.

    Returns:
    None: The function does not return anything but updates the `ephemeris` attribute of the object.
    """
    valid_integrator_types = ["RK45"]
    if integrator_type not in valid_integrator_types:
        raise ValueError(f"Invalid integrator. Must be one of the following: {valid_integrator_types}")

    tot_time = (jd_stop - jd_start) * 24 * 60 * 60  # calculate total time in seconds for the propagation
    # If output_freq is not specified, set it to equal the step_size
    if output_freq is None:
        output_freq = step_size
    else:
        output_freq_steps = max(1, round(output_freq / step_size))  # Calculate the output frequency in steps

    combined_ephemeris = []

    if long_term_sgp4:
        numerical_segment_minutes = 30 # 30 mins of higher fidelity numerical propagation to fit the TLE to
        numerical_segment_seconds = numerical_segment_minutes * 60

        # Check if station keeping is specified
        if isinstance(self.station_keeping, list):
            # Propagate using Kepler from the start date to the end of the station keeping date
            ephemeris_station_keep = kepler_prop(jd_start, self.station_keeping[1], step_size, a=self.sma, e=self.eccentricity, i=self.inc, w=self.argp, W=self.raan, V=self.tran)
            combined_ephemeris += ephemeris_station_keep
            current_jd = self.station_keeping[1]
        else:
            current_jd = jd_start

        # Propagate numerically for the segment time
        next_jd = min(current_jd + numerical_segment_seconds / 86400, jd_stop)
        ephemeris_numerical = numerical_prop(tot_time=(next_jd - current_jd) * 86400, pos=self.cart_state[0], vel=self.cart_state[1], C_d=self.C_d, area=self.characteristic_area, mass=self.mass, JD_time_start=current_jd, integrator_type=integrator_type, force_model=force_model)
        positions_eci, velocities_eci, mjds =  prep_ephemeris_for_tle_fitting(ephemeris_numerical)
        # Fit TLE from numerical ephemeris
        TLE = fit_TLE_to_ephemeris(positions_eci, velocities_eci, mjds)
        line1 = TLE.getLine1()
        line2 = TLE.getLine2()
        tle_string = line1 + '\n' + line2
        # Propagate using SGP4 for the rest of the orbit
        print(f"Propagating {self.rso_name} from {next_jd} to {jd_stop} using SGP4")
        ephemeris_sgp4 = sgp4_prop_TLE(tle_string, next_jd, jd_stop, step_size)

        # Concatenate results
        combined_ephemeris += ephemeris_numerical + ephemeris_sgp4
        print(f"Combined ephemeris length: {len(combined_ephemeris)}")

        # Update ephemeris attribute
        self.ephemeris = np.array(combined_ephemeris)[::output_freq_steps]

    elif self.station_keeping == True:
        # Object will station keep from launch to decay
        self.ephemeris = kepler_prop(jd_start, jd_stop, step_size, a=self.sma, e=self.eccentricity, i=self.inc, w=self.argp, W=self.raan, V=self.tran)
        self.ephemeris = self.ephemeris[::output_freq_steps]

    elif not self.station_keeping:
        # Object will not station keep, propagate using the numerical integrator
        self.ephemeris = numerical_prop(tot_time=tot_time, pos=self.cart_state[0], vel=self.cart_state[1], C_d=self.C_d, area=self.characteristic_area, mass=self.mass, JD_time_start=jd_start, integrator_type=integrator_type, force_model=force_model)
        self.ephemeris = self.ephemeris[::output_freq_steps]

if __name__ == "__main__":
    pass