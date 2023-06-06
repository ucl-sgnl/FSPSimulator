### class for an individual space object
import uuid
import math 
import datetime
import numpy as np
import sgp4
import matplotlib.pyplot as plt
from sgp4.api import Satrec, WGS72
from utils.coords import kep2car, true_to_mean_anomaly, calculate_kozai_mean_motion, expo_simplified, utc_to_jd, tle_parse, tle_convert, sgp4_prop_TLE, write_tle, orbital_period, get_day_of_year_and_fractional_day, TLE_time, jd_to_utc, kepler_prop
import matplotlib.cm as cm

def verify_value(value, impute_function):
    # function to verify that the value is a float and is not None
    # if it is None or not a float or too small, then impute the value using the impute_function
    try:
        value = float(value)
        if value >= 0.1:  # if value is acceptable
            return value
    except (TypeError, ValueError):
        pass
    return impute_function(value)  # pass the value to the impute_function

def verify_angle(value, name, random=False):
    # function to verify that the value is a float and is between 0 and 360
    # if it is None or not a float or out of range, raise a ValueError
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
    # function to verify that the value is a float and is between 0 and 1
    # if it is None or not a float or out of range, raise a ValueError
    try:
        value = float(value)
        if 0 <= value <= 1:  # if value is acceptable
            return value
        raise ValueError('Object eccentricity is negative or more than 1. Please check the input data')
    except (TypeError, ValueError):
        raise ValueError('Object eccentricity is None or not a float. Please check the input data')


class SpaceObject:
    def __init__(self, rso_name=None, rso_type=None, payload_operational_status=None, orbit_type=None, application=None, source=None, 
                 orbital_status_code=None, launch_site=None, mass=None, maneuverable=False, spin_stabilized=False, bstar = None, 
                 object_type = None, apogee=None, perigee=None, radar_cross_section=None, 
                 characteristic_area=None, characteristic_length=None, propulsion_type=None, epoch=None, sma=None, inc=None, 
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
            # print('WARNING: payload_operational_status should be one of the following: {}'.format(possible_operational_status))
            self.payload_operational_status = '?'
            # raise ValueError('payload_operational_status must be one of the following: {}'.format(possible_operational_status)) # raise an error that stops the code
        
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
        # self.orbital_status_code = str(orbital_status_code) #TODO: I also would argue this is not that useful. We have payload_operational_status (especially if focus is just Earth orbit this is totally redundant i believe)
        self.launch_site = str(launch_site)
        #set decay date to 01/01/2999 if is it None or '-' or empty
        self.decay_date = datetime.datetime.strptime(decay_date, '%Y-%m-%d') if (self.decay_date is not None and self.decay_date != '-') else datetime.datetime(2999, 1, 1)
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
    
        self.radar_cross_section = float(radar_cross_section) if radar_cross_section is not None else None #in meters^2
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
                print('WARNING: tle must be a string of the format: "{69 characters}\\n{69 characters}" i.e. 2 lines of 69 characters. Setting tle to None')
                # Example TLE that works: "1 53544U 22101T   23122.20221856  .00001510  00000-0  11293-3 0  9999\n2 53544  53.2176  64.0292 0001100  79.8127 280.2989 15.08842383 38928"
                self.tle = None
            else:
                self.tle = tle

        self.ndot = 0.00000140 #First derivative of mean motion- hardcoded for time being. Minimally impacts propagation
        self.nddot = 0.0 # If none this is set to 0.0. As above, minimally impacts propagation
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
            
        self.altitude = self.perigee #TODO: this needs to get deleted imo (discuss first)
        
        self.get_atmospheric_density(model = "exponential") #exponential here is ~USSA 76

        self.C_d = 2.2 #Drag coefficient
        
        if bstar is None:
            self.bstar = 0.00113295 # impute the median BStar of the SpaceTrack BStars that are below 1000km
            
            #-------- Compute BStar method -------# 
            # Not using as decay is crazy with this method   
            # -(self.C_d * self.characteristic_area * self.atmos_density)/2*self.mass #BStar = Cd * A * rho / 2m. Where Cd is the drag coefficient, A is the cross-sectional area of the satellite, rho is the density of the atmosphere, and m is the mass of the satellite.
            # print("Setting bstar to {}".format(self.bstar))
            # print("Area: {}".format(self.characteristic_area))
            # print("Mass: {}".format(self.mass))
            # print("Density: {}".format(self.atmos_density))
        else:
            self.bstar = bstar
        self.no_kozai = calculate_kozai_mean_motion(a = self.sma, mu = 398600.4418)

        self.sgp4epoch = self.sgp4_epoch() #SGP4 epoch is the number of days since 1949 December 31 00:00 UT

        # This is to store the result of the conversion of the Keplerian elements to cartesian ECI coordinates
        self.cart_state = None #cartesian state vector [x,y,z,u,v,w] to be computed using generate_cart (from keplerian elements)


    def _validate_types(self):
        # function to validate the types and values of the parameters
        possible_operational_status = ['+', '-', 'P', 'B', 'S', 'X', 'D', '?'] #https://celestrak.org/satcat/status.php
        if self.payload_operational_status not in possible_operational_status:
            raise ValueError('payload_operational_status must be one of the following: {}'.format(possible_operational_status))
        
        possible_spin_stabilized = ['y', 'n']
        if self.spin_stabilized not in possible_spin_stabilized:
            raise ValueError('spin_stabilized must be one of the following: {}'.format(possible_spin_stabilized))
        
        possible_maneuverable = ['y', 'n']
        if self.maneuverable not in possible_maneuverable:
            raise ValueError('maneuverable must be one of the following: {}'.format(possible_maneuverable))
        
        possible_launch_sites = ['AFETR','AFWTR','CAS','DLS','ERAS','FRGUI','HGSTR','JSC','KODAK','KSCUT','KWAJ','KYMSC','NSC','PLMSC','RLLB','SEAL',
                                 'SEMLS','SMTS','SNMLP','SRILR','SUBL','SVOBO','TAISC','TANSC','TYMSC','UNK','VOSTO','WLPIS','WOMRA','WRAS','WSC','XICLF','YAVNE','YUN']
        if self.launch_site not in possible_launch_sites:
            raise ValueError('launch_site must be one of the following: {}'.format(possible_launch_sites))
        
        possible_object_types = ['DEB', 'PAY', 'R/B', 'UNK']
        if self.object_type not in possible_object_types:
            raise ValueError('object_type must be one of the following: {}'.format(possible_object_types))
        
        #check that inc, argp, raan and tran are all between 0 and 2pi
        if self.inc < 0 or self.inc > 2*math.pi:
            raise ValueError('inc must be in Radians. Current value not between 0 and 2pi')
        if self.argp < 0 or self.argp > 2*math.pi:
            raise ValueError('argp must be in Radians. Current value not between 0 and 2pi')
        if self.raan < -2*math.pi or self.raan > 2*math.pi:
            raise ValueError('raan must be in Radians. Current value not between 0 and 2pi')
        if self.tran < -2*math.pi or self.tran > 2*math.pi:
            raise ValueError('tran must be in Radians. Current value not between 0 and 2pi')
        if self.eccentricity<0 or self.eccentricity>1:
            raise ValueError('eccentricity must be between 0 and 1')
    
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
            imputed_mass = 0.1646156590 * self.characteristic_length ** (-1.0554951607) # TODO:this basically gives ~0.12 Kgs for a 1.3m length object. Hm...
        else:
            imputed_mass = float(mass) # cast to float in case it is a string
        return imputed_mass

    def decode(self, code):
        # This function decodes all the codes that are associated with an object. Basically a legend associated to each instance of the class
        #"code" is the attribute you want to decode
        code_mapping = {
            'object_type': {
                'DEB': 'Debris',
                'PAY': 'Payload',
                'R/B': 'Rocket Body',
                'UNK': 'Unknown'
            },
            'launch_site': {
                'AFETR': 'Air Force Eastern Test Range, Florida, USA',
                'AFWTR': 'Air Force Western Test Range, California, USA',
                'CAS': 'Canaries Airspace',
                'DLS': 'Dombarovskiy Launch Site, Russia',
                'ERAS': 'Eastern Range Airspace',
                'FRGUI': "Europe's Spaceport, Kourou, French Guiana",
                'HGSTR': 'Hammaguira Space Track Range, Algeria',
                'JSC': 'Jiuquan Space Center, PRC',
                'KODAK': 'Kodiak Launch Complex, Alaska, USA',
                'KSCUT': 'Uchinoura Space Center (Formerly Kagoshima Space Center—University of Tokyo, Japan)',
                'KWAJ': 'US Army Kwajalein Atoll (USAKA)',
                'KYMSC': 'Kapustin Yar Missile and Space Complex, Russia',
                'NSC': 'Naro Space Complex, Republic of Korea',
                'PLMSC': 'Plesetsk Missile and Space Complex, Russia',
                'RLLB': 'Rocket Lab Launch Base, Mahia Peninsula, New Zealand',
                'SEAL': 'Sea Launch Platform (mobile)',
                'SEMLS': 'Semnan Satellite Launch Site, Iran',
                'SMTS': 'Shahrud Missile Test Site, Iran',
                'SNMLP': 'San Marco Launch Platform, Indian Ocean (Kenya)',
                'SRILR': 'Satish Dhawan Space Centre, India (Formerly Sriharikota Launching Range)',
                'SUBL': 'Submarine Launch Platform (mobile)',
                'SVOBO': 'Svobodnyy Launch Complex, Russia',
                'TAISC': 'Taiyuan Space Center, PRC',
                'TANSC': 'Tanegashima Space Center, Japan',
                'TYMSC': 'Tyuratam Missile and Space Center, Kazakhstan (Also known as Baikonur Cosmodrome)',
                'UNK': 'Unknown',
                'VOSTO': 'Vostochny Cosmodrome, Russia',
                'WLPIS': 'Wallops Island, Virginia, USA',
                'WOMRA': 'Woomera, Australia',
                'WRAS': 'Western Range Airspace',
                'WSC': 'Wenchang Satellite Launch Site, PRC',
                'XICLF': 'Xichang Launch Facility, PRC',
                'YAVNE': 'Yavne Launch Facility, Israel',
                'YUN': "Yunsong Launch Site (Sohae Satellite Launching Station), Democratic People's Republic of Korea"
            },
            'operational_status': {
                '+': 'Operational',
                '-': 'Nonoperational',
                'P': 'Partially Operational',
                'B': 'Backup/Standby',
                'S': 'Spare',
                'X': 'Extended Mission',
                'D': 'Decayed',
                '?': 'Unknown'
            },
            # Add other mappings here as/when needed
        }

        for key, value in code_mapping.items():
            if code in value:
                return value[code]
            else:
                print('Code not found in mapping. Returning None')
        return None

    def _compute_catid(self):
        # generate a unique ID for the object. This is the internal ID used in this catalog.
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
    
    def get_atmospheric_density(self, model = "exponential"):
        if model == "exponential":
            # self.atmos_density = 1e-12
            self.atmos_density = expo_simplified(self.altitude) #fix units to kg/m^3 
        #TODO: other density models here when ready (USSA 76 probably only one we need)
        else:
            print("WARNING: Atmospheric density model not recognized. Setting density to 1e-12")
            self.atmos_density = 1e-12 #Placeholder value #in kg/m^3
 
    def build_TLE(self):
        TLE = write_tle(
            int(12345),  # catalog_number -> this is a placeholder value and does not affect the propagation TODO: decide whether we need to impute this correctly
            "U",  # classification -> also a placeholder value TODO: impute correctly.
            int(self.epoch.year % 100),  # launch_year -> also a placeholder value
            int(29),  # launch_number -> also a placeholder value TODO: impute correctly.
            'AN',  # launch_piece -> also a placeholder value TODO: impute correctly.
            #last two digits of launch year
            self.epoch.year % 100,
            self.day_of_year, # epoch_day
            self.ndot,  # first_derivative
            self.nddot,  # second_derivative #Typically this is set to 0
            self.bstar,  # drag_term
            0,  # ephemeris_type: Always 0.
            0,  # element_set_number. No impact on propagation this is just Elset Number.
            self.inc,
            self.raan,
            self.eccentricity,
            self.argp,
            self.meananomaly,
            self.no_kozai,
            0,  # revolution_number
        )
        # assign the new TLE to the object
        self.tle = TLE

    def prop_catobjects(self, jd_start, jd_stop, step_size):
        
        if isinstance(self.station_keeping, list): # if a list is passed to the station keeping attribute, then we assume that the first element is station keeping start date and the second element is the end date)
            # propagate using keplerian from the start date to the end of the station keeping date 
            # TODO: this assumes that the satellite always station keeps from the moment it enters orbit. We should be using self.station_keeping[0]
            # then propagate using SGP4 from the second element in the list to jd_stop
            combined_ephemeris = []
            ephemeris_station_keep = kepler_prop(jd_start, self.station_keeping[1], step_size, a=self.sma, e=self.eccentricity, i=self.inc, w=self.argp, W=self.raan, V=self.tran)
            combined_ephemeris.append(ephemeris_station_keep)
            if self.tle is None:
                self.build_TLE()
            ephemeris_sgp4 = sgp4_prop_TLE(self.station_keeping[1], jd_stop, step_size, self.tle)
            combined_ephemeris.append(ephemeris_sgp4)
            self.ephemeris = np.concatenate(combined_ephemeris) # concatenate the two ephemeris arrays into one and append it to the ephemeris attribute
        elif self.station_keeping == True:
            self.ephemeris = kepler_prop(jd_start, jd_stop, step_size, a=self.sma, e=self.eccentricity, i=self.inc, w=self.argp, W=self.raan, V=self.tran)
        elif self.station_keeping == False or self.station_keeping is None: 
            if self.tle is None:
                self.build_TLE()
            self.ephemeris = sgp4_prop_TLE(self.tle, jd_start, jd_stop, step_size)

def test_sgp4_drag():
    """
    Test the implementation of the SGP4 propagator. In particular looking at the decay rates and the effect of the way BSTAR is calculated in the construction of TLEs
    """

    #navstar81, pulled on 27Apr2023
    test_tle1 = "1 48859U 21054A   23116.83449170 -.00000109  00000-0  00000-0 0  9996\n2 48859  55.3054  18.4561 0008790 213.9679 183.6522  2.00556923 13748"
    #OneWeb20, pulled on 27Apr2023
    test_tle2 = "1 45133U 20008C   23116.69886660  .00000678  00000-0  19052-2 0  9998\n2 45133  87.8784 264.1991 0001687  89.2910 270.8411 13.10377378158066"
    #Starlink70, pulled on 27Apr2023
    test_tle3 = "1 53544U 22101T   23122.20221856  .00001510  00000-0  11293-3 0  9999\n2 53544  53.2176  64.0292 0001100  79.8127 280.2989 15.08842383 38928"

    #dictionary of TLEs
    test_tles = {'navstar81': test_tle1, 'OneWeb20': test_tle2, 'Starlink70': test_tle3}

    prop_start = [datetime.datetime.strptime('2023-05-02 12:45:00', '%Y-%m-%d %H:%M:%S')]
    prop_end = [datetime.datetime.strptime('2028-04-27 01:48:00', '%Y-%m-%d %H:%M:%S')]
    start_jd = utc_to_jd(prop_start)
    end_jd = utc_to_jd(prop_end)
    t_step = 60*60*24 #1 day in seconds

    #test propagation of TLEs in test_tles
    tseries_altitude = [] #for each satellite, this will contain two lists (one for the real and one for the fabricated TLE) of altitude values
    t_series_pos = [] #for each satellite, this will be a list of position differences between the real and fabricated TLEs

    for sat in test_tles:
        print("Propagating TLE for ", sat)
        sat_altitude = []#list of real and fabricated altitudes for this satellite
        sat_pos = []#list of real and fabricated positions for this satellite
 
        ###### SECTION 1: Make TLEs and propagate them ######
        #Get the Keplerian elements from the TLE
        tle_dict = tle_parse(test_tles[sat])
        tle_kepels = tle_convert(tle_dict)
        #Get the epoch from the TLE
        tle_time = TLE_time(test_tles[sat])
        tle_epoch = jd_to_utc(tle_time)
        # make a SpaceObject from the TLE
        tle_epoch_str = str(tle_epoch)
        epoch = tle_epoch_str.replace(' ', 'T')
        test_sat = SpaceObject(sma = tle_kepels['a'], perigee=tle_kepels['a']-6378.137, apogee=tle_kepels['a']-6378.137, eccentricity=tle_kepels['e'], inc = tle_kepels['i'], argp = tle_kepels['arg_p'], raan=tle_kepels['RAAN'], tran=tle_kepels['true_anomaly'], characteristic_area=0.011, mass = 250, epoch = epoch, launch_date='2023-05-02')
        test_sat.prop_catobjects(start_jd[0], end_jd[0], t_step)
        test_sat_ephem = test_sat.ephemeris
        print("made up TLE:", test_sat.tle)
        print("real TLE:", test_tles[sat])
        # test_sat_ephem is list of a tuples of the form [(time, position, velocity), (time, position, velocity), ...]

        test_pos = [x[1] for x in test_sat_ephem]
        test_pos = np.array(test_pos)
        test_times = [x[0] for x in test_sat_ephem]
        test_altitudes = [np.linalg.norm(x)-6378.137 for x in test_pos]
        sat_altitude.append(test_altitudes)
        sat_pos.append(test_pos)

        ###### SECTION 2: Use existing TLEs and propagate them ######
        valid_tle = test_tles[sat]
        valid_tle_ephem = sgp4_prop_TLE(valid_tle, jd_start=start_jd[0], jd_end=end_jd[0], dt=t_step)
        valid_tle_pos = [x[1] for x in valid_tle_ephem]
        valid_tle_pos = np.array(valid_tle_pos)
        valid_tle_times = [x[0] for x in valid_tle_ephem]
        valid_tle_altitudes = [np.linalg.norm(x)-6378.137 for x in valid_tle_pos]
        sat_altitude.append(valid_tle_altitudes)
        sat_pos.append(valid_tle_pos)

        tseries_altitude.append(sat_altitude)
        t_series_pos.append(sat_pos)
        
#plot the altitude of each satellite over time
    for i in range (len(tseries_altitude)): 
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(tseries_altitude[i][0], label = list(test_tles.keys())[i] + 'fabricated')
        ax.plot(tseries_altitude[i][1], label = list(test_tles.keys())[i] + 'real')
        ax.set_xlabel('Time MJD')
        ax.set_ylabel('Altitude (km)')
        ax.set_title('Altitude of Satellites Over Time')
        ax.legend()
        plt.show()

    pass

def test_spaceobject_creation():
    """
    Test that when an object is instantiated, the correct attributes are set.
    In particular, the characteristic area and mass are imputed to some value if they are not specified or 0. 
    """
    object_1 = SpaceObject(  object_type='+', 
                                    payload_operational_status='Active',
                                    application="Unknown", 
                                    operator='Starink', 
                                    mass=0, 
                                    eccentricity=0.1, 
                                    inc=2*np.pi, 
                                    argp=12, 
                                    raan=252, 
                                    source = "Guatemala",
                                    launch_date='2012-10-12', 
                                    decay_date='2013-10-12',
                                    rso_name='bigboisat',
                                    perigee='1200',
                                    apogee='1000',
                                    tle="1 53544U 22101T   23122.20221856  .00001510  00000-0  11293-3 0  9999\n2 53544  53.2176  64.0292 0001100  79.8127 280.2989 15.08842383 38928")
    
if __name__ == "__main__":
    test_sgp4_drag()
    # test_spaceobject_creation()
