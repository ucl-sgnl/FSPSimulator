### class for an individual space object
import uuid
import math 
import datetime
import numpy as np
import sgp4
import matplotlib.pyplot as plt
from sgp4.api import Satrec, WGS72
from src.utils.coords import kep2car, true_to_mean_anomaly, calculate_kozai_mean_motion, expo_simplified, utc_to_jd, tle_parse, tle_convert, sgp4_prop_TLE, build_tle, orbital_period, get_day_of_year_and_fractional_day, TLE_time, jd_to_utc
import matplotlib.cm as cm

class SpaceObject:
    def __init__(self, cospar_id=None, rso_name=None, rso_type=None, payload_operational_status=None, orbit_type=None, application=None, source=None, 
                 orbital_status_code=None, launch_site=None, mass=None, maneuverable=False, spin_stabilized=False, 
                object_type = None, apogee_altitude=None, perigee_altitude=None, radar_cross_section=None, 
                 characteristic_area=None, characteristic_length=None, propulsion_type=None, epoch=None, sma=None, inc=None, 
                 argp=None, raan=None, tran=None, eccentricity=None, operator=None, launch_date=None,  decay_date=None, tle = None, station_keeping=None):
        
        self.CATID = None # This is a computed property based on the GUID and isn't included as a parameter
        self.GUID = self._compute_catid() # This will be set by the system and isn't included as a parameter
        self.launch_date = launch_date
        self.decay_date = decay_date
        self.cospar_id = str(cospar_id)
        self.rso_name = str(rso_name)
        self.rso_type = str(rso_type) if rso_type is not None else None
        self.payload_operational_status = str(payload_operational_status)
        self.object_type = str(object_type)
        self.application = str(application)
        self.operator = str(operator)
        # if any of the following if None then we call the function to compute them]
        if characteristic_area is None:
            self.impute_char_area()
        else: 
            self.characteristic_area = float(characteristic_area)
        if characteristic_length is None:
            self.impute_char_length()
        else:
            self.characteristic_length = float(characteristic_length)
        if mass is None:
            self.impute_mass()
        else:
            self.mass = float(mass)
        self.source = source # http://celestrak.org/satcat/sources.php #TODO: do we really want source? This is just the country of origin, but I think maybe owner is more relevant nowadays
        # self.orbital_status_code = str(orbital_status_code) #TODO: I also would argue this is not that useful. We have payload_operational_status (especially if focus is just Earth orbit this is totally redundant i believe)
        self.launch_site = str(launch_site)
        self.decay_date = datetime.datetime.strptime(decay_date, '%Y-%m-%d') if (self.decay_date is not None and self.decay_date != '-') else None
        self.maneuverable = str(maneuverable)
        self.spin_stabilized = str(spin_stabilized)# I have left spin_stabilized and maneuverable as strings in case we later want to add more options than just True/False (for example different thruster types resulting in different kinds of maneuverability)
        self.apogee_altitude = float(apogee_altitude) #in Km
        self.perigee_altitude = float(perigee_altitude) #in Km
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
        #Variables for Propagation
        self.tle = tle if tle is not None else None
        # ballisitic coefficient. If none this is replaced with 0.00000140 #TODO: find correct value for this
        self.ndot = 0.00000140
        self.nddot = 0.0 # If none this is set to 0.0
        self.sgp4_ephemeris = None # This is where the ephemeris of the SGP4 propagation will be stored
        self.sma = (self.apogee_altitude + self.perigee_altitude)/2 + 6378.137 #in km
        self.orbital_period = orbital_period(self.sma) #in minutes
        self.inc = float(inc)
        self.argp = float(argp) if sma is not None else None
        self.raan = float(raan)
        self.tran = float(tran) if tran is not None else None
        self.eccentricity = float(eccentricity)
        self.tran = np.random.uniform(0, 2*np.pi) if self.tran is None else self.tran
        self.meananomaly = true_to_mean_anomaly(self.tran, self.eccentricity)
            # this can be 0 as it doesn't really change much for the time scale we are looking at
        self.altitude = (self.perigee_altitude+self.apogee_altitude)/2 #TODO: make this not allowed for non-cirular orbits
        # self.atmos_density = 1e-12
        self.get_atmospheric_density(model = "exponential")            #BStar = rho0 #TODO: FIX
        self.C_d = 2.2 #Drag coefficient
        self.bstar = -(self.C_d * self.characteristic_area * self.atmos_density)/2*self.mass #BStar = Cd * A * rho / 2m. Where Cd is the drag coefficient, A is the cross-sectional area of the satellite, rho is the density of the atmosphere, and m is the mass of the satellite.
        self.no_kozai = calculate_kozai_mean_motion(a = self.sma, mu = 398600.4418)
        self.sgp4epoch = self.sgp4_epoch() #SGP4 epoch is the number of days since 1949 December 31 00:00 UT

        #Variables that may be populated by functions during propagation
        self.cart_state = None #cartesian state vector [x,y,z,u,v,w] to be computed using generate_cart (from keplerian elements)
        
        # this cannot be used currently as it is not set up to validate Celestrak/Spacetrack data
        # self._validate_types()

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
    
    def impute_char_length(self):
        """
        This function will impute a characteristic length based on the object type
        """
        if self.object_type == 'PAYLOAD':
            self.characteristic_length = 1.3
        elif self.object_type == "ROCKET BODY":
            self.characteristic_length = 8
        else:
            self.characteristic_length = 1.3

    def impute_char_area(self):
        """
        This function will impute a characteristic area based on the object type
        """
        if self.object_type == 'PAYLOAD':
            self.characteristic_area = 1.3 * 3
        elif self.object_type == "ROCKET BODY":
            self.characteristic_area = 8 * 3
        else:
            self.characteristic_area = 1.3 * 3

    def impute_mass(self):
        """
        This function will impute a mass based on the object type
        """
        # From Alfano, Oltrogge and Sheppard, 2020
        self.mass = 0.1646156590 * self.characteristic_length ** (-1.0554951607)


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
            self.atmos_density = expo_simplified(self.altitude)*1e6 #fix units to kg/m^3 
        #TODO: other density models here when ready (USSA 76 probably only one we need)
        else:
            self.atmos_density = 1e-12 #Placeholder value #in kg/m^3
        
    def prop_catobjects(self, jd_start, jd_stop, step_size):
        #TODO: make it so that we can specify a time window when station keeping is set to be active. Then the rest of the time the object just propagates using SGP4
        if self.station_keeping == True:
            self.ephemeris = kepler_prop(jd_start, jd_stop, step_size, a=self.sma, e=self.eccentricity, i=self.inc, w=self.argp, W=self.raan, V=self.tran)
        else: 
            if self.tle is None:
                #print("No TLE available for this object. Generating one from the catalog data")
                # generate a TLE from the data in the catalog object
                new_TLE = build_tle(
                    int(12345),  # catalog_number -> this is a placeholder value and does not affect the propagation
                    "U",  # classification -> also a placeholder value
                    int(self.epoch.year % 100),  # launch_year -> also a placeholder value
                    int(29),  # launch_number -> also a placeholder value
                    'AN',  # launch_piece -> also a placeholder value
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
                self.tle = new_TLE
            # propagate the TLE and return the ephemeris 
            self.ephemeris = sgp4_prop_TLE(self.tle, jd_start, jd_stop, step_size)

def test_sgp4_prop():

    #navstar81, pulled on 27Apr2023
    test_tle1 = "1 48859U 21054A   23116.83449170 -.00000109  00000-0  00000-0 0  9996\n2 48859  55.3054  18.4561 0008790 213.9679 183.6522  2.00556923 13748"
    #OneWeb20, pulled on 27Apr2023
    test_tle2 = "1 45133U 20008C   23116.69886660  .00000678  00000-0  19052-2 0  9998\n2 45133  87.8784 264.1991 0001687  89.2910 270.8411 13.10377378158066"
    #Starlink70, pulled on 27Apr2023
    test_tle3 = "1 56135U 23046AV  23117.00001157  .00046284  00000-0  35489-3 0  9999\n2 56135  43.0006  94.0710 0001704 259.7544 135.6794 15.71324535  4975"

    #dictionary of TLEs
    test_tles = {'navstar81': test_tle1, 'OneWeb20': test_tle2, 'Starlink70': test_tle3}

    prop_start = [datetime.datetime.strptime('2023-01-24 12:45:00', '%Y-%m-%d %H:%M:%S')]
    prop_end = [datetime.datetime.strptime('2024-04-27 01:48:00', '%Y-%m-%d %H:%M:%S')]
    start_jd = utc_to_jd(prop_start)
    end_jd = utc_to_jd(prop_end)
    t_step = 60*60*24 #1 day in seconds

    #test propagation of TLEs in test_tles
    RMS_dict = {} #for each satellite, this will be a list of RMS position differences between the real and fabricated TLEs
    tseries_3D = [] #for each satellite, this will be a list of 3D position differences between the real and fabricated TLEs
    tseries_altitude = [] #for each satellite, this will contain two lists (one for the real and one for the fabricated TLE) of altitude values
    t_series_pos = [] #for each satellite, this will be a list of position differences between the real and fabricated TLEs
    altitude_diffs = [] #for each satellite, this will be a list of altitude differences between the real and fabricated TLEs
    inclination_diffs = [] #for each satellite, this will be a list of inclination differences between the real and fabricated TLEs

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
        test_sat = SpaceObject(sma = tle_kepels['a'], perigee_altitude=tle_kepels['a']-6378.137, apogee_altitude=tle_kepels['a']-6378.137, eccentricity=tle_kepels['e'], inc = tle_kepels['i'], argp = tle_kepels['arg_p'], raan=tle_kepels['RAAN'], tran=tle_kepels['true_anomaly'], characteristic_area=2.5, mass = 150, epoch = epoch, station_keeping=True)
        test_sat.prop_catobjects(start_jd[0], end_jd[0], t_step)
        test_sat_ephem = test_sat.sgp4_ephemeris
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
        print("valid TLE:", valid_tle)
        valid_tle_pos = [x[1] for x in valid_tle_ephem]
        valid_tle_pos = np.array(valid_tle_pos)
        valid_tle_times = [x[0] for x in valid_tle_ephem]
        valid_tle_altitudes = [np.linalg.norm(x)-6378.137 for x in valid_tle_pos]
        sat_altitude.append(valid_tle_altitudes)
        sat_pos.append(valid_tle_pos)
        tseries_altitude.append(sat_altitude)
        t_series_pos.append(sat_pos)
        ### calculate the RMS position differences time series
        # diff_pos = [np.linalg.norm(test_pos[i] - valid_tle_pos[i]) for i in range(len(test_pos))]
        # diff_pos = np.array(diff_pos)
        # tseries_3D.append(diff_pos)
        # RMS = np.sqrt(np.mean(diff_pos**2))
        # RMS_dict[sat] = RMS
        
        ### calculate the altitude differences time series
        sat_altitude_diffs = [sat_altitude[0][i] - sat_altitude[1][i] for i in range(len(sat_altitude[0]))]
        altitude_diffs.append(sat_altitude_diffs)

#plot the altitude of each satellite over time
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i in range(1):
        ax.plot(tseries_altitude[i][0], label = list(test_tles.keys())[i] + ' real')
        ax.plot(tseries_altitude[i][1], label = list(test_tles.keys())[i] + ' fabricated')
    ax.set_xlabel('Time MJD')
    ax.set_ylabel('Altitude (km)')
    ax.set_title('Altitude of Satellites Over Time')
    ax.legend()
    plt.show()


   #plot altitude differences time series
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i in range(len(altitude_diffs)):
        ax.plot(altitude_diffs[i], label = list(test_tles.keys())[i])
    ax.set_xlabel('Time MJD')
    ax.set_ylabel('Altitude Difference (km)')
    ax.set_title('Altitude Difference Between Real and Fabricated TLEs')    
    ax.legend()
    plt.show()


    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # for i in range(len(t_series_pos)):
    #     ax.scatter(t_series_pos[i][0][:,0], t_series_pos[i][0][:,1], t_series_pos[i][0][:,2], label = list(test_tles.keys())[i] + ' real')
    #     ax.scatter(t_series_pos[i][1][:,0], t_series_pos[i][1][:,1], t_series_pos[i][1][:,2], label = list(test_tles.keys())[i] + ' fabricated')
    # ax.set_xlabel('X (km)')
    # ax.set_ylabel('Y (km)')
    # ax.set_zlabel('Z (km)')
    # ax.set_title('3D Position Difference Between Real and Fabricated TLEs')
    # ax.legend()
    # plt.show()
    
    # print("RMS values for each satellite: ", RMS_dict)

    # for i in range(len(tseries_3D)):
    #     plt.plot(test_times, tseries_3D[i], label = list(test_tles.keys())[i])
    # plt.xlabel('Time (Julian Date)')
    # plt.ylabel('3D Position Difference (km)')
    # plt.title('3D Position Difference Between Real and Fabricated TLEs')
    # plt.legend()
    # plt.show()
    #
    # for i in range(len(tseries_altitude)):
    #     plt.plot(test_times, tseries_altitude[i][0], label = list(test_tles.keys())[i] + ' real')
    #     plt.plot(test_times, tseries_altitude[i][1], label = list(test_tles.keys())[i] + ' fabricated')
    # plt.xlabel('Time (Julian Date)')
    # plt.ylabel('Altitude (km)')
    # plt.title('Altitude of Real and Fabricated TLEs')
    # plt.legend()
    plt.show()

if __name__ == "__main__":
    test_sgp4_prop()
