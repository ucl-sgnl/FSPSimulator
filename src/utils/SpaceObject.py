### class for an individual space object
import uuid
import math 
import datetime
import numpy as np
import sgp4
import matplotlib.pyplot as plt
from sgp4.api import Satrec, WGS72, CustomSatrec
from src.utils.coords import kep2car, trueanom2meananom, calculate_kozai_mean_motion, expo_simplified, utc_to_jd, tle_parse, tle_convert, sgp4_prop_TLE, build_tle
from sgp4.api import Satrec, WGS72, CustomSatrec 
import matplotlib.cm as cm

class SpaceObject:
    def __init__(self, cospar_id=None, rso_name=None, rso_type=None, payload_operational_status=None, orbit_type=None, application=None, source=None, 
                 orbital_status_code=None, launch_site=None, mass=None, maneuverable=False, spin_stabilized=False, 
                 orbital_period=None, object_type = None, apogee_altitude=None, perigee_altitude=None, radar_cross_section=None, 
                 characteristic_area=None, characteristic_length=None, propulsion_type=None, epoch=None, sma=None, inc=None, 
                 argp=None, raan=None, tran=None, eccentricity=None, operator=None, launch_date=None,  decay_date=None, tle = None):
        
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
        self.tle = None
        self.sgp4_ephemeris = None
        self.source = source # http://celestrak.org/satcat/sources.php #TODO: do we really want source? This is just the country of origin, but I think maybe owner is more relevant nowadays
        # self.orbital_status_code = str(orbital_status_code) #TODO: I also would argue this is not that useful. We have payload_operational_status (especially if focus is just Earth orbit this is totally redundant i believe)
        self.launch_site = str(launch_site)
        self.decay_date = datetime.datetime.strptime(decay_date, '%Y-%m-%d') if (self.decay_date is not None and self.decay_date != '-') else None
        self.mass = float(mass) if mass is not None else None #in Kg
        self.maneuverable = str(maneuverable)
        self.spin_stabilized = str(spin_stabilized)# I have left spin_stabilized and maneuverable as strings in case we later want to add more options than just True/False (for example different thruster types resulting in different kinds of maneuverability)
        self.orbital_period = float(orbital_period) if orbital_period is not None else None #in minutes #TODO: I don't think we need this, can be calculated from sma
        self.apogee_altitude = float(apogee_altitude) #in Km
        self.perigee_altitude = float(perigee_altitude) #in Km
        self.radar_cross_section = float(radar_cross_section) if radar_cross_section is not None else None #in meters^2
        self.characteristic_area = float(characteristic_area) if characteristic_area is not None else None #in meters^2
        self.characteristic_length = float(characteristic_length) if characteristic_length is not None else None #in meters
        self.propulsion_type = str(propulsion_type)
        if epoch is None:
            self.epoch = None #epoch must be cast to datetime object and be specified in UTC time in the format: datetime(year-month-day hour:minute:second)
        else:
            self.epoch = datetime.datetime.strptime(epoch, '%Y-%m-%d %H:%M:%S') #in UTC
        self.sma = (self.apogee_altitude + self.perigee_altitude)/2 + 6378.137 #in km

        self.inc = float(inc)
        self.argp = float(argp) if sma is not None else None
        self.raan = float(raan)
        self.tran = float(tran) if tran is not None else None
        self.eccentricity = float(eccentricity)
        if (self.tran is not None):
            self.meananomaly = trueanom2meananom(self.tran, self.eccentricity)

        # These are attributes that are not required to be specified on instantiation, but are to be computed later on
        #mean anomaly to be computed using trueanom2meananom (from true anomaly)
        self.cart_state = None #cartesian state vector [x,y,z,u,v,w] to be computed using generate_cart (from keplerian elements)
        self.C_d = 2.2 #Drag coefficient
        # self.orbit_type = orbit_type #TODO: i think this is redundant on instantiation, we can calculate this from altitude and inclination. I wrote the function just call it here
        self.tle = tle if tle is not None else None
        # this cannot be used currently as it is not set up to validate JSR's catalogue
        #self._validate_types()

        #Calculate TLE components
        altitude = (self.perigee_altitude+self.apogee_altitude)/2 #TODO: make this not allowed for non-cirular orbits
        density = self.get_density(altitude, model = "exponential")
        self.bstar = (self.C_d * self.characteristic_area * density)/2*self.mass #BStar = Cd * A * rho / 2m. Where Cd is the drag coefficient, A is the cross-sectional area of the satellite, rho is the density of the atmosphere, and m is the mass of the satellite.
        self.ndot = 0.0
        self.nddot = 0.0
        self.meananomaly = trueanom2meananom(self.tran, self.eccentricity)
        self.no_kozai = calculate_kozai_mean_motion(a = self.sma, mu = 398600.4418)
        self.satnum = 0 #TODO: update this placeholder value with NORAD ID. Does not affect result of SGP4 propagation
        self.sgp4epoch = self.sgp4_epoch() #SGP4 epoch is the number of days since 1949 December 31 00:00 UT

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
    
    def get_density(self, altitude, model = "exponential"):
        if model == "exponential":
            return expo_simplified(altitude)
        #TODO: other density models here when ready (USSA 76 probably only one we need)
        else:
            return 1e-12 #Placeholder value #in kg/m^3
        
    def sgp4_prop_catobjects(self, jd_start, jd_stop, step_size):
        #TODO: make it so that this is only applied if a TLE is not already available
        if self.tle is None:
            print("No TLE available for this object. Generating one from the catalog data")
            # generate a TLE from the data in the catalog object
            new_TLE = build_tle(
                int(25654),  # catalog_number
                "U",  # classification
                int(22),  # launch_year
                int(10),  # launch_number
                'AAA',  # launch_piece
                self.epoch.year % 100,  # epoch_year
                self.epoch.day,  # epoch_day
                self.ndot,  # first_derivative
                self.nddot,  # second_derivative
                self.bstar,  # drag_term
                0,  # ephemeris_type
                0,  # element_set_number
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
            print("new TLE: ", self.tle)
        # propagate the TLE and return the ephemeris 
        self.sgp4_ephemeris = sgp4_prop_TLE(self.tle, jd_start, jd_stop, step_size)


def test_sgp4_prop():

    startday = [datetime.datetime.strptime('2023-04-26 00:00:00', '%Y-%m-%d %H:%M:%S')]
    enday = [datetime.datetime.strptime('2023-04-26 00:00:00', '%Y-%m-%d %H:%M:%S')]
    start_jd = utc_to_jd(startday)
    end_jd = utc_to_jd(enday)

    # Random 3LE from SpaceTrack SATELLITE
    valid_tle = "SPACEX\n 1 44271U 19029AN  20288.57092606  .06839568  12140-4  13629-2 0  9994\n2 44271  52.9879 323.6967 0003539  53.2438  81.7998 16.30723255 78035 "

    #make my own TLE
    tle_dict = tle_parse(valid_tle)
    tle_kepels = tle_convert(tle_dict)
    test_sat = SpaceObject(sma = tle_kepels['a'], perigee_altitude=tle_kepels['a']-6378.137, apogee_altitude=tle_kepels['a']-6378.137, eccentricity=tle_kepels['e'], inc = tle_kepels['i'], argp = tle_kepels['arg_p'], raan=tle_kepels['RAAN'], tran=tle_kepels['true_anomaly'], characteristic_area=20.0, mass = 100, epoch = '2023-04-26 00:00:00')
    test_sat.sgp4_prop_catobjects(start_jd, end_jd, 600)
    test_sat_ephem = test_sat.sgp4_ephemeris
    print("testsat ephem:", test_sat_ephem)
    ########## TLE (not 3LE) ##########
    valid_tle = " 1 44271U 19029AN  20288.57092606  .06839568  12140-4  13629-2 0  9994\n2 44271  52.9879 323.6967 0003539  53.2438  81.7998 16.30723255 78035 "
    valid_tle_ephem = sgp4_prop_TLE(valid_tle, jd_start=start_jd, jd_end=end_jd, dt=600)
    print("valid tle ephemeris: ", valid_tle_ephem)

    # #plot ephemeris
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # ax.scatter(test_ephemeris[:,0], test_ephemeris[:,1],test_ephemeris[:,2], c=test_date_range, s=1)
    # #force aspect ratio to be 1:1:1
    # ax.set_xlim(-7000, 7000)
    # ax.set_ylim(-7000, 7000)
    # ax.set_zlim(-7000, 7000)
    # ax.legend()
    # #add colorbar
    # m = cm.ScalarMappable(cmap=cm.jet)
    # m.set_array(test_date_range)
    # fig.colorbar(m)
    # plt.show()

if __name__ == "__main__":
    test_sgp4_prop()
# 1 25654U 22010AAA 23026.00000000 -.000021    00000-0 +84175-7 0 00009
# 1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927


# 2 25654 000.9248 005.6496 0003539 000.9293 000.1584 00.07115372 000000
# 2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537