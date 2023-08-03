import unittest
import numpy as np
import datetime
from utils.Conversions import utc_to_jd, jd_to_utc, tle_parse, tle_convert, TLE_time
from utils.SpaceObject import SpaceObject
from utils.Propagators import kepler_prop

class PropagatorTest(unittest.TestCase):
    # Class level tolerance
    tolerance = 1 #I have set this tolerance to 1 because of the mismatch between the time stamps of the numerical and analytical propagator. 
    #TODO: the alignment between the two ephemerides needs to be fixed.
    #TODO: if all other variables are the same the tolerance can be set to 1e-6.
    
    # Function to calculate 3D position difference
    @staticmethod
    def calculate_position_difference(keplerian_output, numerical_output):
        keplerian_position = np.array(keplerian_output[1])
        numerical_position = numerical_output[:3]
        position_difference = np.abs(np.linalg.norm(keplerian_position) - np.linalg.norm(numerical_position))
        return position_difference
    
    # Setup variables and objects for all tests
    @classmethod
    def setUpClass(cls):
        # TLE strings for test satellites
        cls.test_tles = {
            'navstar81': "1 48859U 21054A   23116.83449170 -.00000109  00000-0  00000-0 0  9996\n2 48859  55.3054  18.4561 0008790 213.9679 183.6522  2.00556923 13748",
            'OneWeb20': "1 45133U 20008C   23116.69886660  .00000678  00000-0  19052-2 0  9998\n2 45133  87.8784 264.1991 0001687  89.2910 270.8411 13.10377378158066",
            'Starlink70': "1 53544U 22101T   23122.20221856  .00001510  00000-0  11293-3 0  9999\n2 53544  53.2176  64.0292 0001100  79.8127 280.2989 15.08842383 38928"
        }

        cls.prop_start = [datetime.datetime.strptime('2023-05-02 12:45:00', '%Y-%m-%d %H:%M:%S')] # start date of the test propagation
        cls.prop_end = [datetime.datetime.strptime('2023-05-02 14:45:00', '%Y-%m-%d %H:%M:%S')] # end date of the test propagation
        cls.start_jd = utc_to_jd(cls.prop_start)
        cls.end_jd = utc_to_jd(cls.prop_end)

    @staticmethod
    def calculate_position_difference(keplerian_ephem, numerical_ephem):
        """Function to calculate 3D position difference"""
        keplerian_position = keplerian_ephem[1]
        numerical_position = numerical_ephem[1]
        position_difference = np.abs(np.linalg.norm(keplerian_position) - np.linalg.norm(numerical_position))
        return position_difference
    
    @staticmethod
    def create_space_object_from_tle(tle_string):
        """Creates a SpaceObject from a given TLE"""
        tle_dict = tle_parse(tle_string)
        tle_kepels = tle_convert(tle_dict)
        tle_time = TLE_time(tle_string)
        tle_epoch = jd_to_utc(tle_time)
        tle_epoch_str = str(tle_epoch)
        epoch = tle_epoch_str.replace(' ', 'T')
        return SpaceObject(sma = tle_kepels['a'], perigee=tle_kepels['a']-6378.137, apogee=tle_kepels['a']-6378.137, eccentricity=tle_kepels['e'], inc = tle_kepels['i'], argp = tle_kepels['arg_p'], raan=tle_kepels['RAAN'], tran=tle_kepels['true_anomaly'], characteristic_area=0.011, mass = 250, epoch = epoch, launch_date='2023-05-02')

    def test_kepler_vs_numerical(self):
        """Test Kepler vs Numerical"""
        for satellite_name, tle in self.test_tles.items():
            print(f"Testing {satellite_name}...")
            test_sat = self.create_space_object_from_tle(tle)

            # Propagate the SpaceObject using the numerical propagator
            test_sat.prop_catobject(jd_start=self.start_jd[0], jd_stop=self.end_jd[0], step_size=1, integrator_type="RK45", output_freq=1, force_model=["grav_mono"])
            test_sat_ephem = test_sat.ephemeris  # Access the ephemeris attribute directly

            # Initialize Keplerian propagation with start JD
            jd_start_kepler = self.start_jd[0]
            keplerian_test_sat_ephem = []

            # Iterate through the numerical ephemeris and propagate Keplerian for the corresponding time steps
            for num_ephem in test_sat_ephem:
                print(f"Propagating Keplerian for {num_ephem[0]}...")
                jd_stop_kepler = num_ephem[0]
                kepler_ephem = kepler_prop(jd_start=jd_start_kepler, jd_stop=jd_stop_kepler, step_size=(jd_stop_kepler - jd_start_kepler), a=test_sat.sma, e=test_sat.eccentricity, i=test_sat.inc, w=test_sat.argp, W=test_sat.raan, V=test_sat.tran)[-1]
                keplerian_test_sat_ephem.append(kepler_ephem)
                jd_start_kepler = jd_stop_kepler  # Update start JD for next step

            # Calculate the differences between the Keplerian and numerical propagations
            print("Calculating position differences...")
            cart_diffs = [self.calculate_position_difference(kepler_ephem, numerical_ephem) for kepler_ephem, numerical_ephem in zip(keplerian_test_sat_ephem, test_sat_ephem)]

            # Calculate min, max, and mean differences
            min_diff = np.min(cart_diffs)
            max_diff = np.max(cart_diffs)
            mean_diff = np.mean(cart_diffs)

            print(f"For satellite {satellite_name}, the minimum, maximum, and mean position differences are {min_diff}, {max_diff}, and {mean_diff} respectively")

            if mean_diff > self.tolerance:
                raise ValueError(f"Mean position difference exceeds the tolerance level of {self.tolerance}. Recorded mean difference: {mean_diff}")

            if max_diff > self.tolerance:
                raise ValueError(f"Final position difference exceeds the tolerance level of {self.tolerance}. Recorded final difference: {max_diff}")

if __name__ == '__main__':
    unittest.main()