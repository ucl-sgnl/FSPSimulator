#Only way i could rin this was to bring it into src and run using "python -m src.tests.PropagatorTest"

import numpy as np
import datetime
import matplotlib.pyplot as plt
from utils.Conversions import tle_parse, tle_convert, TLE_time, jd_to_utc, utc_to_jd
from utils.Propagators import sgp4_prop_TLE, numerical_prop, kepler_prop
from utils.SpaceObject import SpaceObject

#Satellites used in validation

#navstar81, pulled on 27Apr2023
test_tle1 = "1 48859U 21054A   23116.83449170 -.00000109  00000-0  00000-0 0  9996\n2 48859  55.3054  18.4561 0008790 213.9679 183.6522  2.00556923 13748"
#OneWeb20, pulled on 27Apr2023
test_tle2 = "1 45133U 20008C   23116.69886660  .00000678  00000-0  19052-2 0  9998\n2 45133  87.8784 264.1991 0001687  89.2910 270.8411 13.10377378158066"
#Starlink70, pulled on 27Apr2023
test_tle3 = "1 53544U 22101T   23122.20221856  .00001510  00000-0  11293-3 0  9999\n2 53544  53.2176  64.0292 0001100  79.8127 280.2989 15.08842383 38928"

#dictionary of TLEs
test_tles = {'navstar81': test_tle1, 'OneWeb20': test_tle2, 'Starlink70': test_tle3}

prop_start = [datetime.datetime.strptime('2023-05-02 12:45:00', '%Y-%m-%d %H:%M:%S')]
prop_end = [datetime.datetime.strptime('2023-05-05 01:48:00', '%Y-%m-%d %H:%M:%S')]
start_jd = utc_to_jd(prop_start)
end_jd = utc_to_jd(prop_end)



# Function to calculate 3D position difference
def calculate_position_difference(keplerian_output, numerical_output):
    keplerian_position = np.array(keplerian_output[1])
    numerical_position = numerical_output[:3]
    position_difference = np.abs(keplerian_position - numerical_position)
    return position_difference


def kepler_vs_numerical_test():
    # Propagate an object using the keplerian propgator
    # Propgate the same orbit using the numerical propagator and a wide range of time steps (10-1000 seconds in steps of 10 seconds)
    # record the RMS of the position and velocity at each time step for the numericall propagated orbits agains the keplerian propagated orbit
    
    step_sizes = np.arange(10, 200, 10)

    for sat in test_tles:
        print("Propagating numerically: ", sat)
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
        
        cart_diffs = []

        for step_size in step_sizes:
            #Propagate the SpaceObject using the keplerian propagator
            test_sat.prop_catobject(jd_start=start_jd[0], jd_stop=end_jd[0], step_size=step_size, integrator_type="RK45", output_freq = 1000, force_model = ["grav_mono"])
            test_sat_ephem = test_sat.ephemeris
            keplerian_test_sat_ephem = kepler_prop(jd_start = start_jd[0], jd_stop=end_jd[0], step_size=step_size, a = tle_kepels['a'], e = tle_kepels['e'], i = tle_kepels['i'], w = tle_kepels['arg_p'], W = tle_kepels['RAAN'], V = tle_kepels['true_anomaly'])
            print("numerical len and keplerian len:", len(test_sat_ephem), len(keplerian_test_sat_ephem))
            cart_diff = []
            for i in range(len(test_sat_ephem)):
                cart_diff.append(calculate_position_difference(keplerian_test_sat_ephem[i], test_sat_ephem[i]))
            cart_diffs.append(cart_diff)

            step_size_in_days = step_size / (24 * 60 * 60) # convert step_size from seconds to days
            jd_time_stamps = np.linspace(start_jd[0], end_jd[0], num=len(cart_diff))

            print("jd_time_stamps shape: ", np.shape(jd_time_stamps))
            print("cart_diff shape: ", np.array(cart_diff).shape)

            #plot the RMS of the position difference against the time stamps for each step size
            plt.plot(jd_time_stamps, np.sqrt(np.sum(np.square(cart_diff), axis=1)), label = f"Step size: {step_size} seconds")
            plt.xlabel("Julian Date")
            plt.ylabel("RMS of position difference (km)")
            plt.title(f"RMS of position difference between keplerian and numerical propagators for {sat}")
            plt.legend()
            plt.show()


        

if __name__ == "__main__":
    kepler_vs_numerical_test()
