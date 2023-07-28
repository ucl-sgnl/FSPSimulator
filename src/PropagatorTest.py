#Only way i could rin this was to bring it into src and run using "python -m src.tests.PropagatorTest"
import numpy as np
import datetime
import matplotlib.pyplot as plt
from utils.Conversions import tle_parse, tle_convert, TLE_time, jd_to_utc, utc_to_jd
from utils.Propagators import sgp4_prop_TLE, numerical_prop, kepler_prop
from utils.SpaceObject import SpaceObject

#Satellites used in validation

#Doing 1 MEO + 2 LEOs
#navstar81, pulled on 27Apr2023
test_tle1 = "1 48859U 21054A   23116.83449170 -.00000109  00000-0  00000-0 0  9996\n2 48859  55.3054  18.4561 0008790 213.9679 183.6522  2.00556923 13748"
#OneWeb20, pulled on 27Apr2023
test_tle2 = "1 45133U 20008C   23116.69886660  .00000678  00000-0  19052-2 0  9998\n2 45133  87.8784 264.1991 0001687  89.2910 270.8411 13.10377378158066"
#Starlink70, pulled on 27Apr2023
test_tle3 = "1 53544U 22101T   23122.20221856  .00001510  00000-0  11293-3 0  9999\n2 53544  53.2176  64.0292 0001100  79.8127 280.2989 15.08842383 38928"

#dictionary of TLEs
test_tles = {'navstar81': test_tle1, 'OneWeb20': test_tle2, 'Starlink70': test_tle3}

prop_start = [datetime.datetime.strptime('2023-05-02 12:45:00', '%Y-%m-%d %H:%M:%S')]
prop_end = [datetime.datetime.strptime('2023-08-02 12:45:00', '%Y-%m-%d %H:%M:%S')]
start_jd = utc_to_jd(prop_start)
end_jd = utc_to_jd(prop_end)

# Function to calculate 3D position difference
def calculate_position_difference(keplerian_output, numerical_output):
    keplerian_position = np.array(keplerian_output[1])
    numerical_position = numerical_output[:3]
    position_difference = np.abs(np.linalg.norm(keplerian_position) - np.linalg.norm(numerical_position))
    return position_difference

def kepler_vs_numerical_test_plot():
    # Propagate an object using the keplerian propgator
    # Propagate the same orbit using the numerical propagator and a wide range of time steps (10-1000 seconds in steps of 10 seconds)
    # Record the RMS of the position and velocity at each time step for the numerically propagated orbits against the keplerian propagated orbit
    
    step_sizes = [10000000]

    for sat in test_tles:
        print("Propagating: ", sat)

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
        cart_diff_means = []

        for step_size in step_sizes:
            #Propagate the SpaceObject using the keplerian propagator
            print(f"Propagating numerically with step size {step_size} seconds")
            test_sat.prop_catobject(jd_start=start_jd[0], jd_stop=end_jd[0], step_size=step_size, integrator_type="RK45", output_freq = step_size, force_model = ["grav_mono"])
            test_sat_ephem = test_sat.ephemeris
            keplerian_test_sat_ephem = kepler_prop(jd_start = start_jd[0], jd_stop=end_jd[0], step_size=step_size, a = tle_kepels['a'], e = tle_kepels['e'], i = tle_kepels['i'], w = tle_kepels['arg_p'], W = tle_kepels['RAAN'], V = tle_kepels['true_anomaly'])
            cart_diff = []
            for i in range(len(test_sat_ephem)):
                cart_diff.append(calculate_position_difference(keplerian_test_sat_ephem[i], test_sat_ephem[i]))

            print(f"Initial position difference at start of prop {step_size} is: {cart_diff[0]}")
            print(f"Final position difference at end of prop {step_size} is: {cart_diff[-1]}")
            print(f"Mean position difference for prop {step_size} is: {np.mean(cart_diff)}")
            print(f"variance of position difference for prop {step_size} is: {np.var(cart_diff)}")
            cart_diffs.append(cart_diff)
            cart_diff_means.append(np.mean(cart_diff))

        for i, step_size in enumerate(step_sizes):
            jd_time_stamps = np.linspace(start_jd[0], end_jd[0], num=len(cart_diffs[i]))

            #plot the RMS of the position difference against the time stamps for each step size
            plt.scatter(jd_time_stamps, cart_diffs[i], label = f"Step size: {step_size} seconds", s=1, alpha=0.2)
        
        plt.xlabel("Julian Date")
        plt.ylabel("3D position difference (km)")
        plt.title(f"Norm of position difference between keplerian and numerical propagators for {sat}")
        plt.yscale('log')  # This line sets the y-axis to a logarithmic scale

        plt.legend()
        plt.show()

if __name__ == "__main__":
    kepler_vs_numerical_test_plot()

import numpy as np
import warnings

def kepler_vs_numerical_test(tolerance):
    # Propagate an object using the keplerian propgator
    # Propgate the same orbit using the numerical propagator and a wide range of time steps (10-1000 seconds in steps of 10 seconds)
    # record the RMS of the position and velocity at each time step for the numericall propagated orbits agains the keplerian propagated orbit

    step_sizes = np.arange(10, 20000, 1000)

    for sat in test_tles:
        print("Propagating numerically: ", sat)

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
            test_sat.prop_catobject(jd_start=start_jd[0], jd_stop=end_jd[0], step_size=step_size, integrator_type="RK45", output_freq = step_size, force_model = ["grav_mono"])
            test_sat_ephem = test_sat.ephemeris
            keplerian_test_sat_ephem = kepler_prop(jd_start = start_jd[0], jd_stop=end_jd[0], step_size=step_size, a = tle_kepels['a'], e = tle_kepels['e'], i = tle_kepels['i'], w = tle_kepels['arg_p'], W = tle_kepels['RAAN'], V = tle_kepels['true_anomaly'])
            cart_diff = []
            for i in range(len(test_sat_ephem)):
                cart_diff.append(calculate_position_difference(keplerian_test_sat_ephem[i], test_sat_ephem[i]))

            cart_diffs.append(cart_diff)
        
        # convert to numpy array for easy calculation
        cart_diffs = np.array(cart_diffs)

        min_diff = np.min(cart_diffs)
        max_diff = np.max(cart_diffs)
        mean_diff = np.mean(cart_diffs)

        print(f"For satellite {sat}, the minimum, maximum, and mean position differences are {min_diff}, {max_diff}, and {mean_diff} respectively")

        if mean_diff > tolerance:
            warnings.warn("Mean position difference exceeds the tolerance level!")
