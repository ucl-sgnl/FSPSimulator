#Only way i could rin this was to bring it into src and run using "python -m src.tests.PropagatorTest"

import numpy as np
import matplotlib.pyplot as plt
from fspsim.utils.Conversions import tle_parse, tle_convert, TLE_time, jd_to_utc
from fspsim.utils.Propagators import sgp4_prop_TLE, numerical_prop

# propagate 3 known TLEs (one LEO, one MEO, one GEO) using SGP4 and RK45 and compare the results in decay rates over 7 days.
def TLEvsNumerical():
    #navstar81, pulled on 27Apr2023
    test_tle1 = "1 48859U 21054A   23116.83449170 -.00000109  00000-0  00000-0 0  9996\n2 48859  55.3054  18.4561 0008790 213.9679 183.6522  2.00556923 13748"
    #OneWeb20, pulled on 27Apr2023
    test_tle2 = "1 45133U 20008C   23116.69886660  .00000678  00000-0  19052-2 0  9998\n2 45133  87.8784 264.1991 0001687  89.2910 270.8411 13.10377378158066"
    #Starlink70, pulled on 27Apr2023
    test_tle3 = "1 53544U 22101T   23122.20221856  .00001510  00000-0  11293-3 0  9999\n2 53544  53.2176  64.0292 0001100  79.8127 280.2989 15.08842383 38928"

    #dictionary of TLEs
    test_tles = {'navstar81': test_tle1, 'OneWeb20': test_tle2, 'Starlink70': test_tle3}
    # for each TLE:
    # prop TLE for 7 days and get the ephemeris. Store in cartesian ECI coordinates
    # get the first ECI value from the ephemeris and use it as the initial state for the numerical propagator
    # prop using numerical prop: Monopole; Monopole + j2; Monopole + j2 + drag; Monopole + j2 + drag + solar radiation pressure; Monopole + j2 + drag + solar radiation pressure + 3BP;
    # plot the RMS over time for each of the numerical propagations against the TLE propagation
    # plot the altitude over time for each of the numerical propagations against the TLE propagation
    for tle in test_tles:
        tle_dict = tle_parse(test_tles[tle])
        tle_kepels = tle_convert(tle_dict)
        #Get the epoch from the TLE
        tle_time_JD = TLE_time(test_tles[tle])
        tle_epoch = jd_to_utc(tle_time_JD)
    
        print("tle epoch:", tle_epoch)
        print("tle time:", tle_time_JD)
        prop_days = 15
        tle_epehem = sgp4_prop_TLE(test_tles[tle], jd_start = tle_time_JD, jd_end = tle_time_JD + prop_days, dt = 20)

        tle_epehem_first_pos = np.array(tle_epehem[0][1])
        tle_epehem_first_vel = np.array(tle_epehem[0][2])
        tle_epehem_first_time = np.array(tle_epehem[0][0])
        # print("tle ephemeris first pos:", tle_epehem_first_pos)
        # print("tle ephemeris first vel:", tle_epehem_first_vel)
        # print("tle ephemeris first time:", tle_epehem_first_time)

        tot_prop_time_in_seconds = prop_days*24*60*60
        integrator_step_size = 20 #seconds
        step_size_in_days = integrator_step_size / (24 * 60 * 60) # convert step_size from seconds to days
        jd_time_stamps = np.arange(tle_time_JD, tle_time_JD + prop_days, step_size_in_days) # create an array of julian dates from jd_start to jd_stop with step size of step_size_in_days

        numerical_ephem = numerical_prop(tot_time = tot_prop_time_in_seconds, h=integrator_step_size, pos = tle_epehem_first_pos, vel= tle_epehem_first_vel,JD_time_stamps = jd_time_stamps, C_d = 2.2, area = 3, mass = 150, integrator_type="RK45")
        print("numerical ephem:", numerical_ephem)
        numerical_alts = ([np.linalg.norm(x[:3])-6378.137 for x in numerical_ephem])
        numerical_jds = jd_time_stamps
        #plot the altitude of each satellite over time
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot([x[0] for x in tle_epehem], [np.linalg.norm(x[1])-6378.137 for x in tle_epehem], label = tle + ' TLE')
        ax.plot(numerical_jds, numerical_alts, label = tle + ' numerical')
        ax.set_xlabel('Time MJD')
        ax.set_ylabel('Altitude (km)')
        ax.set_title('Altitude of Satellites Over Time')
        plt.show()
if __name__ == "__main__":
    TLEvsNumerical()

    ###################

    # def test_numerical_prop_of_TLEs():
#     """
#     Test the implementation of the SGP4 propagator. In particular looking at the decay rates and the effect of the way BSTAR is calculated in the construction of TLEs
#     """

#     #navstar81, pulled on 27Apr2023
#     test_tle1 = "1 48859U 21054A   23116.83449170 -.00000109  00000-0  00000-0 0  9996\n2 48859  55.3054  18.4561 0008790 213.9679 183.6522  2.00556923 13748"
#     #OneWeb20, pulled on 27Apr2023
#     test_tle2 = "1 45133U 20008C   23116.69886660  .00000678  00000-0  19052-2 0  9998\n2 45133  87.8784 264.1991 0001687  89.2910 270.8411 13.10377378158066"
#     #Starlink70, pulled on 27Apr2023
#     test_tle3 = "1 53544U 22101T   23122.20221856  .00001510  00000-0  11293-3 0  9999\n2 53544  53.2176  64.0292 0001100  79.8127 280.2989 15.08842383 38928"

#     #dictionary of TLEs
#     test_tles = {'navstar81': test_tle1, 'OneWeb20': test_tle2, 'Starlink70': test_tle3}

#     prop_start = [datetime.datetime.strptime('2023-05-02 12:45:00', '%Y-%m-%d %H:%M:%S')]
#     prop_end = [datetime.datetime.strptime('2023-05-04 01:48:00', '%Y-%m-%d %H:%M:%S')]
#     start_jd = utc_to_jd(prop_start)
#     end_jd = utc_to_jd(prop_end)

#     for sat in test_tles:
#         print("Propagating numerically: ", sat)
#         sat_altitude = []#list of real and fabricated altitudes for this satellite
#         sat_pos = []#list of real and fabricated positions for this satellite
 
#         ###### SECTION 1: Make TLEs and propagate them ######
#         #Get the Keplerian elements from the TLE
#         tle_dict = tle_parse(test_tles[sat])
#         tle_kepels = tle_convert(tle_dict)
#         #Get the epoch from the TLE
#         tle_time = TLE_time(test_tles[sat])
#         tle_epoch = jd_to_utc(tle_time)
#         # make a SpaceObject from the TLE
#         tle_epoch_str = str(tle_epoch)
#         epoch = tle_epoch_str.replace(' ', 'T')
#         test_sat = SpaceObject(sma = tle_kepels['a'], perigee=tle_kepels['a']-6378.137, apogee=tle_kepels['a']-6378.137, eccentricity=tle_kepels['e'], inc = tle_kepels['i'], argp = tle_kepels['arg_p'], raan=tle_kepels['RAAN'], tran=tle_kepels['true_anomaly'], characteristic_area=0.011, mass = 250, epoch = epoch, launch_date='2023-05-02')
#         test_sat.prop_catobject(jd_start=start_jd[0], jd_stop=end_jd[0], step_size=10, propagator="RK45")
#         test_sat_ephem = test_sat.ephemeris
#         print("test sat ephem:", test_sat_ephem)

#         #3D plots of the position of the satellite over time
#         fig = plt.figure()
#         ax = fig.add_subplot(111, projection='3d')
#         ax.plot([x[0] for x in test_sat_ephem], [x[1] for x in test_sat_ephem], [x[2] for x in test_sat_ephem])
#         ax.set_xlabel('X')
#         ax.set_ylabel('Y')
#         ax.set_zlabel('Z')
#         ax.set_title('Position of Satellite Over Time')
#         plt.show()


# def test_numericalprop_visual():
#     """
#     Test the numerical propagation function by plotting the position of the satellite over time
#     """
#     prop_time = 32400 #seconds
#     step_size = 10 #seconds

#     init_state = (2458849.5, #time
#                 [-7134.4015980671975, -1344.2053503962836, 2616.199171181745,],  #position
#                 [2.7370225687820864, -2.6412753868728953, 6.099437705233797]) # velocity

#     ephemeris = numerical_prop(tot_time = prop_time, h=step_size, pos = init_state[1], vel= init_state[2], cd = 2.2, area = 3, mass = 150)
#     print("begginging of ephemeris:", ephemeris[0:5])
#     print("len of ephemeris:", len(ephemeris))
#     print("shape of ephemeris:", np.shape(ephemeris))
#     #print the norm of the position vector at each time step
#     for i in range(len(ephemeris[0:10])):
#         print(np.linalg.norm(ephemeris[i][0:3]))

#     #plot the position of the satellite over time in 3D
#     fig = plt.figure()
#     ax = fig.add_subplot(111, projection='3d')
#     ax.plot(ephemeris[:,0], ephemeris[:,1], ephemeris[:,2])
#     ax.set_xlabel('X (km)')
#     ax.set_ylabel('Y (km)')
#     ax.set_zlabel('Z (km)')
#     ax.set_title('Position of Satellite Over Time')
#     plt.show()

# def test_rk_vs_sgp4():
#     """
#     Test the implementation of the SGP4 propagator. In particular looking at the decay rates and the effect of the way BSTAR is calculated in the construction of TLEs
#     """

#     prop_start = [datetime.datetime.strptime('2023-05-02 12:45:00', '%Y-%m-%d %H:%M:%S')]
#     prop_end = [datetime.datetime.strptime('2028-04-27 01:48:00', '%Y-%m-%d %H:%M:%S')]
#     start_jd = utc_to_jd(prop_start)
#     end_jd = utc_to_jd(prop_end)
#     t_step = 60*60*24 #1 day in seconds

#     #test propagation of TLEs in test_tles
#     tseries_altitude = [] #for each satellite, this will contain two lists (one for the real and one for the fabricated TLE) of altitude values
#     t_series_pos = [] #for each satellite, this will be a list of position differences between the real and fabricated TLEs

#     for sat in test_tles:
#         print("Propagating TLE for ", sat)
#         sat_altitude = []#list of real and fabricated altitudes for this satellite
#         sat_pos = []#list of real and fabricated positions for this satellite
 
#         ###### SECTION 1: Make TLEs and propagate them ######
#         #Get the Keplerian elements from the TLE
        
#         tle_dict = tle_parse(test_tles[sat])
#         tle_kepels = tle_convert(tle_dict)
#         #Get the epoch from the TLE
#         tle_time = TLE_time(test_tles[sat])
#         tle_epoch = jd_to_utc(tle_time)
#         # make a SpaceObject from the TLE
#         tle_epoch_str = str(tle_epoch)
#         epoch = tle_epoch_str.replace(' ', 'T')
#         test_sat = SpaceObject(sma = tle_kepels['a'], 
#                                perigee=tle_kepels['a']-6378.137, 
#                                apogee=tle_kepels['a']-6378.137, 
#                                eccentricity=tle_kepels['e'], 
#                                inc = tle_kepels['i'], 
#                                argp = tle_kepels['arg_p'], 
#                                raan=tle_kepels['RAAN'], 
#                                tran=tle_kepels['true_anomaly'], 
#                                characteristic_area=4, 
#                                mass = 250, 
#                                epoch = epoch, 
#                                launch_date='2023-05-02')
        
#         test_sat.prop_catobject(start_jd[0], end_jd[0], t_step)
#         test_sat_ephem = test_sat.ephemeris
#         print("made up BStar:", test_sat.bstar)
#         print("made up TLE:", test_sat.tle)
#         print("real TLE:", test_tles[sat])

#         test_pos = [x[1] for x in test_sat_ephem]
#         test_pos = np.array(test_pos)
#         test_times = [x[0] for x in test_sat_ephem]
#         test_altitudes = [np.linalg.norm(x)-6378.137 for x in test_pos]
#         sat_altitude.append(test_altitudes)
#         sat_pos.append(test_pos)

#         # ###### SECTION 2: Use existing TLEs and propagate them ######
#         # valid_tle = test_tles[sat]
#         # valid_tle_ephem = sgp4_prop_TLE(valid_tle, jd_start=start_jd[0], jd_end=end_jd[0], dt=t_step)
#         # valid_tle_pos = [x[1] for x in valid_tle_ephem]
#         # valid_tle_pos = np.array(valid_tle_pos)
#         # valid_tle_times = [x[0] for x in valid_tle_ephem]
#         # valid_tle_altitudes = [np.linalg.norm(x)-6378.137 for x in valid_tle_pos]
#         # sat_altitude.append(valid_tle_altitudes)
#         # sat_pos.append(valid_tle_pos)

#         tseries_altitude.append(sat_altitude)
#         t_series_pos.append(sat_pos)
        
# #plot the altitude of each satellite over time
#     for i in range (len(tseries_altitude)): 
#         fig = plt.figure()
#         ax = fig.add_subplot(111)
#         ax.plot(tseries_altitude[i][0], label = list(test_tles.keys())[i] + 'fabricated')
#         ax.plot(tseries_altitude[i][1], label = list(test_tles.keys())[i] + 'real')
#         ax.set_xlabel('Time MJD')
#         ax.set_ylabel('Altitude (km)')
#         ax.set_title('Altitude of Satellites Over Time')
#         ax.legend()
#         plt.show()
