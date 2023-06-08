import numpy as np
import warnings
import datetime
from pyatmos import coesa76
import matplotlib.pyplot as plt 
from scipy.integrate import solve_ivp

#Local imports
from utils.coords import UTC_step, utc_to_jd, v_rel, tle_parse, tle_convert, TLE_time, jd_to_utc
# from SpaceObject import SpaceObject

# Useful constants
Re = 6378.137  # km
GM = 3.9860043543609598e05  # km^3 / s^2
j2 = 1.082626925638815e-3  # km3/s2

def ussa76_rho(altitude):
    """Return an atmospheric density value from the USSA76 model
    Args:
        altitude (float): altitude in kilometres
    Returns:
        float: value of atmospheric density in Kg/m^3
    """
    atmosphere = coesa76(altitude)
    return atmosphere.rho

def grav_acc(state):
    """
    Calculate the acceleration of the satellite in the ECI frame
    Args:
        state(float/int): state of the satellite
    Returns:
         (float/int): result of acceleration
    """
    r = state[:3]
    acc = -GM * r / np.linalg.norm(r) ** 3
    return acc

def j2_acc(state):
    """
    Calculate the acceleration of the satellite in the ECI frame
    Args:
        state(float/int): state of the satellite
    Returns:
         (float/int): result of acceleration
    """
    r = np.linalg.norm(state[0:3])
    r_norm = np.linalg.norm(r)

    z2 = state[2] ** 2
    r2 = r_norm**2
    tx = state[0] / r_norm * (5 * z2 / r2 - 1)
    ty = state[1] / r_norm * (5 * z2 / r2 - 1)
    tz = state[2] / r_norm * (5 * z2 / r2 - 3)
    a_j2 = (
        1.5 * j2 * GM * Re**2 / r2**2 * np.array([tx, ty, tz])
    )  # calculate the acceleration due to J2
    return a_j2

def rk4_step(f, t, y, h, **kwargs):
    """
        Calculate one Runge Kutta 4th order step
    Args:
        f(function): function to integrate
        t (float/int): t0
        y(float/int): initial state
        h(float/int): step size to be taken
        kwargs: additional arguments to pass to function f
    Returns:
         (float/int): result of integration step
    """

    k1 = f(t, y, **kwargs)
    k2 = f(t + 0.5 * h, y + 0.5 * k1 * h, **kwargs)
    k3 = f(t + 0.5 * h, y + 0.5 * k2 * h, **kwargs)
    k4 = f(t + h, y + k3 * h, **kwargs)

    return y + h / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)

def accelerations (t, state, cd, area, mass):
    """For a single time step, calculates the accelerations for a given instance of an Orbit object

    Args:
        t (float): current time
        state (array): cartesian state vector corresponding to the time stamp. Must contain cartesian position and velocity [x,y,z,u,v,w]
        cd (float): drag coefficient
        area (float): cross-sectional area of the satellite (m^2)
        mass (float): mass of the satellite (kg)

    Returns:
        array: cartesian velocity and the new accelerations in the x,y,z dimensions
    """

    r = state[0:3] #extract the r vector from the state vector
    r_norm = np.linalg.norm(r) #norm of r vector
    #--------------- MONOPOLE ACCELERATION -----------------#
    
    grav_mono=grav_acc(state) #calculate the acceleration due to gravity
    #--------------- J2 PERT ACCELERATION ------------------#
   
    a_j2 = j2_acc(state) #calculate the acceleration due to J2 perturbations
    #--------------- TOTAL GRAVITATIONAL ACCELERATION ------#
    
    grav_a = grav_mono + a_j2
    #--------------- AERO DRAG ACCELERATION --------------#

    alt = r_norm - Re #altitude is the norm of the r vector minus the radius of the earth
    rho = float(ussa76_rho(alt)) #get the density of the air at the altitude from USSA76 model

    # Then we apply the density to the drag force model

    v_rel_atm= v_rel(state)
    v_norm = np.linalg.norm(v_rel_atm) #norm of v vector
    alt = r_norm - Re #altitude is the norm of the r vector minus the radius of the earth
    acc_drag = -0.5 * rho * ((cd*area)/mass)*(1000*(v_rel_atm**2)) # Vallado and Finkleman 2014
    drag_a_vec = acc_drag * (v_rel_atm / v_norm) # Multiply by unit direction vector to apply drag acceleration
    
    #TODO: Add solar radiation pressure acceleration

    #--------------- SUM OF ALL ACCELERATIONS --------------#
    
    a_tot = grav_a + drag_a_vec
    return np.array([state[3],state[4],state[5],a_tot[0],a_tot[1],a_tot[2]])

def stop_propagation(t, y, *args):
    """Event function for solve_ivp to stop integration when altitude is less than 105 km."""
    r = y[:3]
    alt = np.linalg.norm(r) - 6378.137
    return alt - 105  # When this is 0, altitude is 105 km

stop_propagation.terminal = True  # Stop the integration when this event occurs

def numerical_prop(tot_time, pos, vel, C_d, area, mass, h=10, type = "RK45"):
    """
    Numerical Propagation of the orbit

    Args:
        tot_time (float): total propagation time (seconds)
        h (float): time step of the propagation (seconds). Defaults to 10 seconds.
        pos (array): cartesian position vector [x,y,z] (km) at initial conditions
        vel (array): cartesian velocity vector [u,v,w] (km/s) at initial conditions
        C_d (float): drag coefficient
        area (float): cross-sectional area of the satellite (m^2)
        mass (float): mass of the satellite (kg)
        type (str): type of numerical integration to use. Defaults to "RK45".

    Returns:
        array: nested array containing the cartesian state vectors for the propagated orbit at each time step.
    """

    if h > 200:
        warnings.warn(f'The time step of {h} seconds is large. The results may be inaccurate.')

    pos = np.array(pos) #cast to numpy array
    vel = np.array(vel)
    
    x0 = np.concatenate((pos, vel))  # Initial state is position and velocity

    # Call solve_ivp to propagate the orbit
    sol = solve_ivp(accelerations, [0, tot_time], x0, method=type, t_eval=np.arange(0, tot_time, h),
                    args=(C_d, area, mass), events=stop_propagation, rtol=1e-8, atol=1e-8)
    return sol.y.T  # Returns an array where each row is the state at a time

# if __name__ == "__main__":

#     def test_numericalprop_visual():
#         """
#         Test the numerical propagation function by plotting the position of the satellite over time
#         """
#         prop_time = 32400 #seconds
#         step_size = 10 #seconds

#         init_state = (2458849.5, #time
#                     [-7134.4015980671975, -1344.2053503962836, 2616.199171181745,],  #position
#                     [2.7370225687820864, -2.6412753868728953, 6.099437705233797]) # velocity

#         ephemeris = numerical_prop(tot_time = prop_time, h=step_size, pos = init_state[1], vel= init_state[2], cd = 2.2, area = 3, mass = 150)
#         print("begginging of ephemeris:", ephemeris[0:5])
#         print("len of ephemeris:", len(ephemeris))
#         print("shape of ephemeris:", np.shape(ephemeris))
#         #print the norm of the position vector at each time step
#         for i in range(len(ephemeris[0:10])):
#             print(np.linalg.norm(ephemeris[i][0:3]))

#         #plot the position of the satellite over time in 3D
#         fig = plt.figure()
#         ax = fig.add_subplot(111, projection='3d')
#         ax.plot(ephemeris[:,0], ephemeris[:,1], ephemeris[:,2])
#         ax.set_xlabel('X (km)')
#         ax.set_ylabel('Y (km)')
#         ax.set_zlabel('Z (km)')
#         ax.set_title('Position of Satellite Over Time')
#         plt.show()

# test_numericalprop_visual()

# Ultimately I need to be able to go:
# for satellite in SATCAT.Catalogue:
#         satellite.prop_catobject(jd_start, jd_stop, timestep)

# def test_rk_vs_sgp4():
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

#     pass