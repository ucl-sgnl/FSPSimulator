import numpy as np
import warnings
from pyatmos import coesa76
from scipy.integrate import solve_ivp
from sgp4.api import Satrec
import math

#Local imports
from utils.Conversions import v_rel
from utils.AtmosphericDensity import ussa76_rho

# Useful constants
Re = 6378.137  # km
GM = 3.9860043543609598e05  # km^3 / s^2
j2 = 1.082626925638815e-3  # km3/s2

# def ussa76_rho(altitude):
# This one just uses the pyatmos package. Keeping here for reference but have no written a version of this myself to speed up as it is quite slow.
#     """Return an atmospheric density value from the USSA76 model
#     Args:
#         altitude (float): altitude in kilometres
#     Returns:
#         float: value of atmospheric density in Kg/m^3
#     """
#     atmosphere = coesa76(altitude)
#     return atmosphere.rho

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

    #v_rel_atm outputs the relative velocity in km/s. We need to convert it to m/s to use in the drag force model
    v_rel_atm_m = v_rel_atm * 1000
    
    acc_drag = -0.5 * rho * ((cd*area)/mass)*((v_rel_atm_m**2)) # Vallado and Finkleman 2014
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

def numerical_prop(tot_time, pos, vel, C_d, area, mass, h=20, type = "RK45"):
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

    if h > 30:
        warnings.warn(f'The time step of {h} seconds is large. The results may be inaccurate.')

    pos = np.array(pos) #cast to numpy array
    vel = np.array(vel)
    
    x0 = np.concatenate((pos, vel))  # Initial state is position and velocity

    # Call solve_ivp to propagate the orbit
    sol = solve_ivp(accelerations, [0, tot_time], x0, method=type, t_eval=np.arange(0, tot_time, h),
                    args=(C_d, area, mass), events=stop_propagation, rtol=1e-8, atol=1e-8)
    #TODO: I have reduced the tolerance to 1e-4 to get the code to run faster. Need to decide on what is acceptable/necessary here
    return sol.y.T  # Returns an array where each row is the state at a time

def sgp4_prop_TLE(TLE, jd_start, jd_end, dt):

    """Given a TLE, a start time, end time, and time step, propagate the TLE and return the time-series of Cartesian coordinates, and accompanying time-stamps (MJD)
        Note: Simply a wrapper for the SGP4 routine in the sgp4.api package (Brandon Rhodes)
    Args:
        TLE (string): TLE to be propagated
        jd_start (float): start time of propagation in Julian Date format
        jd_end (float): end time of propagation in Julian Date format
        dt (float): time step of propagation in seconds
        alt_series (bool, optional): If True, return the altitude series as well as the position series. Defaults to False.
        
    Returns:
    list: list of lists containing the time-series of Cartesian coordinates, and accompanying time-stamps (MJD)
    
    """

    if jd_start > jd_end:
        print('jd_start must be less than jd_end')
        return

    ephemeris = []
    
    #convert dt from seconds to julian day
    dt_jd = dt/86400

    #split at the new line
    split_tle = TLE.split('\n')

    if len(split_tle) == 3:
        # three line TLE
        s = split_tle[1]
        r = split_tle[2]
    else:
        s = split_tle[0]
        r = split_tle[1]
    fr = 0.0 # precise fraction (SGP4 docs for more info)
    
    #create a satellite object
    satellite = Satrec.twoline2rv(s, r)
    
    time = jd_start
    # for i in range (jd_start, jd_end, dt):
    while time < jd_end:
        # propagate the satellite to the next time step
        # Position is in idiosyncratic True Equator Mean Equinox coordinate frame used by SGP4
        # Velocity is the rate at which the position is changing, expressed in kilometers per second
        error, position, velocity = satellite.sgp4(time, fr)
        if error != 0:
            #print('error: ', error)
            break
        else:
            ephemeris.append([time,position, velocity]) #jd time, pos, vel
            time += dt_jd

    return ephemeris

def kepler_prop(jd_start,jd_stop,step_size,a,e,i,w,W,V):
    """Analytical Propagation. This function calculates the orbit of a
    satellite around the Earth assuming two body-dynamics based
    on Keplerian elements and step size.
    Args:
        jd_start:  starting time in Julian Date
        jd_stop: ending time in Julian Date
        step (int): sampling time of the orbit in seconds
        a (float): semi-major axis in Kilometers
        e (float): eccentricity of the orbit
        i (float): inclination of the orbit in degrees
        w (float): argument of periapsis in degrees
        W (float): longitude of the ascending node in degrees
        V (float): true anomaly in degrees
    Returns:
        ephemeris (list): list in the form [(time, position, velocity), (time, position, velocity), ...]
    """

    # Defining radial distance and semi-latus rectum
    p = a * (1 - (e**2))
    r = p / (1 + e * np.cos(V))
   
    # empty list to store the time step in each iteration
    ephemeris = []
    # calculate the number of steps between jd_start and jd_stop if the step size is step_size seconds) 
    t_diff = jd_stop-jd_start # time difference in Julian Days
    t_diff_secs = 86400 * t_diff
    current_jd = jd_start
    steps = math.ceil(t_diff_secs/step_size) # total number of integer steps
    for i in range(0, steps, 1):
        time = 0 #internal timer
        # Compute the mean motion
        n = np.sqrt(GM / (a**3))
        
        # Compute the eccentric anomaly at t=t0
        cos_Eo = ((r * np.cos(V)) / a) + e
        sin_Eo = (r * np.sin(V)) / (a * np.sqrt(1 - e**2))
        # adding 2 pi for for very small values
        # ensures that Eo stays in the range 0< Eo <2Pi
        Eo = math.atan2(sin_Eo, cos_Eo)
        if Eo < 0.0:
            Eo = Eo + 2 * np.pi
        else:
            Eo = Eo
        # Compute mean anomaly at start point
        Mo = Eo - e * np.sin(Eo)  # From Kepler's equation
        # Compute the mean anomaly at t+newdt
        Mi = Mo + n * time
        # Solve Kepler's equation to compute the eccentric anomaly at t+newdt
        M = Mi
        
        # Initial Guess at Eccentric Anomaly
        # (taken these conditions from
        # Fundamentals of Astrodynamics by Roger.E.Bate)
        r_tol = 1e-7 # relative tolerance
        if M < np.pi:
            E = M + (e / 2)
        if M > np.pi:
            E = M - (e / 2)
        # Initial Conditions
        f = E - e * np.sin(E) - M
        f_prime = 1 - e * np.cos(E)
        ratio = f / f_prime
        # Numerical iteration for ratio compared to level of accuracy wanted
        while abs(ratio) > r_tol:
            f = E - e * np.sin(E) - M
            f_prime = 1 - e * np.cos(E)
            ratio = f / f_prime
            if abs(ratio) > r_tol:
                E = E - ratio
            if abs(ratio) < r_tol:
                break

        Ei = E
        # Compute the gaussian vector component x,y
        x_new = a * (np.cos(Ei) - e)
        y_new = a * ((np.sqrt(1 - e**2)) * (np.sin(Ei)))
        # Compute the in-orbital plane Gaussian Vectors
        # This gives P and Q in ECI components
        P = np.matrix(
            [
                [np.cos(W) * np.cos(w) - np.sin(W) * np.cos(i) * np.sin(w)],
                [np.sin(W) * np.cos(w) + np.cos(W) * np.cos(i) * np.sin(w)],
                [np.sin(i) * np.sin(w)],
            ]
        )
        Q = np.matrix(
            [
                [-np.cos(W) * np.sin(w) - np.sin(W) * np.cos(i) * np.cos(w)],
                [-np.sin(W) * np.sin(w) + np.cos(W) * np.cos(i) * np.cos(w)],
                [np.sin(i) * np.cos(w)],
            ]
        )

        # Compute the position vector at t+newdt
        # We know the inertial vector components along the P and Q vectors.
        # Thus we can project the satellite position onto the ECI basis.
        # calcualting the new x coordiante

        cart_pos_x_new = (x_new * P.item(0)) + (y_new * Q.item(0))
        cart_pos_y_new = (x_new * P.item(1)) + (y_new * Q.item(1))
        cart_pos_z_new = (x_new * P.item(2)) + (y_new * Q.item(2))
        # Compute the range at t+dt
        r_new = a * (1 - e * (np.cos(Ei)))
        # Compute the gaussian velocity components
        cos_Ei = (x_new / a) + e
        sin_Ei = y_new / a * np.sqrt(1 - e**2)
        f_new = (np.sqrt(a * GM)) / r_new
        g_new = np.sqrt(1 - e**2)
       
        cart_vel_x_new = (-f_new * sin_Ei * P.item(0)) + (f_new * g_new * cos_Ei * Q.item(0))  # x component of velocity
        cart_vel_y_new = (-f_new * sin_Ei * P.item(1)) + (f_new * g_new * cos_Ei * Q.item(1))  # y component of velocity
        cart_vel_z_new = (-f_new * sin_Ei * P.item(2)) + (f_new * g_new * cos_Ei * Q.item(2))
        pos = ([cart_pos_x_new, cart_pos_y_new, cart_pos_z_new])
        vel = ([cart_vel_x_new, cart_vel_y_new, cart_vel_z_new])
        
        ephemeris.append([current_jd, pos, vel])

        # Update the JD time stamp
        current_jd = current_jd + step_size
        # Update the internal timer
        time = time + step_size
    return ephemeris

if __name__ == "__main__":
    pass