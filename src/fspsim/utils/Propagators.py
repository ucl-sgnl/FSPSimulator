import numpy as np
import warnings
from sgp4.api import Satrec
from sgp4.api import SGP4_ERRORS
import math
import ssl
ssl._create_default_https_context = ssl._create_unverified_context
from pyatmos import coesa76

# Useful constants
Re = 6378.137  # km
GM_earth = 398600.4418 #km^3/s^2

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
        warnings.warn("jd_start is greater than jd_end. SGP4 Propagation will not be performed.")
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
        error, position, velocity = satellite.sgp4(time, fr)
        if error != 0:
            print("SGP4 ERROR: " + str(error))
            print("Error Message: " + SGP4_ERRORS[error])
            warnings.warn(SGP4_ERRORS[error])
            return ephemeris
        else:
            ephemeris.append([time,position, velocity]) #jd time, pos, vel
            time += dt_jd

    return ephemeris

def kepler_prop(jd_start, jd_stop, step_size, a, e, i, w, W, V, area=None, mass=None, cd=None, drag_decay=False):
    """
    Propagates the orbit of a satellite in a Keplerian orbit, taking drag decay into account.
    Uses Kepler's equation for propagation and Gaussian vectors for coordinate transformation.
    Stops the propagation if altitude drops below 200 km.
    Area, mass, and drag coefficient are required only for drag decay.
    """

    # Constants
    r_tol = 1e-7
    GM_earth = 398600.4415
    i, w, W, V = map(np.deg2rad, [i, w, W, V])
    n = np.sqrt(GM_earth / (a**3))
    
    p = a * (1 - (e**2))
    r = p / (1 + e * np.cos(V))
    cos_Eo = ((r * np.cos(V)) / a) + e
    sin_Eo = (r * np.sin(V)) / (a * np.sqrt(1 - e**2))
    Eo = np.arctan2(sin_Eo, cos_Eo) % (2 * np.pi)
    Mo = Eo - e * np.sin(Eo)
    
    # Gaussian vectors
    P, Q = _compute_gaussian_vectors(W, w, i)
    
    t_diff_secs = 86400 * (jd_stop - jd_start)
    steps = math.ceil(t_diff_secs / step_size)
    
    ephemeris = np.empty((steps, 3, 3))

    for step in range(steps):
        # Compute new eccentric anomaly
        E = _compute_eccentric_anomaly(Mo, e, n, step, step_size, r_tol)
        
        # Calculate position and velocity
        x_new, y_new = a * (np.cos(E) - e), a * np.sqrt(1 - e**2) * np.sin(E)
        cart_pos = x_new * P + y_new * Q
        cart_vel = _compute_velocity(GM_earth, a, e, E, P, Q)

        ephemeris[step, 0] = jd_start + step * step_size / 86400.0
        ephemeris[step, 1] = cart_pos.flatten() #use flatten to convert the column vector to a row vector
        ephemeris[step, 2] = cart_vel.flatten()

        # Check altitude and update semi-major axis
        altitude = np.linalg.norm(cart_pos) - 6378.137
        if altitude <= 200:
            # Reduce the pre-allocated space to the actual data size before breaking
            ephemeris = ephemeris[:step+1]
            break
        
        if drag_decay:
            # Semi-major axis decay due to drag
            a += _sma_drag_decay(cart_pos, cart_vel, cd, area, mass, step_size)

    return ephemeris

def _compute_gaussian_vectors(W, w, i):
    P = np.array([
        [np.cos(W) * np.cos(w) - np.sin(W) * np.cos(i) * np.sin(w)], 
        [np.sin(W) * np.cos(w) + np.cos(W) * np.cos(i) * np.sin(w)], 
        [np.sin(i) * np.sin(w)]
    ])
    
    Q = np.array([
        [-np.cos(W) * np.sin(w) - np.sin(W) * np.cos(i) * np.cos(w)], 
        [-np.sin(W) * np.sin(w) + np.cos(W) * np.cos(i) * np.cos(w)], 
        [np.sin(i) * np.cos(w)]
    ])
    return P, Q

def _compute_eccentric_anomaly(Mo, e, n, step, step_size, r_tol):
    Mi = Mo + n * step * step_size
    E = Mi + (e/2) if Mi < np.pi else Mi - (e/2)
    f = E - e * np.sin(E) - Mi
    f_prime = 1 - e * np.cos(E)
    while abs(f/f_prime) > r_tol:
        E -= f/f_prime
        f = E - e * np.sin(E) - Mi
        f_prime = 1 - e * np.cos(E)
    return E

def _compute_velocity(GM_earth, a, e, E, P, Q):
    f_new = np.sqrt(a * GM_earth) / (a * (1 - e * np.cos(E)))
    g_new = np.sqrt(1 - e**2)
    cos_Ei, sin_Ei = np.cos(E) + e, np.sin(E) * np.sqrt(1 - e**2)
    return (-f_new * sin_Ei * P) + (f_new * g_new * cos_Ei * Q)

def _sma_drag_decay(cart_pos, cart_vel, cd, area, mass, step_size):
    #doing this in meters as rho is in kg/m^3 
    GM_earth_m = 398600.4415e9  # in m^3/s^2
    
    r = np.linalg.norm(cart_pos) * 1000  # Convert to meters
    altitude = r - 6378.137e3  # Convert Earth radius to meters
    altitude_km = altitude / 1000

    MIN_ALTITUDE = 0  # Set to the minimum valid altitude for the model
    MAX_ALTITUDE = 2000  # Set to the maximum valid altitude for the model

    # Check if altitude is within bounds
    if MIN_ALTITUDE <= altitude_km <= MAX_ALTITUDE:
        rho = coesa76(altitude_km).rho
    else:
        return 0 # Don't decay if altitude is out of bounds

    rho = coesa76(altitude_km).rho  # Atmospheric density
    
    vel_ms = np.linalg.norm(cart_vel) * 1000  # Convert to meters/second
    drag_acc = -(cd * area * rho * vel_ms**2) / (2 * mass)  # Drag acceleration in m/s^2
    
    delta_v = drag_acc * step_size
    
    # Semi-major axis decay by vis-viva equation 
    a_initial = 1 / ((2 / r) - (vel_ms**2 / GM_earth_m)) # Initial semi-major axis in meters
    a_final = 1 / ((2 / r) - ((vel_ms + delta_v)**2 / GM_earth_m)) # Final semi-major axis in meters

    delta_a_m = a_final - a_initial
    delta_a_km = delta_a_m / 1000  # Convert to kilometers
    return delta_a_km

if __name__ == "__main__":
    pass