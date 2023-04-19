import numpy as np
import warnings
import math
from astropy import units as u
from astropy.time import Time
from poliastro.bodies import Earth
from poliastro.twobody import Orbit

def kep2car(a, e, i, w, W, V):
    # Suppress the UserWarning for true anomaly wrapping
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        
        # Create an Orbit object from the Keplerian elements
        orbit = Orbit.from_classical(Earth,
                                     a * u.km,
                                     e * u.one,
                                     i * u.rad,
                                     w * u.rad,
                                     W * u.rad,
                                     V * u.rad,
                                     epoch=Time.now())

    # Get the position and velocity vectors in ECI frame
    pos_vec = orbit.r.value
    vel_vec = orbit.v.value

    # Extract the Cartesian coordinates and velocities
    x, y, z = pos_vec
    vx, vy, vz = vel_vec

    return x, y, z, vx, vy, vz

def trueanom2meananom(tran, e):
    """Converts true anomaly to mean anomaly."""
    E = 2 * math.atan(math.sqrt((1 - e) / (1 + e)) * math.tan(tran / 2)) #Eccentric anomaly
    mean_anomaly = E - e * math.sin(E) # Kepler's equation for mean anomaly
    return mean_anomaly

def calculate_kozai_mean_motion(a, mu):
    #TODO: check if mean motion can be used in place of kozai mean motion in the sgp4 propagator
    # I think Kozai mean motion is the mean motion of the osculating orbit, which is not the same as the mean motion of the orbit

    # a: semi-major axis in meters
    # mu: standard gravitational parameter in m^3/s^2 (e.g., Earth's mu = 3.986004418e14)
    
    # Calculate mean motion (no_kozai) in radians/second
    no_kozai_rad_sec = math.sqrt(mu / a**3)

    # Convert mean motion to radians/minute
    no_kozai_rad_min = no_kozai_rad_sec * 60

    return no_kozai_rad_min


def expo_simplified(altitude, alt_type='geometric'):
    #Simple atmospheric density model based on Vallado 2013
    # altitude: altitude in km
    zb = np.array([0., 25, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120,
                   130, 140, 150, 180, 200, 250, 300, 350, 400, 450,
                   500, 600, 700, 800, 900, 1000])

    zb_expand = zb.copy()
    zb_expand[0], zb_expand[-1] = -np.inf, np.inf

    rhob = np.array([1.225, 3.899e-2, 1.774e-2, 3.972e-3, 1.057e-3,
                     3.206e-4, 8.770e-5, 1.905e-5, 3.396e-6, 5.297e-7,
                     9.661e-8, 2.438e-8, 8.484e-9, 3.845e-9, 2.070e-9,
                     5.464e-10, 2.789e-10, 7.248e-11, 2.418e-11,
                     9.518e-12, 3.725e-12, 1.585e-12, 6.967e-13,
                     1.454e-13, 3.614e-14, 1.170e-14, 5.245e-15, 3.019e-15])

    ZS = np.array([7.249, 6.349, 6.682, 7.554, 8.382, 7.714, 6.549,
                   5.799, 5.382, 5.877, 7.263, 9.473, 12.636, 16.149,
                   22.523, 29.740, 37.105, 45.546, 53.628, 53.298,
                   58.515, 60.828, 63.822, 71.835, 88.667, 124.64, 181.05, 268.00])

    if alt_type == 'geopotential':
        altitude = altitude * (1 - 0.0026373 * altitude / (1 + 0.0026373 * altitude))

    ind = np.where((altitude - zb_expand) >= 0)[0][-1]
    rho = rhob[ind] * np.exp(-(altitude - zb[ind]) / ZS[ind])

    return rho


