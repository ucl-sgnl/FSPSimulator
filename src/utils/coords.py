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

