import numpy as np
import warnings
import math
import datetime
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


def tle_parse(tle_string):

    """
    Parses a TLE string (e.g. as provided by SpaceTrack.org) and returns all the data in a dictionary.
    Args:
        tle_string (string): TLE string to be parsed
    Returns:
        tle_dict (dict): dictionary of all the data contained in the TLE string
    """

    # This function takes a TLE string and returns a dictionary of the TLE data
    tle_lines = tle_string.split('\n')
    tle_dict = {}
    line_one, line_two, line_three = tle_lines[0],tle_lines[1], tle_lines[2]
    
    #Parse the first line
    tle_dict['line number'] = line_one[0]
    tle_dict['satellite catalog number'] = line_one[2:7]
    tle_dict['classification'] = line_one[7]
    tle_dict['International Designator(launch year)'] = line_one[9:11] 
    tle_dict['International Designator (launch num)'] = line_one[11:14]
    tle_dict['International Designator (piece of launch)'] = line_one[14:17]
    tle_dict['epoch year'] = line_one[18:20]
    tle_dict['epoch day'] = line_one[20:32]
    tle_dict['first time derivative of mean motion(ballisitc coefficient)'] = line_one[33:43]
    tle_dict['second time derivative of mean motion(delta-dot)'] = line_one[44:52]
    tle_dict['bstar drag term'] = line_one[53:61]
    tle_dict['ephemeris type'] = line_one[62]
    tle_dict['element number'] = line_one[63:68]
    tle_dict['checksum'] = line_one[68:69]

    #Parse the second line (ignore the line number, satellite catalog number, and checksum)
    tle_dict['inclination'] = line_two[8:16]
    tle_dict['right ascension of the ascending node'] = line_two[17:25]
    tle_dict['eccentricity'] = line_two[26:33]
    tle_dict['argument of perigee'] = line_two[34:42]
    tle_dict['mean anomaly'] = line_two[43:51]
    tle_dict['mean motion'] = line_two[52:63]
    tle_dict['revolution number at epoch'] = line_two[63:68]

    return tle_dict

def tle_convert(tle_dict, display=False):
    """
    Converts a TLE dictionary into the corresponding keplerian elements
    
    Args:
        tle_dict (dict): dictionary of TLE data as provided by the tle_parse function

    Returns:
        keplerian_dict(dict): dictionary containing Keplerian elements
    """

    # Standard gravitational parameter for the Earth
    GM = 398600.4415 * (1e3)**3 # m^3/s^2

    # Convert RAAN from degrees to radians
    RAAN = np.radians(float(tle_dict['right ascension of the ascending node']))
    
    # Convert argument of perigee from degrees to radians
    arg_p = np.radians(float(tle_dict['argument of perigee']))
    
    # Convert mean motion from revolutions per day to radians per second
    mean_motion = float(tle_dict['mean motion']) * (2 * np.pi / 86400)
    
    # Compute the period of the orbit in seconds
    period = 2 * np.pi / mean_motion
    
    # Compute the semi-major axis
    n = mean_motion # mean motion in radians per second
    a = (GM / (n ** 2)) ** (1/3) / 1000 # in km
    
    # Convert mean anomaly from degrees to radians
    M = np.radians(float(tle_dict['mean anomaly']))
    
    # Extract eccentricity as decimal value
    e = float("0." + tle_dict['eccentricity'])
    
    # Convert inclination from degrees to radians
    inclination = np.radians(float(tle_dict['inclination']))
    
    # Initial Guess at Eccentric Anomaly
    if M < np.pi:
        E = M + (e / 2)
    else:
        E = M - (e / 2)

    # Numerical iteration for Eccentric Anomaly
    f = lambda E: E - e * np.sin(E) - M
    fp = lambda E: 1 - e * np.cos(E)
    E = np.float64(E)
    r_tol = 1e-8 # set the convergence tolerance for the iteration
    max_iter = 50 # set the maximum number of iterations allowed
    for it in range(max_iter):
        f_value = f(E)
        fp_value = fp(E)
        E_new = E - f_value / fp_value
        if np.abs(E_new - E) < r_tol:
            E = E_new
            break
        E = E_new
    else:
        raise ValueError("Eccentric anomaly did not converge")
        
    eccentric_anomaly = E

    # Compute True Anomaly
    true_anomaly = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(eccentric_anomaly / 2),
                                  np.sqrt(1 - e) * np.cos(eccentric_anomaly / 2))

    # Dictionary of Keplerian elements
    keplerian_dict = {'a': a, 'e': e, 'i': inclination, 'RAAN': RAAN, 'arg_p': arg_p, 'true_anomaly': np.degrees(true_anomaly)}
    if display == True:
        print("Keplerian Elements:")
        print("a = {:.2f} km".format(keplerian_dict['a']))
        print("e = {:.2f}".format(keplerian_dict['e']))
        print("i = {:.2f} deg".format(np.degrees(keplerian_dict['i'])))
        print("RAAN = {:.2f} deg".format(np.degrees(keplerian_dict['RAAN'])))
        print("arg_p = {:.2f} deg".format(np.degrees(keplerian_dict['arg_p'])))
        print("true_anomaly = {:.2f} deg".format(np.degrees(keplerian_dict['true_anomaly'])))

    return keplerian_dict

def TLE_time(TLE):
    """Find the time of a TLE in julian day format"""
    #find the epoch section of the TLE
    epoch = TLE[18:32]
    #convert the first two digits of the epoch to the year
    year = 2000+int(epoch[0:2])
    
    # the rest of the digits are the day of the year and fractional portion of the day
    day = float(epoch[2:])
    #convert the day of the year to a day, month, year format
    date = datetime.datetime(year, 1, 1) + datetime.timedelta(day - 1)
    #convert the date to a julian date
    jd = (date - datetime.datetime(1858, 11, 17)).total_seconds() / 86400.0 + 2400000.5
    return jd