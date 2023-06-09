import numpy as np
import warnings
import math
import datetime
import glob
from astropy import units as u
from astropy.time import Time
from poliastro.bodies import Earth
from poliastro.twobody import Orbit
from poliastro.frames import Planes
from jplephem.spk import SPK

Re = 6378.137 #km Earth's equatorial radius

def utc_to_jd(time_stamps):
    """Converts UTC time to Julian Date.
    Args: 
        time_stamps(list): list of UTC times datetime.datetime objects that need to be converted to JD.
    Returns:
        jd(list): is a list of Julian Date times.
    """
    UTC_string = []
    for i in range(0,len(time_stamps),1):
        try:
            UTC_string.append(time_stamps[i].strftime('%Y-%m-%d %H:%M:%S'))
        except:
            time_str = datetime.datetime.strptime(time_stamps, '%Y-%m-%d')
            UTC_string.append(time_str.strftime('%Y-%m-%d %H:%M:%S'))

    t = Time(UTC_string, format='iso', scale='utc') #astropy time object
    jd = t.to_value('jd', 'long') #convert to jd

    jd_vals = []
    for i in range (0, len(jd),1):
        jd_vals.append(float(jd[i]))
    
    return jd_vals

def jd_to_utc(jd):
    """Converts Julian Date to UTC time tag(datetime object) using Astropy"""
    #convert jd to astropy time object
    time = Time(jd, format='jd', scale='utc')
    #convert astropy time object to datetime object
    utc = time.datetime
    return utc

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

def car2kep(x, y, z, u, v, w, deg=False, arg_l=False):
    """Convert cartesian to keplerian elements.

    Args:
        x (float): x position in km
        y (float): y position in km
        z (float): z position in km
        u (float): x velocity in km/s
        v (float): y velocity in km/s
        w (float): z velocity in km/s
        deg (bool, optional): If True, return angles in degrees. If False, return angles in radians. Defaults to False.
        arg_l (bool, optional): If True, return argument of latitude in degrees. If False, return argument of latitude in radians. Defaults to False.

    Returns:
        tuple: a, e, i, w, W, V, arg_lat
    """
    #TODO: add argument of latitude
    #make the vectors in as astropy Quantity objects
    r = [x, y, z] * u.km
    v = [u, v, w] * u.km / u.s

    #convert to cartesian
    orb = Orbit.from_vectors(Earth, r, v, plane=Planes.EARTH_EQUATOR)

    #convert to keplerian
    if deg == True:

        a = orb.a.value
        e = orb.ecc.value
        i = np.rad2deg(orb.inc.value)
        w = np.rad2deg(orb.raan.value)
        W = np.rad2deg(orb.argp.value)
        V = np.rad2deg(orb.nu.value)
        # arg_lat = np.rad2deg(orb.arglat.value)

        if arg_l == True:
            return a, e, i, w, W, V
        elif arg_l == False:
            return a, e, i, w, W, V

    elif deg == False:
        a = orb.a.value
        e = orb.ecc.value
        i = orb.inc.value
        w = orb.raan.value
        W = orb.argp.value
        V = orb.nu.value
        # arg_lat = orb.arg_lat.value

        if arg_l == True:
            return a, e, i, w, W, V

    return a, e, i, w, W, V                

def true_to_eccentric_anomaly(true_anomaly, eccentricity):
    cos_E = (eccentricity + math.cos(true_anomaly)) / (1 + eccentricity * math.cos(true_anomaly))
    sin_E = math.sqrt(1 - cos_E**2)
    if true_anomaly > math.pi:
        sin_E = -sin_E
    eccentric_anomaly = math.atan2(sin_E, cos_E)
    return eccentric_anomaly

def calculate_energy_from_semi_major_axis(a, m):
    """Calculate the total mechanical energy of an object in orbit around the Earth, given its semi-major axis (a) and mass (m)."""
    G = 6.67430e-11  # Gravitational constant, units: m^3 kg^-1 s^-2
    M_earth = 5.97219e24  # Mass of Earth, units: kg
    energy = -G * M_earth * m / (2 * a)
    return energy

def eccentric_to_mean_anomaly(eccentric_anomaly, eccentricity):
    mean_anomaly = eccentric_anomaly - eccentricity * math.sin(eccentric_anomaly)
    return mean_anomaly

def true_to_mean_anomaly(true_anomaly, eccentricity):
    #TODO: i dont know where I have lost a minus sign but I made the output here the absolute value to deal wtih it. Needs fixing.
    eccentric_anomaly = true_to_eccentric_anomaly(true_anomaly, eccentricity)
    mean_anomaly = eccentric_to_mean_anomaly(eccentric_anomaly, eccentricity)
    return np.abs(mean_anomaly)

def orbit_classify(altitude):
    # Classifies the orbit based on the altitude
    #TODO: expand to include all possible altitudes. Classifications must be set to match JSR/Celestrak orbit classifications
    # if less than 160 km, it is deorbiting
    if altitude < 160:
        orbit_class = 'Deorbiting'
    # if between 160 and 2000 km, it is low earth orbit
    elif altitude >= 160 and altitude < 2000:
        orbit_class = 'LEO'
    # between 2000 and 10000 km, it is medium earth orbit
    elif altitude >= 2000 and altitude < 10000:
        orbit_class = 'MEO'
    # between 10000 and 36000 km, it is GEO
    elif altitude >= 35000 and altitude < 36000:
        orbit_class = 'GEO'
    else:
        orbit_class = 'Other'
    return orbit_class

def tle_parse(tle_2le):
    #NOTE: does not work on 3LE strings
    """
    Parses a 2LE string (e.g. as provided by Celestrak) and returns all the data in a dictionary.
    Args:
        tle_2le (string): 2LE string to be parsed
    Returns:
        2le_dict (dict): dictionary of all the data contained in the TLE string
    """

    # This function takes a TLE string and returns a dictionary of the TLE data
    tle_lines = tle_2le.split('\n')
    tle_dict = {}
    line_one, line_two = tle_lines[0],tle_lines[1]
    
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
    
    ### Args:
        - tle_dict (dict): dictionary of TLE data as provided by the tle_parse function

    ### Returns:
        - keplerian_dict(dict): dictionary containing Keplerian elements
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
    true_anomaly = (2 * np.arctan2(np.sqrt(1 + e) * np.sin(eccentric_anomaly / 2),
                               np.sqrt(1 - e) * np.cos(eccentric_anomaly / 2))) % (2 * np.pi)

    # Dictionary of Keplerian elements
    keplerian_dict = {'a': a, 'e': e, 'i': inclination, 'RAAN': RAAN, 'arg_p': arg_p, 'true_anomaly': (true_anomaly)}
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
    #take the first tow numerical digits of the epoch (deal with white )
    #convert the first two digits of the epoch to the year
    year = int(2000)+int(epoch[0:2])
    # the rest of the digits are the day of the year and fractional portion of the day
    day = float(epoch[2:])
    #convert the day of the year to a day, month, year format
    date = datetime.datetime(year, 1, 1) + datetime.timedelta(day - 1)
    #convert the date to a julian date
    jd = (date - datetime.datetime(1858, 11, 17)).total_seconds() / 86400.0 + 2400000.5
    return jd

def orbital_period(semi_major_axis):
    G = 6.67430 * 10**(-11)  # Gravitational constant in m^3 kg^-1 s^-2
    a = semi_major_axis*1000      # Semi-major axis in Km
    M = 5.972 * 10**24       # Mass of Earth in kg
    
    # Calculate orbital period using Kepler's Third Law
    orbital_period_seconds = 2 * math.pi * math.sqrt(a**3 / (G * M))
    # convert orbital period to minutes
    orbital_period_minutes = orbital_period_seconds / 60
    
    return orbital_period_minutes

def generate_cospar_id(launch_year, launch_number, launch_piece):
    """
    Generates a COSPAR ID for a spacecraft.

    Args:
    launch_year (int): The launch year of the spacecraft.
    launch_number (int): The launch number of the spacecraft.
    launch_piece (str): The piece of the launch.

    Returns:
    str: The generated COSPAR ID.
    """
    # Format the launch year by taking the last two digits
    year_str = str(launch_year)

    # Format the launch number to have a minimum width of 3 digits, padded with zeros
    launch_number_str = f"{launch_number:03d}"
    # chop the launch number to 3 digits if it is longer than 3 digits
    if len(launch_number_str) > 3:
        launch_number_str = launch_number_str[-3:]

    # Combine the formatted parts to create the COSPAR ID
    cospar_id = f"{year_str}-{launch_number_str}{launch_piece}"

    return cospar_id

def tle_exponent_format(value):
    # Format a value in scientific notation used in TLEs
    if value == 0:
        return "00000-0"

    exponent = int('{:e}'.format(value).split('e')[-1])
    mantissa = '{:.5e}'.format(abs(value)).replace('.', '')[:5]
    sign = '-' if value < 0 else '+'
    exponent_sign = '-' if exponent < 0 else '+'
    exponent_str = '{}{}'.format(exponent_sign, abs(exponent)).rstrip('0')
    formatted_value = "{}{}{}".format(sign, mantissa, exponent_str)
    return formatted_value

def tle_checksum(line):
    checksum = 0
    for char in line[:-1]:
        if char.isdigit():
            checksum += int(char)
        elif char == '-':
            checksum += 1
    return checksum % 10

def write_tle(catalog_number, classification, launch_year, launch_number, launch_piece,
              epoch_year, epoch_day, first_derivative, second_derivative, drag_term,
              ephemeris_type, element_set_number,
              inclination, raan, eccentricity, arg_perigee, mean_anomaly, mean_motion,
              revolution_number):

    formatted_drag_term = tle_exponent_format(drag_term)
    #drop the sign 
    formatted_drag_term = formatted_drag_term[1:]
    formatted_second_derivative = tle_exponent_format(second_derivative)
    first_derivative= "0.00000140" #TODO: unhardcode this
    l1_col1 = ' 1'
    l1_col2 = ' '
    l1_col3_7 = '{:05d}'.format(catalog_number)
    l1_col8 = classification
    l1_col9_10 = ' '
    l1_col10_11 = '{:02d}'.format(launch_year)
    l1_col12_14 = '{:03d}'.format(launch_number)
    l1_col15_17 = '{:3s}'.format(launch_piece)
    l1_col18 = ' '
    l1_col19_20 = '{:02d}'.format(epoch_year)
    l1_col21_32 = '{:012.8f}'.format(epoch_day)
    l1_col33 = ' '
    l1_col34_43 = '{:11.8s}'.format(first_derivative)
    l1_col44 = ' '
    l1_col45_52 = '{:8.8s}'.format(formatted_second_derivative)
    l1_col54_61 = '{:8.8s}'.format(formatted_drag_term)
    l1_col62 = ' '
    l1_col63 = '{:01d}'.format(ephemeris_type)
    l1_col64 = ' '
    l1_col65_68 = '{:04d}'.format(element_set_number)
    l1_upto_68 = l1_col1 + l1_col2 + l1_col3_7 + l1_col8 + l1_col9_10 + l1_col10_11 + l1_col12_14 + l1_col15_17 + l1_col18 + l1_col19_20 + l1_col21_32 + l1_col33 + l1_col34_43 + l1_col44 + l1_col45_52  + l1_col54_61 + l1_col62 + l1_col63 + l1_col64 + l1_col65_68
    l1_col69 = tle_checksum(l1_upto_68)
    l1 = l1_upto_68 + str(l1_col69)

    l2_col1 = '2'
    l2_col2 = ' '
    l2_col3_7 = '{:05d}'.format(catalog_number)
    l2_col8 = ' '
    l2_col9_16 = '{:8.4f}'.format(np.rad2deg(inclination))
    l2_col17 = ' '
    l2_col18_25 = '{:8.4f}'.format(np.rad2deg(raan))
    l2_col26 = ' '
    l2_col27_33 = '{:7.7s}'.format(str(eccentricity)[2:]) #remove the 0. from the beginning
    l2_col34 = ' '
    l2_col35_42 = '{:8.4f}'.format(np.rad2deg(arg_perigee))
    l2_col43 = ' '
    l2_col44_51 = '{:8.4f}'.format(np.rad2deg(mean_anomaly))
    l2_col52 = ' '
    l2_col53_63 = '{:11.8f}'.format(mean_motion)
    l2_col65_68 = '{:5d}'.format(revolution_number)
    l2_upto_68 = l2_col1 + l2_col2 + l2_col3_7 + l2_col8 + l2_col9_16 + l2_col17 + l2_col18_25 + l2_col26 + l2_col27_33 + l2_col34 + l2_col35_42 + l2_col43 + l2_col44_51 + l2_col52 + l2_col53_63  + l2_col65_68
    l2_col69 = tle_checksum(l2_upto_68)
    l2 = l2_upto_68 + str(l2_col69)

    #combine l1 and l2 with a new line character
    newtle = l1 + '\n' + l2

    return newtle

def get_day_of_year_and_fractional_day(epoch):
    start_of_the_year = datetime.datetime(epoch.year, 1, 1)
    day_of_year = (epoch - start_of_the_year).days + 1
    seconds_since_midnight = (epoch - epoch.replace(hour=0, minute=0, second=0, microsecond=0)).total_seconds()
    fractional_day = seconds_since_midnight / 86400
    return day_of_year + fractional_day

def UTC_step(date_list, steps,h):
    """ Make an array of datetime strings in UTC format

    Args:
        UTC (list): UTC time described in the following format: [year,month,day,hour,minute,second]
        steps (int): number of steps in the propagation
        h (float): step size the propagation uses in seconds

    Returns:
        UTC_array (list): list of datetime strings in UTC format
    """
    UTC_times = []
   # convert step size of integrator (seconds) to minutes
    dt = datetime.datetime(date_list[0],date_list[1],date_list[2],date_list[3],date_list[4],date_list[5]) 
    for i in range (0,steps,1): #increment by step size of propagator
        dt = dt + datetime.timedelta(seconds=h)
        UTC_times.append(dt)
    return UTC_times

def v_rel(state):
    """Calcualte the speed of the satellite relative to the Earth's atmosphere assuming it co-rotates with the Earth

    Args:
        state (array-like): state vector of the satellite in ECI coordinates

    Returns:
        v_rel: velocity of the satellite with respect to the Earth's atmosphere in Km/s
    """
    r = np.array([state[0], state[1], state[2]])
    v = np.array([state[3], state[4], state[5]])  # velocity vector
    atm_rot = np.array([0, 0, 72.9211e-6])  # rotation vector of the atmosphere in radians in ECI coordinates
    # speed of satellite relative to atmosphere
    v_rel = v - np.cross(atm_rot, r)
    return v_rel

def earth_sun_vec(jd, unit = True):
    """
    Calculates the Earth-Sun vector in ECI coordinates.
    Args:
        jd (float): Julian Date for which we want the Earth-Sun vector.
        unit (bool): if True, returns the unit vector. If False, returns the vector itself. Default is True.
    Returns:
        earth_sun_vec (array): Earth-Sun vector in ECI coordinates.
    """
    bsp_file = 'src/data/JPL_ephemerides/de421.bsp' #select the first (and only) one
    kernel = SPK.open(bsp_file) # Load the planetary SPK kernel file. This file must be in the working directory.
    position = kernel[0,10].compute(jd) # Solar System Barycenter -> Sun vector for a given mjd
    position -= kernel[0,3].compute(jd) # Solar System Barycenter -> Earth Barycenter vector 
    position -= kernel[3,399].compute(jd) # Earth Barycenter -> Earth vector
    
    if unit == False:
        return position #output in Km
    
    if unit == True:
        return position/np.linalg.norm(position)
    print("length of earth-sun vector is: ", np.linalg.norm(position))

def earth_moon_vec(jd, unit = True):
    """
    Calculates the Earth-Moon vector in ECI coordinates.
    Args:
        jd (float): Julian Date for which we want the Earth-Moon vector.
        unit (bool): if True, returns the unit vector. If False, returns the vector itself. Default is True.
    Returns:
        earth_moon_vec (array): Earth-Moon vector in ECI coordinates.
    """
    bsp_file = 'src/data/JPL_ephemerides/de421.bsp' #select the first (and only) one
    kernel = SPK.open(bsp_file) # Load the planetary SPK kernel file. This file must be in the working directory.
    position = kernel[0,10].compute(jd) # Solar System Barycenter -> Sun vector for a given mjd
    position -= kernel[0,3].compute(jd) # Solar System Barycenter -> Earth Barycenter vector 
    position -= kernel[3,301].compute(jd) # Earth Barycenter -> Moon vector
    
    if unit == False:
        return position #output in Km
    
    if unit == True:
        return position/np.linalg.norm(position)
    print("length of earth-moon vector is: ", np.linalg.norm(position))

def probe_sun_vec(r, jd, unit = False):
    """
    Calculates the probe-Sun vector in ECI coordinates.
    Args:
        r (array): probe position in ECI coordinates.
        jd (float): Julian Date for which we want the probe-Sun vector.
        unit (bool): if True, returns the unit vector. If False, returns the vector itself. Default is False.
    Returns:
        probe_sun_vec (array): Probe-Sun vector in ECI coordinates.
    """
    s = earth_sun_vec(jd, unit = False)
    p = (s-r) #probe-sun vector = earth-sun vector - earth-probe vector
    if unit == False:
        return p #output in Km
    
    if unit == True:
        return p/np.linalg.norm(p)
    print("length of probe-sun vector is: ", np.linalg.norm(p))
    
def probe_moon_vec(r, jd, unit = False):
    """
    Calculates the probe-Moon vector in ECI coordinates.
    Args:
        r (array): probe position in ECI coordinates.
        jd (float): Julian Date for which we want the probe-Moon vector.
        unit (bool): if True, returns the unit vector. If False, returns the vector itself. Default is False.
    Returns:
        probe_moon_vec (array): Probe-Moon vector in ECI coordinates.
    """
    m = earth_moon_vec(jd, unit = False)
    p = (m-r) #probe-moon vector = earth-moon vector - earth-probe vector
    if unit == False:
        return p #output in Km
    
    if unit == True:
        return p/np.linalg.norm(p)
    print("length of probe-moon vector is: ", np.linalg.norm(p))