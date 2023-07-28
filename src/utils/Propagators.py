import numpy as np
import warnings
from scipy.integrate import solve_ivp
from sgp4.api import Satrec
import math

#Local imports
from utils.Conversions import v_rel, probe_sun_vec, probe_moon_vec, earth_sun_vec, earth_moon_vec
from utils.AtmosphericDensity import ussa76_rho, jb08_rho

# Useful constants
Re = 6378.137  # km
GM_earth = 398600.4418 #km^3/s^2
GM_moon = 4902.7779 #km^3/s^2
GM_sun = 132712440041.9394 #km^3/s^2
j2 = 1.082626925638815e-3  # km3/s2
c = 299792.458 #km/s
P_sun = 1367 #W/m^2 (power of the sun at 1 AU)
AU_km = 149597870.7 #km (1 AU in km)
RSun = 695700 #km (radius of the sun)

def monopole_earth_grav_acc(state):
    """
    Calculate the acceleration of the satellite in the ECI frame
    Args:
        state(float/int): state of the satellite
    Returns:
         (float/int): result of acceleration
    """
    r = state[:3] # distance to the center of mass of the Earth
    acc = -GM_earth * r / np.linalg.norm(r) ** 3 # calculate the acceleration due to gravity
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
        1.5 * j2 * GM_earth * Re**2 / r2**2 * np.array([tx, ty, tz])
    )  # calculate the acceleration due to J2
    return a_j2

def monopole_sun_grav_acc(state, jd_time):
    """
    Calculate the acceleration due to the Sun's gravity acting on the satellite in the ECI frame.
    Args:
        state(float/int): State of the satellite.
        jd_time (float): Julian date time.
    Returns:
        (float/int): Resultant acceleration due to the Sun's gravity.
    """
    probe_sun_vector = probe_sun_vec(r= state[:3], jd = jd_time, unit=False) # distance from the probe to the center of mass of the Sun
    earth_sun_vector = earth_sun_vec(jd = jd_time, unit=False) # distance from the Earth to the center of mass of the Sun
    sun_acc = GM_sun * ((probe_sun_vector / np.linalg.norm(probe_sun_vector) ** 3) -  (earth_sun_vector / np.linalg.norm(earth_sun_vector) ** 3))# calculate the acceleration due to gravity
    return sun_acc

def monopole_moon_grav_acc(state, jd_time):
    """
    Calculate the acceleration due to the Moon's gravity acting on the satellite in the ECI frame.
    Args:
        state(float/int): State of the satellite.
        jd_time (float): Julian date time.
    Returns:
        (float/int): Resultant acceleration due to the Moon's gravity.
    """
    probe_moon_vector = probe_moon_vec(r= state[:3], jd = jd_time, unit=False) # distance from the probe to the center of mass of the Moon
    earth_moon_vector = earth_moon_vec(jd = jd_time, unit=False) # distance from the Earth to the center of mass of the Moon
    moon_acc = GM_moon * ((probe_moon_vector / np.linalg.norm(probe_moon_vector) ** 3) -  (earth_moon_vector / np.linalg.norm(earth_moon_vector) ** 3))# calculate the acceleration due to gravity
    return moon_acc

def solar_shadow_function(state, jd_time):
    #TODO: this needs some proper testing to make sure it works. Seems reasonable but have not rigorously tested it.
    """
    Calculate whether the satellite is in the umbra, penumbra, or full phase of the Earth's shadow.
    Args:
        state(numpy.array): State of the satellite.
        jd_time (float): Julian date time.
    Returns:
        (str): "sun", "umbra", or "penumbra"
    """

    probe_sun_vector = probe_sun_vec(state[:3], jd_time, unit=False)
    s = state[:3] - probe_sun_vector # vector from the satellite to the center of mass of the Sun

    a = np.arcsin(RSun/ np.linalg.norm(probe_sun_vector - state[:3])) # apparent radius of the occulted body (the Sun)
    b = np.arcsin(Re / np.linalg.norm(state[:3])) # apparent radius of the occulting body (the Earth)
    c = np.arccos(np.dot(-state[:3], probe_sun_vector) / (np.linalg.norm(state[:3]) * np.linalg.norm(probe_sun_vector))) # Angle between the vectors to the Sun and to the Earth

    if a + b <= c:
        return "sun"  # Satellite is in full sunlight
    elif a + c <= b:
        return "umbra"  # Satellite is in total shadow (umbra)
    else:
        return "penumbra"  # Satellite is in partial shadow (penumbra)

def aero_drag_acc(state, cd, area, mass, density_model, jd):
    """
    Calculate the acceleration due to atmospheric drag acting on the satellite in the ECI frame.
    Args:
        state(float/int): State of the satellite.
        cd(float): Drag coefficient of the satellite.
        area(float): Cross-sectional area of the satellite.
        mass(float): Mass of the satellite.
        density_model(str): The atmospheric density model to use ("ussa76" or "jb08").
        jd(float): Julian date time.
    Returns:
        (float/int): Resultant acceleration due to atmospheric drag.
    """
    possible_density_models = ["ussa76", "jb08"]
    if density_model not in possible_density_models:
        raise ValueError(f"Invalid density model. Must be one of {possible_density_models}")
    #--------------- AERO DRAG ACCELERATION --------------#
    r = state[0:3] #extract the r vector from the state vector
    r_norm = np.linalg.norm(r) #norm of r vector

    alt = r_norm - Re #altitude is the norm of the r vector minus the radius of the earth
    if density_model == "ussa76":
        rho = float(ussa76_rho(alt))
    elif density_model == "jb08":
        rho = float(jb08_rho(alt, jd))
   
    # Then we apply the density to the drag force model
    v_rel_atm= v_rel(state)
    v_norm = np.linalg.norm(v_rel_atm) #norm of v vector

    #v_rel_atm outputs the relative velocity in km/s. We need to convert it to m/s to use in the drag force model
    acc_drag =  -0.5 * rho * ((cd*area)/mass)*(1000*(v_rel_atm**2)) # Vallado and Finkleman (2014)
    drag_a_vec = acc_drag * (v_rel_atm / v_norm) # Multiply by unit direction vector to apply drag acceleration

    return drag_a_vec

def srp_acc(mass, area, state, jd_time, cr=1):
    """
    Calculate the acceleration due to solar radiation pressure (SRP) acting on the satellite in the ECI frame.
    Args:
        mass(float): Mass of the satellite. in kg
        area(float): Cross-sectional area of the satellite. in m^2
        state(float/int): State of the satellite. in km
        jd_time (float): Julian date time. in JD
        cr(float, optional): Coefficient of reflectivity. Defaults to 1 (perfect reflection). 
    Returns:
        (float/int): Resultant acceleration due to SRP.
    """
    probe_sun_vector = probe_sun_vec(state[0:3], jd_time, unit=False) #distance to the center of mass of the Sun
    probe_sun_vec_m = probe_sun_vector * 1000 #convert to metres
    AU_m = AU_km * 1000 #convert AU to metres

    shadow = solar_shadow_function(state, jd_time) #calculate whether the satellite is in the umbra, penumbra, or full phase of the Earth's shadow
    if shadow == "sun":
        Power_Sun = 4.56e-6 #N/m^-2 (power of the sun at 1 AU)
        probe_sun_norm = np.linalg.norm(probe_sun_vec_m)
        a_srp_vec = -Power_Sun*cr*(area/mass)*(probe_sun_vec_m)/(probe_sun_norm**3)*(AU_m**2) #Canonball SRP from Montenbruck and Gill (2000)
        # units
        a_srp_vec = a_srp_vec/1000 #convert from km/s^2 to m/s^2
    elif shadow == "umbra":
        a_srp_vec = np.array([0,0,0])
    elif shadow == "penumbra":
        a_srp_vec = np.array([0,0,0]) #TODO: this is not correct. Need to calculate the SRP acceleration in the penumbra. For now just set to zero.
    return a_srp_vec

def accelerations (t, state, cd, area, mass, jd_time, force_model=["all"]):
    """For a single time step, calculates the accelerations for a given state vector at a given time.

    Args:
        t (float): current time step (seconds)
        state (array): cartesian state vector corresponding to the time stamp. Must contain cartesian position and velocity [x,y,z,u,v,w]
        cd (float): drag coefficient
        area (float): cross-sectional area of the satellite (m^2)
        mass (float): mass of the satellite (kg)
        force_models (list of str): list of force models to be included in the calculation. Possible values include "grav_mono", "a_j2", "drag_aero", "a_srp", "all". Default is ["all"].

    Returns:
        array: cartesian velocity and the new accelerations in the x,y,z dimensions
    """

    a_tot = np.zeros(3)

    #--------------- MONOPOLE ACCELERATION -----------------#
    if "all" in force_model or "grav_mono" in force_model:
        grav_mono=monopole_earth_grav_acc(state)
        a_tot += grav_mono

    #--------------- J2 PERT ACCELERATION ------------------#
    if "all" in force_model or "j2" in force_model:
        a_j2 = j2_acc(state)
        a_tot += a_j2

    #--------------- SUN GRAVITY ------------------------#
    if "all" in force_model or "sun_grav" in force_model:
        sun_grav_mono = monopole_sun_grav_acc(state, jd_time)
        a_tot += sun_grav_mono

    #--------------- MOON GRAVITY -----------------------#
    if "all" in force_model or "moon_grav" in force_model:
        moon_grav_mono = monopole_moon_grav_acc(state, jd_time)
        a_tot += moon_grav_mono

    #--------------- AERO DRAG ACCELERATION --------------#
    if "all" in force_model or "drag_aero" in force_model:
        drag_aero_vec = aero_drag_acc(state, cd, area, mass, density_model="ussa76", jd=jd_time) 
        a_tot += drag_aero_vec

    #--------------- SOLAR RADIATION PRESSURE ACCELERATION --------------#
    if "all" in force_model or "srp" in force_model:
        a_srp_vec = srp_acc(mass, area, state, jd_time, cr=1)
        a_tot += a_srp_vec

    return np.array([state[3],state[4],state[5],a_tot[0],a_tot[1],a_tot[2]])

def stop_propagation(t, y, *args):
    """Event function for solve_ivp to stop integration when altitude is less than 105 km."""
    r = y[:3]
    alt = np.linalg.norm(r) - 6378.137
    return alt - 105  # When this is 0, altitude is 105 km

stop_propagation.terminal = True  # Stop the integration when this event occurs

def numerical_prop(tot_time, pos, vel, C_d, area, mass, JD_time_stamps, h, integrator_type, force_model):
    """
    Numerical Propagation of the orbit

    Args:
        tot_time (float): total propagation time (seconds)
        pos (array): cartesian position vector [x,y,z] (km) at initial conditions
        vel (array): cartesian velocity vector [u,v,w] (km/s) at initial conditions
        C_d (float): drag coefficient
        area (float): cross-sectional area of the satellite (m^2)
        mass (float): mass of the satellite (kg)
        JD_time_stamps (array): array of Julian Date time stamps for each of the time steps
        h (float): time step of the propagation (seconds).
        integrator_type (str): type of numerical integration to use.

    Returns:
        array: nested array containing the cartesian state vectors for the propagated orbit at each time step.
    """

    if h > 30:
        warnings.warn(f'The time step of {h} seconds is large. The results may be inaccurate.')

    pos = np.array(pos) #cast to numpy array
    vel = np.array(vel)
    
    x0 = np.concatenate((pos, vel))  # Initial state is position and velocity

    # Call solve_ivp to propagate the orbit
    sol = solve_ivp(lambda t, state: accelerations(t, state, C_d, area, mass, np.interp(t, np.arange(0, tot_time, h), JD_time_stamps),force_model), [0, tot_time], x0, method=integrator_type, t_eval=np.arange(0, tot_time, h),
                        events=stop_propagation, rtol=1e-13, atol=1e-11)
    #NOTE: these tolerances have been selected to keep the accuracy of the solver within ~1cm after 6 months.
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
    print("split_tle: ", split_tle)

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
        step_size (int): sampling time of the orbit in seconds
        a (float): semi-major axis in Kilometers
        e (float): eccentricity of the orbit
        i (float): inclination of the orbit in degrees
        w (float): argument of periapsis in degrees
        W (float): longitude of the ascending node in degrees
        V (float): true anomaly in degrees
    Returns:
        ephemeris (list): list in the form [(time, position, velocity), (time, position, velocity), ...]
    """
    # Convert angles to radians
    i, w, W, V = map(np.deg2rad, [i, w, W, V])

    # Defining radial distance and semi-latus rectum
    p = a * (1 - (e**2))
    r = p / (1 + e * np.cos(V))

    # empty list to store the time step in each iteration
    ephemeris = []

    # calculate the number of steps between jd_start and jd_stop if the step size is step_size seconds) 
    t_diff = jd_stop-jd_start  # time difference in Julian Days
    t_diff_secs = 86400 * t_diff
    current_jd = jd_start
    steps = math.ceil(t_diff_secs/step_size)  # total number of integer steps
    for step in range(steps):
        # Compute the mean motion
        n = np.sqrt(GM_earth / (a**3))

        # Compute the eccentric anomaly at t=t0
        cos_Eo = ((r * np.cos(V)) / a) + e
        sin_Eo = (r * np.sin(V)) / (a * np.sqrt(1 - e**2))
        Eo = np.arctan2(sin_Eo, cos_Eo)
        Eo = Eo if Eo >= 0 else Eo + 2*np.pi

        # Compute mean anomaly at start point
        Mo = Eo - e * np.sin(Eo)

        # Compute the mean anomaly at t+newdt
        Mi = Mo + n * (step * step_size)

        # Solve Kepler's equation to compute the eccentric anomaly at t+newdt
        M = Mi

        # Initial Guess at Eccentric Anomaly
        if M < np.pi:
            E = M + (e / 2)
        if M >= np.pi:
            E = M - (e / 2)

        # Iterative solution for E
        r_tol = 1e-7  # relative tolerance
        f = E - e * np.sin(E) - M
        f_prime = 1 - e * np.cos(E)
        ratio = f / f_prime
        while abs(ratio) > r_tol:
            f = E - e * np.sin(E) - M
            f_prime = 1 - e * np.cos(E)
            ratio = f / f_prime
            E -= ratio

        Ei = E

        # Compute the gaussian vector component x,y
        x_new = a * (np.cos(Ei) - e)
        y_new = a * (np.sqrt(1 - e**2) * np.sin(Ei))
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
        f_new = (np.sqrt(a * GM_earth)) / r_new
        g_new = np.sqrt(1 - e**2)
       
        cart_vel_x_new = (-f_new * sin_Ei * P.item(0)) + (f_new * g_new * cos_Ei * Q.item(0))  # x component of velocity
        cart_vel_y_new = (-f_new * sin_Ei * P.item(1)) + (f_new * g_new * cos_Ei * Q.item(1))  # y component of velocity
        cart_vel_z_new = (-f_new * sin_Ei * P.item(2)) + (f_new * g_new * cos_Ei * Q.item(2))
        pos = ([cart_pos_x_new, cart_pos_y_new, cart_pos_z_new])
        vel = ([cart_vel_x_new, cart_vel_y_new, cart_vel_z_new])
        
        ephemeris.append([current_jd, pos, vel])
        # Update the JD time stamp
        current_jd = current_jd + step_size/86400.0  # convert seconds to days
    return ephemeris

if __name__ == "__main__":
    pass