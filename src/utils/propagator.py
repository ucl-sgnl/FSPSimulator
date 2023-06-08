import numpy as np
import warnings
import datetime
from pyatmos import coesa76

#Local imports
from src.utils.coords import UTC_step, utc_to_jd, v_rel

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

def rk4_step(f, t, y, h):
    """
        Calculate one Runge Kutta 4th order step
    Args:
        f(function): function to integrate
        t (float/int): t0
        y(float/int): initial state
        h(float/int): step size to be taken
    Returns:
         (float/int): result of integration step
    """

    k1 = f(t, y)
    k2 = f(t + 0.5 * h, y + 0.5 * k1 * h)
    k3 = f(t + 0.5 * h, y + 0.5 * k2 * h)
    k4 = f(t + h, y + k3 * h)

    return y + h / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)

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

def propagation(self, tot_time, h, x0, y0, z0, u0, v0, w0, model,atmosphere, gravity, plot=True):
    """
    Propagates the entire Orbit object

    Args:
        tot_time (float): total propagation time (seconds)
        h (float): time step of the propagation (seconds)
        x0 (float): cartesian x-component of position (Km)
        y0 (_type_): cartesian y-component of position (Km)
        z0 (_type_): cartesian z-component of position (Km)
        u0 (_type_): cartesian x-component of velocity (Km/s)
        v0 (_type_): cartesian y-component of velocity (Km/s)
        w0 (_type_): cartesian z-component of velocity (Km/s)
        model (string): attitude model type to be used in the calculation of aerodynamic drag acceleration 
        atmosphere (string): atmosphere model to be used in the calcualtion of aerodynamic drag acceleration 
        gravity (string): gravity model type to be used in the calculation of gravitational acceleration
        plot (bool, optional): Whether the propagation orbit should be plotted in 3D. Defaults to True.

    Returns:
        array: nested array containing the cartesian state vectors for the propagated orbit at each time step.
    """
    if h > 200:
        warnings.warn(f'The time step of {h} seconds is large. The results may be inaccurate.')
    
    steps = int(tot_time/h)
    ets = np.zeros((steps,1))
    states = np.zeros((steps,6))
    states[0] = [x0, y0, z0, u0, v0, w0]
    
    self.UTC_times = UTC_step(self.UTC_TIME, steps=steps, h=h) # List of UTC times of each step
    self.jd_times = utc_to_jd(self.UTC_times) # Convert the UTC time to JD
    self.time_list  = np.linspace(0, tot_time, steps) # List of time steps in seconds from 0 to tot_time

    for step in range(steps - 1):
        step = step #this is where the step counter is updated
        time += h #updating the time counter accordingly
        r_vec = np.array([states[step,0], states[step,1], states[step,2]]) #position vector
        alt = np.linalg.norm(r_vec) - 6378.137 #altitude
        if alt < 105:
            break
        else:
            if self.prop_method == "rk4":
                states[step + 1] = rk4_step(self.accelerations, ets[step], states[step], h)

    return states

def accelerations (self,state):
    """For a single time step, calculates the accelerations for a given instance of an Orbit object

    Args:
        t (float): time stamp (seconds)
        state (array): cartesian state vector corresponding to the time stamp. Must contain cartesian position and velocity [x,y,z,u,v,w]

    Returns:
        array: cartesian velocity and the new accelerations in the x,y,z dimensions
    """

    r = state[:3] #extract the r vector from the state vector
    v = state[3:6] #extract the v vector from the state vector
    r_norm = np.linalg.norm(r) #norm of r vector

    #--------------- MONOPOLE ACCELERATION -----------------#
    grav_mono=grav_acc(state) #calculate the acceleration due to gravity
    
    #--------------- J2 PERT ACCELERATION ------------------#
    a_j2 = j2_acc(state) #calculate the acceleration due to J2 perturbations

    #--------------- TOTAL GRAVITATIONAL ACCELERATION ------#
    grav_a = grav_mono + a_j2
    
    self.grav_accs.append(grav_a) #append the gravitational acceleration to the list of accelerations

    #--------------- AERO DRAG ACCELERATION --------------#

    alt = r_norm - Re #altitude is the norm of the r vector minus the radius of the earth
    rho = float(ussa76_rho(alt)) #get the density of the air at the altitude from USSA76 model

    # Then we apply the density to the drag force model

    v_rel_atm= v_rel(state)
    v_norm = np.linalg.norm(v_rel_atm) #norm of v vector
    alt = r_norm - Re #altitude is the norm of the r vector minus the radius of the earth
    acc_drag = -0.5 * rho * ((self.cd*self.area)/self.mass)*(1000*(v_rel_atm**2)) # Vallado and Finkleman 2014
    drag_a_vec = acc_drag * (v_rel_atm / v_norm) # Multiply by unit direction vector to apply drag acceleration
    
    #TODO: Add solar radiation pressure acceleration

    #--------------- SUM OF ALL ACCELERATIONS --------------#
    
    a_tot = grav_a + drag_a_vec
    return np.array([state[3],state[4],state[5],a_tot[0],a_tot[1],a_tot[2]])

def propagate(self, prop_time, step_size , model, atmosphere, gravity, plot=True):
    """ Propagates the orbit for a given time period with a given step size, applying the chosen attitude, gravity and atmosphere model.

    Args:
        prop_time (float): time to propagate the orbit for (seconds)
        step_size (float): time step size (seconds)
        model (str): attitude model type. One of '1D' or '3D_elastic'
        atmosphere (str): atmosphere model type. One of 'JB08' or 'USSA76' or None
        gravity (str): gravity model type. One of 'Monopole' or 'Monopole+J2'
        plot (bool): whether to plot the results or not. Default is True.

    Returns:
        tuple: (ephemeris, density, altitude, drag vector, julian day time stamp)

    Example:
        >>> Orbit.propagate(prop_time=86400, step_size=60, model='3D_elastic', atmosphere='JB08', gravity='Monopole+J2')
    """

    self.ephemeris = self.propagation(tot_time = prop_time, h = step_size, x0 = self.latest_state[0], y0 = self.latest_state[1], z0 = self.latest_state[2], u0 = self.latest_state[3], v0 = self.latest_state[4], w0 = self.latest_state[5], plot = plot, model = model, atmosphere = atmosphere, gravity = gravity)
    
    rhos_avg = []
    alts_avg = []
    drags_avg = []
    # jd_times_avg = []

    #take the average for every 4 values in the rhos list to get a value for every step. Then slice the list to get the same length as the ephemeris
    for i in range(0,len(self.rhos),4):
        rhos_avg.append(np.mean(self.rhos[i:i+4]))
        alts_avg.append(np.mean(self.alts[i:i+4]))
        drags_avg.append(np.mean(self.drags[i:i+4]))
        # jd_times_avg.append(np.mean(self.jd_times[i:i+4]))
    
    return self.ephemeris, rhos_avg, alts_avg, drags_avg, self.jd_times