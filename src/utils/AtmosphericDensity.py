import numpy as np
#all of this is basically yeeted from pyatmos source code.
# I have just stripped out the parts I don't need to make it faster.

def lapse_tp(t_lower, p_lower, lr, h_lower, h_upper):
    '''
    Calculate the temperature and pressure at a given geopotential altitude above base of a specific layer.
    The temperature is computed by the linear interpolation with the slope defined by lapse rates. 
    The pressure is computed from the hydrostatic equations and the perfect gas law. 
    See the detailed documentation in http://www.pdas.com/hydro.pdf

    Usage:
    [t_upper, p_upper] = lapse_tp(t_lower, p_lower, lr, h_lower, h_upper)

    Inputs:
    t_lower -> [float] temperature[K] at the lower boundary of the subset in the specific layer
    p_lower -> [float] pressure[Pa] at the lower boundary of the subset in the specific layer
    lr -> [float] lapse rate[K/m] for the specific layer
    h_lower -> [float] geopotential altitude[m] at the lower boundary of the subset in the specific layer
    h_upper -> [float] geopotential altitude[m] at the upper boundary of the subset in the specific layer
    
    Outputs:
    t1 -> [float] temperature[K] at the upper boundary of the subset in the specific layer
    p1 -> [float] pressure[Pa] at the upper boundary of the subset in the specific layer
    
    Reference: Public Domain Aeronautical Software(http://www.pdas.com/atmos.html) 
               https://gist.github.com/buzzerrookie/5b6438c603eabf13d07e
    '''
    M0 = 28.9644 # mean molecular weight, [kg/kmol]
    R_air = 8314.32/M0
    g0 = 9.80665 # standard gravity acceleration, [m/s^2]

    if lr == 0:
        t_upper = t_lower
        p_upper = p_lower * np.exp(-g0 / R_air / t_lower * (h_upper - h_lower)*1e3)
    else:
        t_upper = t_lower + lr * (h_upper - h_lower)
        p_upper = p_lower * (t_upper / t_lower) ** (-g0 / (lr/1e3) / R_air)

    return t_upper, p_upper

def ussa76(h):
    '''
    Implements the U.S. Standard Atmosphere 1976(USSA76) up to 86km. 
    The standard atmosphere is defined as a set of layers by specified geopotential altitudes and lapse rates.
    The temperature is computed by linear interpolation with the slope defined by the lapse rate. 
    The pressure is computed from the hydrostatic equations and the perfect gas law; the density follows from the perfect gas law. 

    Usage:
    [rho, T, P, C, eta, Kc] = ussa76(h)

    Inputs:
    h -> [float] geopotentail altitude, [km]

    Outputs:
    rho -> [float] density at a given altitude, [kg/m^3]
    T -> [float] temperature ..., [K]
    P -> [float] pressure ..., [Pa]
    C -> [float] speed of sound ..., [m/s]
    eta -> [float] dynamic viscosity ..., [kg/m/s]
    Kc -> [float] thermal conductivity ..., [J/(m*s*K)]
    
    Note: the geometric altitude should be in [-0.611,86] km, otherwise the output will be extrapolated for those input altitudes.

    Reference: 
        U.S. Standard Atmosphere, 1976, U.S. Government Printing Office, Washington, D.C. 
        Public Domain Aeronautical Software(http://www.pdas.com/atmos.html) 
        https://gist.github.com/buzzerrookie/5b6438c603eabf13d07e
        https://ww2.mathworks.cn/help/aerotbx/ug/atmosisa.
        http://www.braeunig.us/space/atmmodel.htm#USSA1976
    '''
    T0 = 288.15 # temperature, [K]
    p0 = 101325 # pressure, [Pa]
    h0 = 0 # height, [km] 
    M0 = 28.9644 # mean molecular weight, [kg/kmol]
    R_air = 8314.32/M0 # gas constant for dry air, [J/kg/K]
    gamma = 1.4 # specific heat ratio, namely the ratio of constant-pressure to constant-volume specific heat

    # the lower atmosphere below 86km is separated into seven layers 
    geopotential_alt = [-np.inf, 11, 20, 32, 47, 51, 71, np.inf] # Geopotential altitudes above MSL, [km]

    lr = np.array([-6.5, 0, 1, 2.8, 0, -2.8, -2]) # Lapse rate, [K/km]   

    for i in range(len(lr)):
        if h <= geopotential_alt[i+1]:
            T, P = lapse_tp(T0, p0, lr[i], h0, h)
            break
        else:
            # if altitudes are greater than the first several layers, then it has to integeate these layers first.
            T0, p0 = lapse_tp(T0, p0, lr[i], h0, geopotential_alt[i+1])
            h0 = geopotential_alt[i+1]

    # density
    rho = P / (R_air * T)

    # # speed of sound
    # C = np.sqrt(gamma * R_air * T)

    # # dynamic viscosity by Sutherland's law
    # eta = 1.458e-6*T**1.5/(T+110.4) 

    # # thermal conductivity
    # Kc = 2.64638e-3 * T ** 1.5 / (T + 245.4 * (10 ** (-12.0 / T))) 

    return rho

def ussa76_rho(alts):
    '''
    Implements the U.S. Committee on Extension to the Standard Atmosphere(COESA 1976).
    This version only accepts geometric altitudes (in km).

    Inputs:
    alts -> [float list/array] geometric altitudes, [km]

    Outputs:
    rhos -> [float array] densities at a set of given altitudes, [kg/m^3]
    '''
    # Convert alts to at least 1-D array
    alts = np.atleast_1d(alts)

    # Base altitude for the COESA 1976, [km].
    zb = np.array([86, 91, 100, 110, 120, 150, 200, 300, 500, 750, np.inf])

    # load the coefficients used to approximate density and pressure above 86km 
    data = np.load('src/data/density_coeffs/coesa76_coeffs.npz')
    rho_coeffs, _ = data['rho'], data['p']

    # Test if altitudes are inside valid range
    # check_altitude(alts,(-0.611,1e3),'warning')  

    rhos = np.zeros_like(alts)
    for i, z in enumerate(alts):
        if z <= zb[0]: 
            rho = ussa76(z)
        else:
            ind = np.where((z - zb) >= 0)[0][-1]
            poly_rho = np.poly1d(rho_coeffs[ind])
            rho = np.exp(poly_rho(z))

        rhos[i] = rho

    return rhos