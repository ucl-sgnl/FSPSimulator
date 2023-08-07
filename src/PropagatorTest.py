# #Only way i could rin this was to bring it into src and run using "python -m src.tests.PropagatorTest"
import numpy as np
import datetime
import matplotlib.pyplot as plt
import sp3
import astropy.time
import astropy.coordinates
import pandas as pd
import astropy.units as u
from utils.Conversions import car2kep
# Convert to GCRS 
from astropy import units as u
from astropy.coordinates import CartesianDifferential, CartesianRepresentation, GCRS, ITRS
from astropy.time import Time
import unittest
from utils.Conversions import tle_parse, tle_convert, TLE_time, jd_to_utc, utc_to_jd, car2kep, ecef2eci_astropy
from utils.Propagators import sgp4_prop_TLE, numerical_prop, kepler_prop
from utils.SpaceObject import SpaceObject
from typing import List, Tuple

def SP3_file_validate():
    sp3.cddis.username = "charlesc"
    sp3.cddis.password = "0!EQVG8aDRWE"

    # Generate dates from 2022-01-01 to 2023-01-01 at 24-hour intervals
    dates = pd.date_range(start="2022-01-01", end="2022-01-02", freq='3H')
    #convert the start and end dates to JD
    jd_start = astropy.time.Time(dates[0]).jd
    jd_end = astropy.time.Time(dates[-1]).jd

    obstimes = astropy.time.Time(dates)

    # Get position in ITRS coordinate system
    sp3_itrs = sp3.itrs(
        id=sp3.NoradId("22824"), # Stella NORAD ID
        obstime=obstimes,
        download_directory="src/tests/SP3",
    )

    # Extract position and velocity data
    positions = sp3_itrs.cartesian.xyz.value.T  # Transpose to get shape (N, 3) #these are in meters
    #convert to km
    positions = positions/1000
    velocities = sp3_itrs.velocity.d_xyz.value.T  # Transpose to get shape (N, 3) #these are in kilometers per second already

    # Convert ECEF to ECI for each timestamp using ecef2eci_astropy
    mjds = obstimes.jd - 2400000.5  # Convert from JD to MJD

    eci_coords = [ecef2eci_astropy(pos, vel, mjd) for pos, vel, mjd in zip(positions, velocities, mjds)]

    print("eci_coords: ", eci_coords)

    print("first set of eci_coords: ", eci_coords[0])
    print("first position: ", eci_coords[0])
    print("first velocity: ", eci_coords[1])

    # Extract positions and velocities in the ECI frame
    positions_eci = np.squeeze(np.array([coord[0] for coord in eci_coords]))  # in km
    velocities_eci = np.squeeze(np.array([coord[1] for coord in eci_coords]))  # in km/s
    # Compute velocity magnitudes
    velocity_magnitudes = np.linalg.norm(velocities_eci, axis=1)
    #compute altitude
    altitudes_km = np.linalg.norm(positions_eci, axis=1) - 6378.137
    print("altitudes_km: ", altitudes_km)
    print("velocity_magnitudes: ", velocity_magnitudes)

    # Take the first position and convert it to Keplerian elements
    print("first position: ", positions_eci[0])
    print("first velocity: ", velocities_eci[0])

    posx = positions_eci[0][0]+3.6
    posy = positions_eci[0][1]+0.87
    posz = positions_eci[0][2]+1.49

    velx = velocities_eci[0][0]+0.00066
    vely = velocities_eci[0][1]+0.0005
    velz = velocities_eci[0][2]+0.0018

    print("posx: ", posx)
    print("posy: ", posy)
    print("posz: ", posz)
    print("velx: ", velx)
    print("vely: ", vely)
    print("velz: ", velz)
    
    a, e, i, w, W, V = car2kep(posx, posy, posz, velx, vely, velz)
    print(f"Keplerian elements: a={a}, e={e}, i={i}, w={w}, W={W}, V={V}")

    stella_sim = SpaceObject(sma = a, perigee=804, apogee=812, eccentricity=e, inc = np.rad2deg(i), argp = np.rad2deg(w), raan=np.rad2deg(W), tran=np.rad2deg(V), characteristic_area=2.5, mass = 48, epoch = "2022-01-01 00:00:00", launch_date='2022-01-01')
    print("jd start: ", jd_start)
    print("jd end: ", jd_end)

    # stella_sim.prop_catobject(jd_start=jd_start, jd_stop=jd_end, step_size=1000, output_freq=10800, integrator_type="RK45", force_model = ["all"])
    # stella_ephem = stella_sim.ephemeris
    # np.save("src/tests/sim_ephemeris/stella_ephem_3bp.npy", stella_ephem)

    #load starlette_ephem from .npy file
    stella_ephem = np.load("src/tests/sim_ephemeris/stella_ephem_3bp.npy") #this one does not have lunisolar perturbations
    print("stella_ephem: ", stella_ephem)

    #now the position difference between the SP3 and the starlette_ephem
    sp3pos = positions_eci[:-1]
    sp3vel = velocities_eci[:-1]
    sp3pos = sp3pos.reshape(-1, 3)
    sp3vel = sp3vel.reshape(-1, 3)

    stellapos = stella_ephem[:,0:3]
    stellavel = stella_ephem[:,3:6]
    a_stel, e_stel, i_stel, w_stel, W_stel, V_stel = car2kep(stellapos[0][0], stellapos[0][1], stellapos[0][2], stella_ephem[0][3], stella_ephem[0][4], stella_ephem[0][5])
    print("keplerian elements of first pos of stella_sim: ", a_stel, e_stel, i_stel, w_stel, W_stel, V_stel)

    #position difference at the first time step
    print("sp3pos[0]: ", sp3pos[0])
    print("starlettepos[0]: ", stellapos[0])
    print("sp3vel[0]: ", sp3vel[0])
    print("starlettevel[0]: ", stellavel[0])
    sp3_r = np.linalg.norm(sp3pos[0])
    stella_r = np.linalg.norm(stellapos[0])
    cart_pos_diff_t0 = np.linalg.norm(sp3_r - stella_r)
    print("pos_diff at start: ", cart_pos_diff_t0)

    # # Plot it over time
    obstimes_short = obstimes[:-1]
    #3D scatter plot of both positions
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.scatter(sp3pos[:,0], sp3pos[:,1], sp3pos[:,2], c='r', marker='o', s=1, label='SP3')
    ax.scatter(stellapos[:,0], stellapos[:,1], stellapos[:,2], c='b', marker='o', s=1, label='Starlette')
    ax.set_xlabel('X [km]')
    ax.set_ylabel('Y [km]')
    ax.set_zlabel('Z [km]')
    ax.legend()
    plt.show()

    #plot the altitude of the SP3 and starlette_ephem over time
    sp3alt = np.linalg.norm(sp3pos, axis=1) - 6378.137
    stellaalt = np.linalg.norm(stellapos, axis=1) - 6378.137   
    plt.plot(obstimes_short.datetime, sp3alt, label='SP3')
    plt.plot(obstimes_short.datetime, stellaalt, label='Starlette')
    plt.legend()
    plt.xlabel('Time')
    plt.ylabel('Altitude [km]')
    plt.show()

    #plot the velocity magnitude of the SP3 and starlette_ephem over time
    sp3vel_mag = np.linalg.norm(sp3vel, axis=1)
    stellavel_mag = np.linalg.norm(stellavel, axis=1)
    plt.plot(obstimes_short.datetime, sp3vel_mag, label='SP3')
    plt.plot(obstimes_short.datetime, stellavel_mag, label='Starlette')
    plt.legend()
    plt.xlabel('Time')
    plt.ylabel('Velocity magnitude [km/s]')
    plt.show()

    #3D/cartesian difference between the two positions over time
    cart_pos_diffs_all = []
    for i in range(len(sp3pos)):
        sp3_r = np.linalg.norm(sp3pos[i])
        stella_r = np.linalg.norm(stellapos[i])
        cart_pos_diff = np.linalg.norm(sp3_r - stella_r)
        print("cart_pos_diff: ", cart_pos_diff)
        print("sp3pos[i]: ", sp3pos[i])
        print("stellapos[i]: ", stellapos[i])
        cart_pos_diffs_all.append(cart_pos_diff)

    plt.plot(obstimes_short.datetime, cart_pos_diffs_all)
    plt.xlabel('Time')
    plt.ylabel('Cartesian position difference [km]')
    plt.show()


if __name__ == '__main__':
    SP3_file_validate()