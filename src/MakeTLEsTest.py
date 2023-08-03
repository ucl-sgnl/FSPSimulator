#import SP3 orbit
#fit a TLE to it
#compare the TLE orbit positions to the SP3 orbit

import numpy as np
import datetime
import matplotlib.pyplot as plt
import sp3
import astropy.time
import astropy.coordinates
import pandas as pd
import astropy.units as u
from utils.Conversions import car2kep
from astropy import units as u
from astropy.coordinates import CartesianDifferential, CartesianRepresentation, GCRS, ITRS
from astropy.time import Time
import unittest
from utils.Conversions import tle_parse, tle_convert, TLE_time, jd_to_utc, utc_to_jd, car2kep, ecef2eci_astropy
from utils.Propagators import sgp4_prop_TLE, numerical_prop, kepler_prop
from utils.SpaceObject import SpaceObject
from typing import List, Tuple

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

print("first sp3 position eci: ", eci_coords[0])
print("first sp3 velocity eci: ", eci_coords[1])

# Extract positions and velocities in the ECI frame
positions_eci = np.squeeze(np.array([coord[0] for coord in eci_coords]))  # in km
velocities_eci = np.squeeze(np.array([coord[1] for coord in eci_coords]))  # in km/s

####### Now we have times, positions and velocities for SP3 orbit, fit a TLE to them using Orekit ######

import orekit
from org.orekit.bodies import CelestialBodyFactory
from org.orekit.data import DataProvidersManager, ZipJarCrawler
from org.orekit.estimation.leastsquares import BatchLSEstimator
from org.orekit.estimation.measurements import ObservableSatellite, Position
from org.orekit.forces.gravity.potential import GravityFieldFactory
from org.orekit.orbits import KeplerianOrbit, PositionAngle
from orekit.pyhelpers import setup_orekit_curdir, download_orekit_data_curdir
from org.orekit.frames import FramesFactory
from org.orekit.utils import Constants
from org.hipparchus.geometry.euclidean.threed import Vector3D
from org.orekit.time import AbsoluteDate
from org.orekit.utils import PVCoordinates
from org.orekit.propagation import SpacecraftState
from org.orekit.frames import FramesFactory
from org.orekit.bodies import CelestialBodyFactory
from org.orekit.orbits import CartesianOrbit
# Initialize Orekit
download_orekit_data_curdir()
orekit.initVM()
setup_orekit_curdir()

# Your positions, velocities, and dates in MJD format
positions = [[7000000.0, 1000000.0, 2000000.0], [7100000.0, 900000.0, 1900000.0]]
velocities = [[2000.0, 3000.0, 1000.0], [2100.0, 2900.0, 900.0]]
dates_mjd = [58923, 58924]  # Modify with your actual MJD dates

# GCRF frame (geocentric celestial reference frame)
gcrf = FramesFactory.getGCRF()
frame = FramesFactory.getICRF()

# Create a list to store the SpacecraftState objects
from java.util import ArrayList
spacecraft_states = ArrayList()

# Iterate through the provided data to create the SpacecraftState objects
for position, velocity, date_mjd in zip(positions, velocities, dates_mjd):
    # Convert MJD to Orekit AbsoluteDate
    date_orekit = AbsoluteDate(AbsoluteDate.MODIFIED_JULIAN_EPOCH, date_mjd*86400.0)

    # Create PVCoordinates from position and velocity
    pv_coordinates = PVCoordinates(Vector3D(*position), Vector3D(*velocity))

    # Create a Cartesian Orbit object
    orbit = CartesianOrbit(pv_coordinates, frame, date_orekit, Constants.EIGEN5C_EARTH_MU)

    # Create a SpacecraftState object and add to the list
    state = SpacecraftState(orbit)
    spacecraft_states.add(state)

#Keplerian Elements
a = 7000e3
e = 0.001
i = float(np.deg2rad(98.0)) #conversion to float required for orekit
pa = float(np.deg2rad(42.0))
raan = float(np.deg2rad(42.0))
ma = float(np.deg2rad(42.0))

#Satellite information
satellite_number = 99999
classification = 'X'
launch_year = 2018
launch_number = 42
launch_piece = 'F'
ephemeris_type = 0
element_number = 999
revolution_number = 100

#Numerical Propagator parameters
dt = 60.0  # s, period at which the spacecraft states are saved to fit the TLE

prop_min_step = 0.001 # s
prop_max_step = 300.0 # s
prop_position_error = 10.0 # m

# Estimator parameters
estimator_position_scale = 1.0 # m
estimator_convergence_thres = 1e-3
estimator_max_iterations = 25
estimator_max_evaluations = 35

#spacecraft properties
sc_mass = 400.0 # kg
sc_cross_section = 0.3 # m2
cd_drag_coeff = 2.0 
cr_radiation_pressure = 1.0

#date determines solar activity and drag and fitting duration is mean elements being fitted to keplerian elements
from datetime import datetime
date_start = datetime(2019, 1, 1)
fitting_duration_d = 1 # days

from org.orekit.orbits import KeplerianOrbit, PositionAngle
from org.orekit.utils import Constants as orekit_constants
from orekit.pyhelpers import datetime_to_absolutedate

from org.orekit.frames import FramesFactory, ITRFVersion
from org.orekit.utils import IERSConventions
gcrf = FramesFactory.getGCRF()
teme = FramesFactory.getTEME()
itrf = FramesFactory.getITRF(IERSConventions.IERS_2010, False)

from org.orekit.models.earth import ReferenceEllipsoid
wgs84_ellipsoid = ReferenceEllipsoid.getWgs84(itrf)

from org.orekit.bodies import CelestialBodyFactory
moon = CelestialBodyFactory.getMoon()
sun = CelestialBodyFactory.getSun()

from org.orekit.orbits import KeplerianOrbit, PositionAngle
from org.orekit.utils import Constants as orekit_constants
from orekit.pyhelpers import datetime_to_absolutedate
date_start_orekit = datetime_to_absolutedate(date_start)
keplerian_orbit = KeplerianOrbit(a, e, i, pa, raan, ma, PositionAngle.MEAN, 
                                 gcrf, date_start_orekit, orekit_constants.EIGEN5C_EARTH_MU)

print("keplerian orbit:", keplerian_orbit)

from org.orekit.propagation.analytical.tle import TLE
mean_motion = float(np.sqrt(orekit_constants.EIGEN5C_EARTH_MU / np.power(a, 3)))
mean_motion_first_derivative = 0.0
mean_motion_second_derivative = 0.0
b_star_first_guess = 1e-5  # Does not play any role, because it is a free parameter when fitting the TLE

tle_first_guess = TLE(satellite_number, 
                        classification,
                        launch_year,
                        launch_number,
                        launch_piece,
                        ephemeris_type,
                        element_number,
                        date_start_orekit,
                        mean_motion,
                        mean_motion_first_derivative, 
                        mean_motion_second_derivative,
                        e,
                        i,
                        pa,
                        raan,
                        ma,
                        revolution_number,
                        b_star_first_guess)

print("tle first guess:", tle_first_guess)
print("spacecraft states to fit TLE to:", spacecraft_states)


from org.orekit.propagation.conversion import TLEPropagatorBuilder, FiniteDifferencePropagatorConverter
from org.orekit.propagation.analytical.tle import TLEPropagator
threshold = 1.0  # "absolute threshold for optimization algorithm", but no idea about its impact
tle_builder = TLEPropagatorBuilder(tle_first_guess, PositionAngle.MEAN, 1.0)
fitter = FiniteDifferencePropagatorConverter(tle_builder, threshold, 1000)
fitter.convert(spacecraft_states, False, 'BSTAR')  # Setting BSTAR as free parameter
tle_propagator = TLEPropagator.cast_(fitter.getAdaptedPropagator())
tle_fitted = tle_propagator.getTLE()

print("fitted TLE:", tle_fitted)