import numpy as np
import pandas as pd
import astropy.time
import orekit
from astropy.time import Time
from org.orekit.orbits import CartesianOrbit, PositionAngle
from org.orekit.propagation import SpacecraftState
from org.orekit.time import AbsoluteDate
from org.orekit.propagation.analytical.tle import TLE
from org.orekit.utils import Constants, PVCoordinates
from org.hipparchus.geometry.euclidean.threed import Vector3D
from org.orekit.frames import FramesFactory
from org.orekit.propagation.conversion import TLEPropagatorBuilder, FiniteDifferencePropagatorConverter
from org.orekit.propagation.analytical.tle import TLEPropagator
from orekit.pyhelpers import setup_orekit_curdir, download_orekit_data_curdir
from java.util import ArrayList
from typing import List
from utils.Conversions import car2kep

def initialize_orekit():
    """Initializes Orekit."""
    download_orekit_data_curdir()
    orekit.initVM()
    setup_orekit_curdir()

def create_spacecraft_states(positions: List[List[float]], velocities: List[List[float]], dates_mjd: List[float]) -> ArrayList:
    """Creates a list of SpacecraftState objects based on positions, velocities, and dates.

    Args:
        positions (List[List[float]]): Positions.
        velocities (List[List[float]]): Velocities.
        dates_mjd (List[float]): Dates in Modified Julian Date format.

    Returns:
        ArrayList: List of SpacecraftState objects.
    """
    frame = FramesFactory.getICRF()
    spacecraft_states = ArrayList()
    for position, velocity, date_mjd in zip(positions, velocities, dates_mjd):
        date_orekit = AbsoluteDate(AbsoluteDate.MODIFIED_JULIAN_EPOCH, date_mjd * 86400.0)
        pv_coordinates = PVCoordinates(Vector3D(*position), Vector3D(*velocity))
        orbit = CartesianOrbit(pv_coordinates, frame, date_orekit, Constants.EIGEN5C_EARTH_MU)
        state = SpacecraftState(orbit)
        spacecraft_states.add(state)
    return spacecraft_states

def fit_tle_to_spacecraft_states(spacecraft_states: ArrayList, satellite_number: int, classification: str,
                                 launch_year: int, launch_number: int, launch_piece: str, ephemeris_type: int,
                                 element_number: int, date_start_orekit: AbsoluteDate, mean_motion: float,
                                 mean_motion_first_derivative: float, mean_motion_second_derivative: float, e: float,
                                 i: float, pa: float, raan: float, ma: float, revolution_number: int,
                                 b_star_first_guess: float) -> TLE:
    """Fits a TLE to a given list of SpacecraftState objects.

    Args:
        spacecraft_states (ArrayList): Spacecraft states.
        satellite_number (int): Satellite number.
        classification (str): Classification.
        launch_year (int): Launch year.
        launch_number (int): Launch number.
        launch_piece (str): Launch piece.
        ephemeris_type (int): Ephemeris type.
        element_number (int): Element number.
        date_start_orekit (AbsoluteDate): Start date in Orekit format.
        mean_motion (float): Mean motion.
        mean_motion_first_derivative (float): Mean motion first derivative.
        mean_motion_second_derivative (float): Mean motion second derivative.
        e (float): Eccentricity.
        i (float): Inclination.
        pa (float): Argument of perigee.
        raan (float): Right ascension of the ascending node.
        ma (float): Mean anomaly.
        revolution_number (int): Revolution number.
        b_star_first_guess (float): B* drag term first guess.

    Returns:
        TLE: Fitted TLE.
    """
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
    threshold = 1.0 #TODO: what is the physical meaning of this threshold?
    tle_builder = TLEPropagatorBuilder(tle_first_guess, PositionAngle.MEAN, 1.0)
    fitter = FiniteDifferencePropagatorConverter(tle_builder, threshold, 1000)
    fitter.convert(spacecraft_states, False, 'BSTAR')
    tle_propagator = TLEPropagator.cast_(fitter.getAdaptedPropagator())
    return tle_propagator.getTLE()

def generate_dates(mjds: list) -> pd.DatetimeIndex:
    """Generates dates from a given list of Modified Julian Dates (MJDs).

    Args:
        mjds (list): List of Modified Julian Dates.

    Returns:
        pd.DatetimeIndex: Generated dates.
    """
    utc_times = [Time(mjd, format="mjd", scale="utc").datetime for mjd in mjds]
    return pd.DatetimeIndex(utc_times)

def fit_TLE_to_ephemeris(positions_eci, velocities_eci, jds):

    a, e, i, pa, raan, ma = car2kep(*positions_eci[0], *velocities_eci[0])
    e = float(np.deg2rad(e))
    i = float(np.deg2rad(i))
    pa = float(np.deg2rad(pa))
    raan = float(np.deg2rad(raan))
    ma = float(np.deg2rad(ma))

    mjds = jds - 2400000.5
    dates = generate_dates(mjds)
    obstimes = astropy.time.Time(dates)
    mjds = obstimes.jd - 2400000.5

    # Create spacecraft states
    spacecraft_states = create_spacecraft_states(positions_eci, velocities_eci, mjds)
    ## Placeholder parameters.
    ## TODO: I believe none of these are actually used in the propagtion itself but they are required to make a TLE
    ## TODO: we must double check that none of these are used in the propagation itself.
    satellite_number = 99999 #TODO: this is going to have to be changed to self.norad id
    classification = 'U' #
    launch_year = 2021
    launch_number = 42
    launch_piece = 'A'
    ephemeris_type = 0
    element_number = 999
    mean_motion = 14.322151  # Example mean motion in revs per day
    mean_motion_first_derivative = 0.0
    mean_motion_second_derivative = 0.0
    revolution_number = 12345

    date_start_orekit = AbsoluteDate(AbsoluteDate.MODIFIED_JULIAN_EPOCH, mjds[0]*86400.0)  # set the start date of the TLE
    b_star_first_guess = 1e-5 # doesn't matter what this is set to, it will be fit to the spacecraft states

    # Call the function to fit TLE
    fitted_tle = fit_tle_to_spacecraft_states(spacecraft_states, satellite_number, classification,
                                            launch_year, launch_number, launch_piece, ephemeris_type,
                                            element_number, date_start_orekit, mean_motion,
                                            mean_motion_first_derivative, mean_motion_second_derivative,
                                            e, i, pa, raan, ma, revolution_number, b_star_first_guess)

    print("Fitted TLE:", fitted_tle)

if __name__ == "__main__":
    initialize_orekit() #TODO: this is going to have to be activated in the SpaceCatalogue class if sgp4 is used so that it is not reinitialized every time
    #test example
    # Generate dates
    mjds = [60066.6128, 60067.3451]
    positions_eci = [[-2537.205, 6342.796, 0.0], [-2437.205, 6542.796, 123.0]]
    velocities_eci = [[-1.937, -0.72, 7.361], [-2.937, -5.72, 0.361]]
    fit_TLE_to_ephemeris(positions_eci, velocities_eci, mjds)
    