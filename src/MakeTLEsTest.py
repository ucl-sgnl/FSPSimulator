import numpy as np
import pandas as pd
import sp3
import astropy.time
import astropy.units as u
import orekit
from astropy.coordinates import CartesianDifferential, CartesianRepresentation, GCRS, ITRS
from astropy.time import Time
from org.orekit.bodies import CelestialBodyFactory
from org.orekit.data import DataProvidersManager, ZipJarCrawler
from org.orekit.estimation.leastsquares import BatchLSEstimator
from org.orekit.estimation.measurements import ObservableSatellite, Position
from org.orekit.forces.gravity.potential import GravityFieldFactory
from org.orekit.orbits import CartesianOrbit, KeplerianOrbit, PositionAngle
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
from typing import List, Tuple
from utils.Conversions import car2kep, ecef2eci_astropy

sp3.cddis.username = "charlesc"
sp3.cddis.password = "0!EQVG8aDRWE"


def generate_dates(start: str, end: str, freq: str = '3H') -> pd.DatetimeIndex:
    """Generates dates within a given range at specified intervals.

    Args:
        start (str): Start date.
        end (str): End date.
        freq (str, optional): Frequency of intervals. Defaults to '3H'.

    Returns:
        pd.DatetimeIndex: Generated dates.
    """
    return pd.date_range(start=start, end=end, freq=freq)


def get_sp3_itrs_data(obstimes: Time, norad_id: str, download_directory: str) -> Tuple[np.array, np.array]:
    """Gets the SP3 ITRS data (positions and velocities).

    Args:
        obstimes (Time): Observation times.
        norad_id (str): NORAD ID.
        download_directory (str): Directory to download data.

    Returns:
        Tuple[np.array, np.array]: Positions and velocities.
    """
    sp3_itrs = sp3.itrs(
        id=sp3.NoradId(norad_id),
        obstime=obstimes,
        download_directory=download_directory,
    )
    positions = sp3_itrs.cartesian.xyz.value.T / 1000
    velocities = sp3_itrs.velocity.d_xyz.value.T
    return positions, velocities


def convert_ecef_to_eci(positions: np.array, velocities: np.array, mjds: np.array) -> Tuple[np.array, np.array]:
    """Converts ECEF coordinates to ECI.

    Args:
        positions (np.array): Positions in ECEF.
        velocities (np.array): Velocities in ECEF.
        mjds (np.array): Modified Julian Dates.

    Returns:
        Tuple[np.array, np.array]: Positions and velocities in ECI.
    """
    eci_coords = [ecef2eci_astropy(pos, vel, mjd) for pos, vel, mjd in zip(positions, velocities, mjds)]
    positions_eci = np.squeeze(np.array([coord[0] for coord in eci_coords]))
    velocities_eci = np.squeeze(np.array([coord[1] for coord in eci_coords]))
    return positions_eci, velocities_eci


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
    gcrf = FramesFactory.getGCRF()
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
    threshold = 1.0
    tle_builder = TLEPropagatorBuilder(tle_first_guess, PositionAngle.MEAN, 1.0)
    fitter = FiniteDifferencePropagatorConverter(tle_builder, threshold, 1000)
    fitter.convert(spacecraft_states, False, 'BSTAR')
    tle_propagator = TLEPropagator.cast_(fitter.getAdaptedPropagator())
    return tle_propagator.getTLE()


if __name__ == "__main__":
    #test example
    # Generate dates
    dates = generate_dates(start="2022-01-01", end="2022-01-02")
    obstimes = astropy.time.Time(dates)
    # Get SP3 data
    positions, velocities = get_sp3_itrs_data(obstimes, norad_id="22824", download_directory="src/tests/SP3")
    # Convert ECEF to ECI
    mjds = obstimes.jd - 2400000.5
    positions_eci, velocities_eci = convert_ecef_to_eci(positions, velocities, mjds)
    # Initialize Orekit
    initialize_orekit()
    # Create spacecraft states
    spacecraft_states = create_spacecraft_states(positions_eci, velocities_eci, mjds)
    ## Example parameters
    spacecraft_states = spacecraft_states  # Assuming this is already defined in your code
    satellite_number = 99999
    classification = 'U'
    launch_year = 2021
    launch_number = 42
    launch_piece = 'A'
    ephemeris_type = 0
    element_number = 999
    date_start_orekit = AbsoluteDate(AbsoluteDate.MODIFIED_JULIAN_EPOCH, 58923*86400.0)  # Example date
    mean_motion = 14.322151  # Example mean motion in revs per day
    mean_motion_first_derivative = 0.0
    mean_motion_second_derivative = 0.0
    e = 0.001
    i = float(np.deg2rad(98.0))
    pa = float(np.deg2rad(42.0))
    raan = float(np.deg2rad(42.0))
    ma = float(np.deg2rad(42.0))
    revolution_number = 12345
    b_star_first_guess = 1e-5

    # Call the function to fit TLE
    fitted_tle = fit_tle_to_spacecraft_states(spacecraft_states, satellite_number, classification,
                                            launch_year, launch_number, launch_piece, ephemeris_type,
                                            element_number, date_start_orekit, mean_motion,
                                            mean_motion_first_derivative, mean_motion_second_derivative,
                                            e, i, pa, raan, ma, revolution_number, b_star_first_guess)

    print("Fitted TLE:", fitted_tle)
