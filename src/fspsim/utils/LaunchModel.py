import math
from datetime import datetime
from datetime import timedelta
import os
import json 
import numpy as np
from tqdm import tqdm
from datetime import timedelta
from ..utils.Conversions import orbital_period
from ..utils.SpaceObject import SpaceObject
from ..utils.Formatting import future_constellations_csv_handler

def import_configuration_json(filename):
    with open(filename) as f:
        data = json.load(f)
    return data

LEO_launchers = {
    'Falcon 9': {
        'capacity': 22.8,
        'cost_per_kg': 2700
    },
    'Falcon Heavy': {
        'capacity': 63.8,
        'cost_per_kg': 1400
    },
    'Electron': {
        'capacity': 0.3,
        'cost_per_kg': 25000
    },
    'Terran 1': {
        'capacity': 1.2,
        'cost_per_kg': 6000
    },
    'Vega': {
        'capacity': 1.5,
        'cost_per_kg': 20000
    },
    'Ariane 6': {
        'capacity': 21.6,
        'cost_per_kg': 8000
    },
    'Soyuz-2': {
        'capacity': 7.8,
        'cost_per_kg': 4000
    },
    'Delta IV Heavy': {
        'capacity': 28.4,
        'cost_per_kg': 15000
    },
    'Atlas V': {
        'capacity': 18.8,
        'cost_per_kg': 9000
    },
    'New Glenn': {
        'capacity': 45,
        'cost_per_kg': 2200
    }
}

def global_launch_schedule(sub_constellation_metadata_dicts, monthly_ton_capacity, launches_start_date, rocket = "Falcon 9"):
    """
    Determines the launch dates for various sub-constellations based on a greedy algorithm, limited by
    global montly launch capacity (tons/month), the date to start scheduling launches, and the type of rocket used.

    :param sub_constellation_metadata_dicts: List of metadata for each sub-constellation.
    :type sub_constellation_metadata_dicts: list[dict]
    :param monthly_ton_capacity: The maximum weight capacity available for launching in tons per month.
    :type monthly_ton_capacity: float
    :param launches_start_date: The initial date to begin scheduling launches.
    :type launches_start_date: str
    :param rocket: The type of rocket used for launches, defaults to "Falcon 9".
    :type rocket: str, optional
    :raises ValueError: If the provided rocket is not in the list of LEO launchers.
    :return: A dictionary of launch dates for each sub-constellation.
    :rtype: dict[str, list[str]]
    """
    if rocket not in LEO_launchers:
        raise ValueError(f'rocket must be one of the following: {list(LEO_launchers.keys())}')

    total_mass = sum(int(sub_constellation['_mass']) * int(sub_constellation['N']) for sub_constellation in sub_constellation_metadata_dicts)
    max_ton_per_launch = LEO_launchers[rocket]['capacity']
    monthly_launches_frequency = math.ceil(monthly_ton_capacity / max_ton_per_launch)
    months_required = math.ceil(total_mass / 1000 / monthly_ton_capacity)

    sub_constellation_launches_required = {
        sub_constellation['_soname']: math.ceil((sub_constellation['_mass'] * sub_constellation['N']) / 1000 / monthly_ton_capacity) 
        for sub_constellation in sub_constellation_metadata_dicts
    }

    # Precompute all possible launch dates
    launch_bank = [
        (datetime.strptime(launches_start_date, "%Y-%m-%d") + timedelta(days=i*30) + timedelta(days=j*30/monthly_launches_frequency)).strftime("%Y-%m-%d") 
        for i in range(months_required) 
        for j in range(monthly_launches_frequency)
    ]

    # Sort the constellations by the number of launches required
    sorted_constellations = sorted(sub_constellation_launches_required.items(), key=lambda item: item[1], reverse=True)

    # Allocate launch dates to constellations
    sub_constellation_launch_dates = {}
    for constellation, num_launches in sorted_constellations:
        sub_constellation_launch_dates[constellation] = launch_bank[:num_launches]
        launch_bank = launch_bank[num_launches:]

    return sub_constellation_launch_dates


def create_subconstellation_Space_Objects(N, i, h, _soname, _application, _owner, launch_schedule, _mass, _area, _length, _maneuverable, _propulsion):
    """
    Generate a list of SpaceObjects that are in a Walker-Delta constellation.
    This maximizes geometric coverage of the Earth given the number of satellites and the inclination and altitude of the orbit.
    Based on the parameters passed in the launch file.

    :param N: Number of satellites in the constellation.
    :type N: int
    :param i: Inclination of the satellite orbit.
    :type i: float
    :param h: Altitude of the satellite orbit.
    :type h: float
    :param _soname: Name of the sub-constellation.
    :type _soname: str
    :param _application: Application for which the satellite is used.
    :type _application: str
    :param _owner: Owner or operator of the satellite.
    :type _owner: str
    :param launch_schedule: Schedule of satellite launches.
    :type launch_schedule: list[str]
    :param _mass: Mass of the satellite.
    :type _mass: float
    :param _area: Characteristic area of the satellite.
    :type _area: float
    :param _length: Characteristic length of the satellite.
    :type _length: float
    :param _maneuverable: Indicates if the satellite is maneuverable.
    :type _maneuverable: str
    :param _propulsion: Type of propulsion used in the satellite.
    :type _propulsion: str
    :raises ValueError: If N is less than 1.
    :return: A list of space objects for the sub-constellation.
    :rtype: list[SpaceObject]
    """
    # Check that N is always >= 1 
    if N < 1:
        raise ValueError('N (number of satellites in constellation) must be >= 1')

    # Constants and initializations
    apogee_alt = perigee_alt = h  # assuming circular orbit
    R_earth = 6371  # Earth's mean radius in km
    a = R_earth + h  # semi-major axis in km

    P = math.ceil(math.sqrt(N))  # number of planes
    S = math.ceil(N / P)  # number of satellites per plane
    delta_Omega = 360 / P
    delta_omega = 360 / S
    delta_TRAN = 360 / N
    period = orbital_period(a)

    if any(math.isnan(val) for val in [a, P, S, delta_Omega, delta_omega, delta_TRAN, period]):
        print("ERROR: NaN value in FSP launch entry")
        return

    # Repeat and truncate launch_schedule 
    launch_schedule = launch_schedule * (N // len(launch_schedule)) + launch_schedule[:N % len(launch_schedule)]

    # Pre-compute decay_schedule
    dates_np = np.array([datetime.strptime(date, "%Y-%m-%d") for date in launch_schedule])
    decay_schedule = [(date + timedelta(days=365*5)).strftime("%Y-%m-01") for date in dates_np]

    if _maneuverable == 'y':
        station_keeping = True
    else:
        station_keeping = False
    subconstellation_Space_Objects = [
        SpaceObject(
            object_type="PAY",
            payload_operational_status="+",
            application=_application,
            operator=_owner,
            apogee=apogee_alt,
            perigee=perigee_alt,
            mass=_mass,
            maneuverable='y',
            spin_stabilized='n',
            characteristic_area=_area,
            characteristic_length=_length,
            propulsion_type=_propulsion,
            sma=a,
            eccentricity=0,
            inc=np.deg2rad(i),
            argp=np.deg2rad((n % S) * delta_omega),
            raan=np.deg2rad(math.floor((n-1) / S) * delta_Omega),
            tran=np.deg2rad((n-1) * delta_TRAN),
            launch_site="JSC",
            launch_date=launch_schedule[n-1],
            decay_date=decay_schedule[n-1],
            rso_name=f"{_soname}_{n}",
            station_keeping=station_keeping,
        ) 
        for n in range(1, N+1)
    ]

    return subconstellation_Space_Objects

def Prediction2SpaceObjects(satellite_predictions_csv, simsettings):
    """
    Generates a list SpaceObjects for each of the sub-constellation in the prediction data CSV file.
    Each object gets assigned a launch and decay date based on the launch model parameters and the parameters 
    specified in the prediction data CSV file.

    :param satellite_predictions_csv: CSV file with the satellite predictions data.
    :type satellite_predictions_csv: str
    :param simsettings: Loaded JSON file with the simulation settings (contains launch model parameters).
    :type simsettings: dict
    :return: A list of space objects generated from the prediction data.
    :rtype: list[SpaceObject]
    """
    all_space_objects = []
    # create a list of dictionaries containing the metadata for each sub-constellation
    # if the user has selected a file use this
    if satellite_predictions_csv:
        metadata_dicts = future_constellations_csv_handler(file_path=satellite_predictions_csv)
    else:
        lauch_file_direc = 'src/fspsim/data/prediction_csv'
        metadata_dicts = future_constellations_csv_handler(file_path=os.path.join(lauch_file_direc, satellite_predictions_csv))

    #global_launch_schedule(sub_constellation_metadata_dicts=metadata_dicts)
    sub_constellation_launch_dates = global_launch_schedule(
                                                            sub_constellation_metadata_dicts = metadata_dicts, 
                                                            monthly_ton_capacity=float(simsettings['monthly_ton_capacity']), 
                                                            launches_start_date = simsettings['launch_start_date']   
                                                         )                                                     
    
    for dict in tqdm(metadata_dicts): #TODO: use the JSON file to set the agressivity and max launch rate etc.
        subconstellation_Space_Objects = create_subconstellation_Space_Objects(N=int(dict['N']), i = float(dict['i']), h=float(dict['h']), _soname=dict['_soname'], _application = dict['_application'], _owner= dict['_owner'], launch_schedule = sub_constellation_launch_dates[dict['_soname']], _mass=dict['_mass'], _area=dict['_area'], _length=dict['_length'], _maneuverable= dict['_maneuverable'], _propulsion=dict['_propulsion'])
        all_space_objects.extend(subconstellation_Space_Objects)
    return all_space_objects

if __name__ == '__main__':
    def get_path(*args):
        return os.path.join(os.getcwd(), *args)

    in_file = 'src/fspsim/data/prediction_csv/04_04_23_fsp.csv'
    settings = json.load(open(get_path('src/fspsim/data/specify_simulation/testsim.json'), 'r')) 
    Prediction2SpaceObjects(in_file, settings)