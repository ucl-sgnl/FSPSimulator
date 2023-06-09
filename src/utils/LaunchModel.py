import os
import math
import pandas as pd
from datetime import datetime
from datetime import timedelta
import os
import json 
import csv
import random
import numpy as np
from datetime import timedelta
from utils.Conversions import orbit_classify, orbital_period, generate_cospar_id
from utils.SpaceObject import SpaceObject

def import_configuration_json(filename):
    with open(filename) as f:
        data = json.load(f)
    return data

def calculate_form_factor(form_factor_str):
    """reads a string describing the form factor of satellties in a sub constellation and 
    returns characteristic length and area of to populate the SpaceObject class metadata
    return error if form factor is not a string

    Args:
        form_factor_str (string): string describing the form factor of satellties

    Raises:
        ValueError: Form factor must be a string

    Returns:
        tuple: characteristic length, characteristic area
    """
    if not isinstance(form_factor_str, str):
        raise ValueError('Form factor must be a string')
    
    if 'u' in form_factor_str: #DIVIDE by 10 since 1U = 10cm^3
        _length = float(form_factor_str.split('u')[0]) /10 
        _area = float(form_factor_str.split('u')[0]) /10
    elif 'U' in form_factor_str:
        _length = float(form_factor_str.split('U')[0]) /10
        _area = float(form_factor_str.split('U')[0]) /10
    # if there is the character '*', extract all the numbers. The largest value is the characteristic length and the characteristic area is all the two largest values multiplied together
    elif '*' in form_factor_str:
        #split the values by the character '*' and make them in ascending order
        form_factor_str = sorted([float(i) for i in form_factor_str.split('*')])
        _length = form_factor_str[-1]
        _area = form_factor_str[-1] * form_factor_str[-2]

    return _length, _area
    
    #####======= THIS NEEDS TO BE ITS OWN FUNCTION =========#####
    #####======= DECIMATE() =========#####
    # # create a list of the indices of the rows in the dataframe
    # indices = list(sat_df.index)
    # # randomly shuffle the indices
    # np.random.shuffle(indices)
    # # calculate the number of rows to drop
    # num_rows_to_drop = math.ceil(len(indices) * startup_failure_rate)
    # # drop the rows
    # sat_df = sat_df.drop(indices[:num_rows_to_drop])
    # #print the number of rows remaining
    # print('number of rows remaining after startup failure rate applied:',sat_df.shape[0])
    #####======= THIS NEEDS TO BE ITS OWN FUNCTION =========#####

def satellite_metadata(file_path):
    """Reads in the excel file predictions and generates a dicitonary of metadate for each subconsetllation

    Args:
        file_path (str): path to the csv file with all the satellite data

    Returns:
        _type_: _description_
    """
    #TODO: Outline exact units for each of the parameters 
    #function that reads in the excel file predictions and generates a dicitonary of metadate for each subconsetllation
    # this subconstellation metadata will then be used to generate the correct SpaceObject instances later on

    sat_df = pd.read_csv(file_path, sep=',') # this is the csv file with all the satellite data

    # slice the dataframe to remove all the rows that have missing values or nan values in the 'Number of sats', 'Inclination', 'Altitude', 'Sub-Constellation', 'Mission type/application', 'Mass(kg)', 'Form Factor' columns
    sat_df = sat_df.dropna(subset=['Number of sats', 'Inclination', 'Altitude', 'Sub-Constellation', 'Mission type/application', 'Mass(kg)', 'Form Factor', 'Maneuverable','Propulsion'])
    #print the number of rows remaining
    print('number of rows with all data required in prediction CSV:',sat_df.shape[0])

    metadata_dicts = [] # this is a list of dictionaries that will contain the metadata for each row (sub-constellation)
    for index, row in sat_df.iterrows():
        subconstellation_dict = {} # this is a dictionary that will contain the metadata for each sub-constellation
        N = row['Number of sats']
        i = row['Inclination']
        h = row['Altitude']
        _maneuverable = row['Maneuverable']
        _propulsion = row['Propulsion']
        operator_sc_name = row['Sub-Constellation']
        #strip the white spaces from the operator_sc_name and replace the spaces with underscores
        _soname = operator_sc_name.replace(' ','_')
        _application = row['Mission type/application']
        _mass = row['Mass(kg)']
        form_factor = row['Form Factor']
        _owner = row['Mega-Constellation']
        # in the form factor row, extract the number before the letter u and multiply it by 10 for the length and 100 for the area
        # the length and area are the number before the letter u times 10 and 100 respectively for a cubesat
        _length, _area = calculate_form_factor(form_factor)
        
        # populate the subconstellation_dict with the metadata
        subconstellation_dict['N'] = N
        subconstellation_dict['i'] = i
        subconstellation_dict['h'] = h
        subconstellation_dict['_soname'] = _soname
        subconstellation_dict['_owner'] = _owner
        subconstellation_dict['_application'] = _application
        subconstellation_dict['_mass'] = _mass
        subconstellation_dict['_maneuverable'] = _maneuverable
        subconstellation_dict['_propulsion'] = _propulsion
        # convert from cm to m
        subconstellation_dict['_length'] = _length
        # convert from cm^2 to m^2
        subconstellation_dict['_area'] = _area
        # append the subconstellation_dict to the metadata_dicts list
        metadata_dicts.append(subconstellation_dict)

    return metadata_dicts

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

def global_launch_schedule(sub_constellation_metadata_dicts, settings, monthly_ton_capacity, launches_start_date, rocket = "Falcon 9"):
    """Given a set of sub-constellation metadata dictionaries, 
    this function generates a global launch schedule for all the sub-constellations

    Args:
        sub_constellation_metadata_dicts (array-like): list of dictionaries containing the metadata for each sub-constellation as generated by the satellite_metadata() function
        monthly_ton_capacity (int, optional): monthly launch capacity in Tons. Defaults to 40.#TODO: add consideration for LEO, MEO, GEO launch capabilities
        launches_start_date (str, optional): Beggining of the launch schedule. Defaults to "2023-05-01".

    Returns:
        dict: launch_schedule_dict
    """

    # monthly ton capacity taken from: https://www.spacexstats.xyz/#payloads-upmass-per-year
    possible_launchers = list(LEO_launchers.keys())
    if rocket not in possible_launchers:
        raise ValueError(f'rocket must be one of the following: {possible_launchers}')
    #calculate total mass of all constellations to put in orbit
    total_mass = 0
    total_cost = 0
    for sub_constellation in sub_constellation_metadata_dicts:
        total_mass += int(sub_constellation['_mass']) * int(sub_constellation['N']) #need to make sure the units here are always kg
        sub_constellation['total_mass'] = int(sub_constellation['_mass']) * int(sub_constellation['N'])
        # print(f'total mass of:',sub_constellation['_soname'],':',(sub_constellation['_mass'] * sub_constellation['N'])/1000,'(T)')

        total_cost += (sub_constellation['_mass'] * sub_constellation['N']) * LEO_launchers[rocket]['cost_per_kg']
        sub_constellation['total_cost'] = (sub_constellation['_mass'] * sub_constellation['N']) * LEO_launchers[rocket]['cost_per_kg']
        # print(f"cost in Millions of USD:",total_cost/1000000)

    # calculate the number of launches required to put all the satellites in orbit based on the monthly ton capacity
    total_mass_tons = total_mass / 1000 # convert from kg to tons
    max_ton_per_launch = LEO_launchers[rocket]['capacity']
    monthly_launches_frequency = math.ceil(monthly_ton_capacity / max_ton_per_launch) 
    print('max number of launches per month:',monthly_launches_frequency)
    
    # now how many months are required to launch the total mass of satellites
    months_required = math.ceil(total_mass_tons / monthly_ton_capacity)
    print('number of months required to put all satellites in orbit:',months_required)

    # dictionary of the number of launches required for each sub constellation
    # create a dictionary of the number of launches required for each sub constellation 
    sub_constellation_launches_required = {}
    for sub_constellation in sub_constellation_metadata_dicts:
        constellation_mass = sub_constellation['_mass'] * sub_constellation['N']
        constellation_launch_cost = constellation_mass * LEO_launchers[rocket]['cost_per_kg']
        constellation_mass_tons = constellation_mass / 1000
        no_launches_required = math.ceil(constellation_mass_tons / monthly_ton_capacity)
        #TODO: instead of printing here write this to a csv file with the sub constellation name, mass, number of launches required and cost
        # print(f"{sub_constellation['_soname']}, mass: {constellation_mass_tons} T, #Launches: {no_launches_required}, M$: {constellation_launch_cost/1000000} M")
        sub_constellation_launches_required[sub_constellation['_soname']] = no_launches_required

    # now we need to calculate the launch dates for each sub constellation
    # make a list of all launches and give them dates dates based on the months_required and monthly_launches_frequency and the launch_start_date. 
    # we will attribute the launches to the sub constellations in the next step
    launch_bank = [] #bank of available launches
    for i in range(months_required):
        for j in range(monthly_launches_frequency):
            # use datetime and timedelta to calculate the launch dates
            launch_date = datetime.strptime(launches_start_date, "%Y-%m-%d") + timedelta(days=i*30) + timedelta(days=j*30/monthly_launches_frequency) # 
            launch_bank.append(launch_date)
    
    # now we need to attribute the launches to the sub constellations
    # to attribute the launches to the sub constellations, we will use a kind of greedy algorithm (i think its not really a greedy algorithm)
    # we will start with the sub constellation that requires the most launches and assign it the first launches in the launch bank
    # go through sub_constellation_launches_required and find the sub_constellations with the most launches required
    # then assign the first launches in the launch bank to that sub constellation 
    # then remove the launches from the launch bank and the lanches required from the sub_constellation_launches_required
    # then repeat until all launches have been assigned to a sub constellation

    # dictionary to store the launch dates for each sub constellation
    sub_constellation_launch_dates = {sub_constellation['_soname']: [] for sub_constellation in sub_constellation_metadata_dicts}

    # loop until all launches have been assigned to a sub constellation
    while len(launch_bank) > 0 and len(sub_constellation_launches_required) > 0:
        # find the sub constellation with the most launches required
        max_launches_sub_constellation = max(sub_constellation_launches_required, key=sub_constellation_launches_required.get)

        # assign the next launch in the launch bank to the sub constellation with the most launches required
        launch_date = launch_bank.pop(0)
        sub_constellation_launch_dates[max_launches_sub_constellation].append(launch_date.strftime("%Y-%m-%d"))

        # decrement the number of launches required for the selected sub constellation
        sub_constellation_launches_required[max_launches_sub_constellation] -= 1

        # remove the sub constellation from sub_constellation_launches_required if all its launches have been assigned
        if sub_constellation_launches_required[max_launches_sub_constellation] == 0:
            sub_constellation_launches_required.pop(max_launches_sub_constellation)

    # return the sub constellation launch dates
    return sub_constellation_launch_dates

def create_subconstellation_Space_Objects(N, i, h, _soname, _application, _owner, launch_schedule, _mass, _area, _length, _maneuverable, _propulsion):
    
    #check that N is always >= 1 
    if N < 1:
        raise ValueError('N (number of satellites in constellation) must be >= 1')
    #create a list of SpaceObject instances for each satellite in the a subconstellation metadata dictionary
    subconstellation_Space_Objects = []
    apogee_alt = h  # assuming circular orbit
    perigee_alt = h  # assuming circular orbit
    launch_site = "JSC"

    soname = _soname
    sotype = "PAY" # As we are only dealing with satellites this is the only option (https://celestrak.org/satcat/satcat-format.php)
    operationalstatus = "+" # assuming they are all operational on launch this should also be the only option (https://celestrak.org/satcat/status.php)
    orbit_type = orbit_classify(h) #calculated based on altitude
    application = _application #this comes directly from the CSV file
    owner = _owner  #this comes directly from the CSV file
    orbit_info_status = "EA0" #All Earth orbiting so all "EA0"

    mass = _mass # this comes directly from the CSV file
    maneuverable = 'y' #this comes directly from the CSV file
    spin_stabilized = 'n' #setting to always no ->  all LEO so no spin stabilization
    area = _area  #computed based on form factor
    length = _length  #computed based on form factor
    propulsion = _propulsion  #this comes directly from the CSV file

    R_earth = 6371  # Earth's mean radius in km
    a = R_earth + h #semi-major axis in km
    ###### TODO: the code below assumed a certain orbital geometry. 
    ###### This needs to be unhardcoded to allow us to specify different constellation geometries that will will result in different keplerian elements
    ###### This needs to go into a function
    P = math.ceil(math.sqrt(N)) #number of planes
    S = math.ceil(N / P) #number of satellites per plane
    delta_Omega = 360 / P #angle between planes
    delta_omega = 360 / S #angle between satellites in the same plane -> argument of perigee
    #now calculate the true anomaly(nu) of each satellite in the plane
    delta_TRAN = 360 / N #angle between satellites in the same plane -> true anomaly
    period = orbital_period(a) #orbital period in minutes
    #if any of the above are nan stop and return an error
    if math.isnan(a) or math.isnan(P) or math.isnan(S) or math.isnan(delta_Omega) or math.isnan(delta_omega) or math.isnan(delta_TRAN) or math.isnan(period):
        print("ERROR: NaN value in FSP launch entry")
        return

    launch_entry = []
    launch_schedule = launch_schedule * (N // len(launch_schedule)) + launch_schedule[:N % len(launch_schedule)] #This will repeat the launch_schedule list enough times to cover all satellites and then take the remainder if there's still a difference.
    # TODO: fix this so that it does repeat the launch schedule in this way but instead does it sequentially
    #take each instance in the launch schedule, make it into a datetime obejct, increment the dates by 5 years and add it to the decay schedule
    decay_schedule = [(datetime.strptime(date, "%Y-%m-%d") + timedelta(days=365*5)).strftime("%Y-%m-01") for date in launch_schedule]
    for n in range(1, N+1, 1): # Start iteration at 1, this way, the calculation for Omega_n won't be negative because n-1 will always be positive.
        a_n = a #semi-major axis in km
        e_n = 0 #eccentricity
        i_n = i #inclination in degrees
        omega_n = (n % S) * delta_omega #argument of perigee in degrees
        Omega_n = math.floor((n-1) / S) * delta_Omega #longitude/right ascension of ascending node in degrees
        TRAN_n = (n-1) * delta_TRAN #true anomaly in degrees
        #if tran is nan then set it to 0
        if math.isnan(TRAN_n):
            TRAN_n = 0
        
        launch_year = int(launch_schedule[n-1][:4]) # index with n-1 because the launch_schedule list is indexed from 0
        launch_number = int(n)
        launch_piece = "A" # TODO: unhardcode this -> decide on a way to determine this
        object_cospar_id = generate_cospar_id(launch_year, launch_number, launch_piece)
        subconstellation_Space_Objects.append(SpaceObject(object_type=sotype, payload_operational_status=operationalstatus, 
                                                          application=application, operator=owner, 
                                                          apogee=apogee_alt, perigee=perigee_alt,
                                                          mass=mass, maneuverable=maneuverable, spin_stabilized=spin_stabilized, 
                                                        characteristic_area=area, characteristic_length=length, 
                                                          propulsion_type=propulsion, sma=a_n, eccentricity=e_n, inc=np.deg2rad(i_n), argp=np.deg2rad(omega_n), 
                                                          raan=np.deg2rad(Omega_n), tran=np.deg2rad(TRAN_n), launch_site=launch_site, launch_date=launch_schedule[n-1], 
                                                          decay_date=decay_schedule[n-1] ,rso_name=soname))

    return subconstellation_Space_Objects

def remove_items_from_dict(dict, failure_rate):
    #TODO: this needs work. leaving it for now as just want to run sim
    # calculate the number of items to remove
    companies_to_exclude = ["Starlink", "OneWeb", "E-Space"]
    output_company_list = []
    companies_to_remove = []

    # firstly remove all the constellations we want to exclude from above list
    for company in dict:
        if company["_owner"] in companies_to_exclude:
            output_company_list.append(company)
        else:
            companies_to_remove.append(company)

    # remove the total number of companies already excluded from the total number of companies
    num_items_to_remove = round(int(len(dict) - len(output_company_list)) * failure_rate) # this is incorrect use of rate. need to fix

    # randomly select items to remove
    indices_to_remove = random.sample(range(len(companies_to_remove)), num_items_to_remove)
    remaining_companies = [v for i, v in enumerate(companies_to_remove) if i not in indices_to_remove]
    
    # add start ups to list of companies excluded
    output_company_list.extend(remaining_companies)
    
    return output_company_list

def apply_settings_at_organisation_level(metadata_dicts, failure_rate , remove_operators):
    # Firstly, remove the startups
    metadata_dicts_settings_applied = remove_items_from_dict(metadata_dicts, failure_rate = failure_rate)

    # Ability to remove specific operators
    if remove_operators:
        operators_to_remove = [operator.strip() for operator in remove_operators.split(",")]
        metadata_dicts_settings_applied = [dictionary for dictionary in metadata_dicts_settings_applied if dictionary["_owner"] not in operators_to_remove]
    else:
        metadata_dicts_settings_applied = metadata_dicts_settings_applied

    return metadata_dicts_settings_applied

# def apply_settings_at_constellation_level(metadata_dicts, settings):
#     list_of_dicts_settings_applied = []
#     for metadata_dict in metadata_dicts:
#         num_launches = np.maximum(np.random.normal(float(settings['satellite_failure_rate']), 0.25, 10)[0], 1) # ensure not negative
#         metadata_dict["N"] = round(num_launches * metadata_dict["N"])
#         list_of_dicts_settings_applied.append(metadata_dict)

#     return list_of_dicts_settings_applied

def Prediction2SpaceObjects(satellite_predictions_csv, simsettings):
    """Generate instances of the SpaceObject class for each of the satellties in the prediction data excel file

    Args:
        in_csv_path (_type_): CSV file with the satellite predictions data
        simsettings_json (_type_): loaded JSON file with the simulation settings to be applied (launch model parameters contained here)
    """
    operators_to_remove = simsettings["remove_operators"]
    all_space_objects = []
    # create a list of dictionaries containing the metadata for each sub-constellation
    metadata_dicts = satellite_metadata(file_path=satellite_predictions_csv)

    # need to be able to policy to the constellation companies, then at constellation level (i.e number of satellites that have failed)
    metadata_dicts = apply_settings_at_organisation_level(metadata_dicts,failure_rate = 0, remove_operators =  operators_to_remove)
    # metadata_dicts = apply_settings_at_constellation_level(metadata_dicts, policy)

    #global_launch_schedule(sub_constellation_metadata_dicts=metadata_dicts)
    sub_constellation_launch_dates = global_launch_schedule(
                                                            sub_constellation_metadata_dicts = metadata_dicts, 
                                                            settings=simsettings,
                                                            monthly_ton_capacity=float(simsettings['monthly_ton_capacity']), 
                                                            launches_start_date = simsettings['launch_start_date']   
                                                         )                                                     
    
    for dict in metadata_dicts: #TODO: use the JSON file to set the agressivity and max launch rate etc.
        subconstellation_Space_Objects = create_subconstellation_Space_Objects(N=int(dict['N']), i = float(dict['i']), h=float(dict['h']), _soname=dict['_soname'], _application = dict['_application'], _owner= dict['_owner'], launch_schedule = sub_constellation_launch_dates[dict['_soname']], _mass=dict['_mass'], _area=dict['_area'], _length=dict['_length'], _maneuverable= dict['_maneuverable'], _propulsion=dict['_propulsion'])
        all_space_objects.extend(subconstellation_Space_Objects)
    return all_space_objects

if __name__ == '__main__':
    def get_path(*args):
        return os.path.join(os.getcwd(), *args)

    in_file = 'src/data/prediction_csv/04_04_23_fsp.csv'
    settings = json.load(open(get_path('src/data/specify_simulation/testsim.json'), 'r')) 
    Prediction2SpaceObjects(in_file, settings)