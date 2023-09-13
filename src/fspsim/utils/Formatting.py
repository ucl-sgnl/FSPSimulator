import pandas as pd
import sys
import os 
import pickle

def calculate_form_factor(form_factor_str):
    """
    Convert the form factor specified in the predictions csv to a characteristic length and area.
    For CubeSats, the area can be specified by writing the number of "U"s in the form factor string. e.g. "3U" or "6u" (the case does not matter)
    For other satellites the form factor must be three numbers (floats or ints), separated by stars. e.g. "1.2*2.3*4.5" or "1*2*3"

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
    elif '*' in form_factor_str:
        form_factor_str = sorted([float(i) for i in form_factor_str.split('*')])
        _length = form_factor_str[-1]
        _area = form_factor_str[-1] * form_factor_str[-2]
    else:
        raise ValueError('Form factor must be of the form "3U" or "1*2*3"')

    return _length, _area

def future_constellations_csv_handler(file_path):
    """
    Checks that the user supplied Future Constellation CSV is in the correct format for the simulation

    :param file_path: File Path of the CSV
    :type file_path: str
    :return: Dictionary of the constellations in a format the fspsim can read.
    :rtype: dict
    """
    try:
        sat_df = pd.read_csv(file_path, sep=',')
    except FileNotFoundError:
        print("Future Constellation File has not been found")
        sys.exit()
    except:
        print("Cannot read Future Constellation csv")
        sys.exit() 

    expected_headers = ['Number of sats', 'Inclination', 'Altitude', 'Sub-Constellation', 'Mission type/application', 'Mass(kg)', 'Form Factor', 'Maneuverable','Propulsion']
    
    # if any of the expected headers are missing, stop the simulation
    if not set(expected_headers).issubset(set(sat_df.columns)):
        print("The Future Constellation file is not complete. Please fill any missing values.\n Stopping simulation.")
        sys.exit()
        
    sat_df = sat_df.dropna(subset=expected_headers)

    # Apply functions to entire columns
    sat_df['Sub-Constellation'] = sat_df['Sub-Constellation'].str.replace(' ', '_')
    sat_df[['_length', '_area']] = sat_df['Form Factor'].apply(calculate_form_factor).to_list()

    sat_df = sat_df.rename(columns={
        'Mega-Constellation': '_owner',
        'Number of sats': 'N',
        'Inclination': 'i',
        'Altitude': 'h',
        'Sub-Constellation': '_soname',
        'Mission type/application': '_application',
        'Mass(kg)': '_mass',
        'Maneuverable': '_maneuverable',
        'Propulsion': '_propulsion'
    })
    
    # Convert DataFrame to list of dictionaries
    metadata_dicts = sat_df.to_dict('records')

    return metadata_dicts

def get_path(*args):
    return os.path.join(os.getcwd(), *args)

def load_pickle(file_path):
    with open(get_path(file_path), 'rb') as f:
        return pickle.load(f)

def dump_pickle(file_path, data):
    with open(get_path(file_path), 'wb') as f:
        pickle.dump(data, f)