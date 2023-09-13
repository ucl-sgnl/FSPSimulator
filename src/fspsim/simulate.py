import os
import json
import pickle
from tqdm import tqdm
import traceback
from .utils.SpaceCatalogue import SpaceCatalogue, check_json_file
from .utils.Conversions import utc_to_jd
from .utils.Formatting import future_constellations_csv_handler

future_constellations_dict = {}

def get_path(*args):
    return os.path.join(os.getcwd(), *args)

def load_pickle(file_path):
    with open(get_path(file_path), 'rb') as f:
        return pickle.load(f)

def dump_pickle(file_path, data):
    with open(get_path(file_path), 'wb') as f:
        pickle.dump(data, f)

def propagate_space_object(args):
    space_object, jd_start, jd_stop, step_size, output_freq = args
    # Execute the prop_catobject method on the space object
    #print(f"Propagating {space_object.rso_name}...")
    try:
        space_object.prop_catobject(jd_start=jd_start, jd_stop=jd_stop, step_size=step_size, output_freq=output_freq)
    except Exception as e:
        print(f"An error occurred while propagating {space_object.rso_name}: {e}")
        #print the error location
        print(traceback.format_exc())
    return

def set_future_constellations(constellations) -> bool: 
    """Incorporates user-specified launch predictions into the simulation. 
    This will check that it is in the correct format
    and can be used by the simulation.

    :param constellations: The path of the specified csv file
    :type constellations: str
    :return: If it is in the correct format.
    :rtype: bool
    """
    future_constellations_dict = future_constellations_csv_handler(constellations)

    if not future_constellations_dict:
        print("There are no future constellations in your provided dictionary")
        return False
    
    return True

def run_sim(settings: json, future_constellations_file: str = None) -> None:
    """
    Propagates a list of space objects over a specified time range.

    The function retrieves space objects from a space catalogue (SATCAT), propagates each object
    up to a the specified time. This will save the output locally but will also return a list of propagated
    space ojects. A progress bar is displayed to track the propagation process.

    Parameters:
        - settings (dict): A dictionary containing the following key-value pairs:
        - "sim_start_date": UTC start date for simulation (str)
        - "sim_end_date": UTC end date for simulation (str)
        - "integrator_step_size": Time step size for integrator (str, converted to int)
        - "output_frequency": Output frequency (str, converted to int)
        - "scenario_name": Name of the scenario (str)
        - "integrator_type": Type of integrator to use (str)
        - "sgp4_long_term": Boolean indicating if SGP4 long term propagation should be used (str, converted to bool)
        - "force_model": Force model settings (can be various types, depending on implementation)

    Returns:
    SATCAT Catalogue with updated ephemeris of the locations.
    Results are saved to pickle files in the directory: 'src/fspsim/data/results/propagated_catalogs/'.

    Notes:
    - Each batch consists of a set number of space objects (defined by the batch_size variable).
    - The pickle files are saved with a filename pattern: '<scenario_name>_batch_<batch_number>.pickle'.
    - The propagation function used is 'propagate_space_object' (not defined in the provided code snippet).
    """
    if future_constellations_file:
        # User defined fsp launch files
        # check if valid 
        if future_constellations_csv_handler(future_constellations_file):
            SATCAT = SpaceCatalogue(settings=settings, future_constellations=future_constellations_file)
        else:
            return
    else:
        # Use default fsp launch files
        SATCAT = SpaceCatalogue(settings=settings)
    jd_start = float(utc_to_jd(settings["sim_start_date"]))
    jd_stop = float(utc_to_jd(settings["sim_end_date"]))
    step_size = int(settings["integrator_step_size"])
    output_freq = int(settings["output_frequency"])
    scenario_name = str(settings["scenario_name"])

    # only select every 1000 objects in SATCAT.Catalogue
    # SATCAT.Catalogue = SATCAT.Catalogue[::1000] 

    # Create a progress bar
    pbar = tqdm(total=len(SATCAT.Catalogue), desc="Propagating")
    propagated_objects = []

    while SATCAT.Catalogue:
        # Propagate space objects
        if not SATCAT.Catalogue:
            break

        space_object = SATCAT.Catalogue.pop(0)
        propagate_space_object((space_object, jd_start, jd_stop, step_size, output_freq))
        propagated_objects.append(space_object)
        pbar.update(1)  # Update the progress bar

        save_path = os.path.join(f'src/fspsim/data/results/propagated_catalogs/', f"{scenario_name}.pickle")
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        dump_pickle(save_path, propagated_objects)

    pbar.close() # Close the progress bar
    print(f"Simulation complete. Results saved to batches in: {get_path(f'src/fspsim/data/results/propagated_catalogs/')}")
    return SATCAT

if __name__ == '__main__':
    sims = os.listdir(get_path('src/fspsim/data/specify_simulations/'))
    for sim in sims:
        if sim.endswith('.json'):
            print(f"Running simulation: {sim}")
            settings = json.load(open(get_path(f'src/fspsim/data/specify_simulations/{sim}'), 'r'))
            check_json_file(settings)#check if the json file is filled out correctly
            # provide my own launch file # charles you might have to change this path
            future_constellations_file = r'C:\Users\IT\Documents\UCL\FSPSimulator\src\fspsim\data\prediction_csv\oneweb_starlink.csv'
            run_sim(settings, future_constellations_file)
            print(f"Simulation {sim} complete")