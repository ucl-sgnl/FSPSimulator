import datetime
import sys
import os
import json
import pickle
from tqdm import tqdm
from multiprocessing import Pool, cpu_count
from fspsim.utils.SpaceCatalogue import SpaceCatalogue, check_json_file
from fspsim.utils.Conversions import utc_to_jd
import logging
import importlib.resources
import pandas as pd
from dotenv import set_key

def get_path(*args):
    return os.path.join(os.getcwd(), *args)

def load_pickle(file_path):
    with open(get_path(file_path), 'rb') as f:
        return pickle.load(f)

def dump_pickle(file_path, data):
    with open(get_path(file_path), 'wb') as f:
        pickle.dump(data, f)

def set_spacetrack_login(username, password):
    set_key(".env", "SPACETRACK_USERNAME", username)
    set_key(".env", "SPACETRACK_PASSWORD", password)

def propagate_space_object(args):
    space_object, jd_start, jd_stop, step_size, output_freq, integrator_type= args
    # Execute the prop_catobject method on the space object
    space_object.prop_catobject(jd_start=jd_start, jd_stop=jd_stop, step_size=step_size, output_freq= output_freq, integrator_type="RK45")
    return space_object

def return_baseline_catalogue():
    """
    Returns a pickle file of the current catalogue of space objects
    """
    with importlib.resources.open_binary('fspsim.data.catalogue', 'SATCAT_before_prop.pickle') as resource:
        print(resource)
        return pickle.loads(resource)
    
def return_baseline_csv():
    with importlib.resources.open_text('fspsim.data.catalogue', 'Cleaned_All_catalogue_latest.csv') as file:
        df = pd.read_csv(file)
    return df
 
def run_parallel_sim(settings):
    """
    The main entry point into the FSP Simulator. This function will run the simulation based on the settings provided in the json file.
    This function will define the simulation settings and run the simulation based on the SpaceCatalogue class instantiation.
    The SpaceCatalogue class will create a catalogue(list) of SpaceObjects based on the settings provided.
    This list of SpaceObjects is just metadata (empty ephemerides) until we propagate the SpaceObjects using the prop_catobject method.   
    ### Exports
    - List of merged space objects to 'src/fspsim/data/external/active_jsr_spacetrack.csv'
    """
    # first check the json file
    if check_json_file(settings) == False:
        sys.exit()

    SATCAT = SpaceCatalogue(settings = settings)
    jd_start = float(utc_to_jd(settings["sim_start_date"])[0])
    jd_stop = float(utc_to_jd(settings["sim_end_date"])[0])
    step_size = int(settings["integrator_step_size"]) # in seconds
    output_freq = int(settings["output_frequency"]) # in seconds
    scenario_name = settings["scenario_name"] # this will be used to name the output pickle file
    integrator_type = settings["integrator_type"]  # must be one of "RK45", "RK23", "DOP853", "Radau", "BDF", "LSODA"

    logging.info("Number of space_object in catalogue specified: ", len(SATCAT.Catalogue))
    logging.info(f"Propagating SpaceObjects and saving state vectors every {settings['output_frequency']} seconds...")

    print("Stopping for now because we need to test the code")
    # decayed_before_start = 0
    # for space_object in SATCAT.Catalogue:
    #     if space_object.decay_date < datetime.datetime.strptime(settings["sim_start_date"], '%Y-%m-%d'):
    #         SATCAT.Catalogue.remove(space_object)
    #         decayed_before_start += 1
    # logging.info("# sats decayed before sim start date: ", decayed_before_start)

    # logging.info("Propagating space objects in parallel...")
    # iterable = [(space_object, jd_start, jd_stop, step_size, output_freq, integrator_type) for space_object in SATCAT.Catalogue]
    # with Pool(processes=cpu_count()) as pool:
    #     with tqdm(total=len(iterable)) as pbar:
    #         results = []
    #         for result in pool.imap_unordered(propagate_space_object, iterable):
    #             results.append(result)
    #             pbar.update()

    # SATCAT.Catalogue = results
    # logging.info("Exporting results...")
    # dump_pickle(f'src/fspsim/data/results/propagated_catalogs/{scenario_name}.pickle', SATCAT)
    # logging.info(f"Simulation complete. Results saved to: {get_path(f'src/fspsim/data/results/propagated_catalogs/{scenario_name}.pickle')}")

if __name__ == '__main__':
    sims = os.listdir(get_path('src/fspsim/data/specify_simulations/'))
    for sim in sims:
        if sim.endswith('.json'):
            print(f"Running simulation: {sim}")
            settings = json.load(open(get_path(f'src/fspsim/data/specify_simulations/{sim}'), 'r'))
            check_json_file(settings)#check if the json file is filled out correctly
            run_parallel_sim(settings)
            print(f"Simulation {sim} complete")