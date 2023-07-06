import datetime
import sys
import os
import json
import pickle
from tqdm import tqdm
from multiprocessing import Pool, cpu_count

from utils.SpaceCatalogue import SpaceCatalogue, check_json_file
from utils.Conversions import utc_to_jd

def get_path(*args):
    return os.path.join(os.getcwd(), *args)

def load_pickle(file_path):
    with open(get_path(file_path), 'rb') as f:
        return pickle.load(f)

def dump_pickle(file_path, data):
    with open(get_path(file_path), 'wb') as f:
        pickle.dump(data, f)

def propagate_space_object(args):

    space_object, jd_start, jd_stop, step_size, output_freq, integrator_type= args
    # Execute the prop_catobject method on the space object
    space_object.prop_catobject(jd_start=jd_start, jd_stop=jd_stop, step_size=step_size, output_freq= output_freq, integrator_type="RK45")

    return space_object

def run_parallel_sim(settings):
    # This function will define the simulation settings and run the simulation based on the SpaceCatalogue class instantiation.
    # The SpaceCatalogue class will create a catalogue(list) of SpaceObjects based on the settings provided.
    # This list of SpaceObjects is just metadata (empty ephemerides) until we propagate the SpaceObjects using the prop_catobject method.

    SATCAT = SpaceCatalogue(settings = settings)
    jd_start = float(utc_to_jd(settings["sim_start_date"])[0])
    jd_stop = float(utc_to_jd(settings["sim_end_date"])[0])
    step_size = int(settings["integrator_step_size"]) # in seconds
    output_freq = int(settings["output_frequency"]) # in seconds
    scenario_name = settings["scenario_name"] # this will be used to name the output pickle file
    integrator_type = settings["integrator_type"]  # must be one of "RK45", "RK23", "DOP853", "Radau", "BDF", "LSODA"

    print("Number of space_object in catalogue specified: ", len(SATCAT.Catalogue))
    print(f"Propagating SpaceObjects and saving state vectors every {settings['output_frequency']} seconds...")

    decayed_before_start = 0
    for space_object in SATCAT.Catalogue:
        if space_object.decay_date < datetime.datetime.strptime(settings["sim_start_date"], '%Y-%m-%d'):
            SATCAT.Catalogue.remove(space_object)
            decayed_before_start += 1
    print("# sats decayed before sim start date: ", decayed_before_start)

    #TODO: testing
    #slice SATCAT.Catalogue to select n elements evenly spaced throughout the list
    SATCAT.Catalogue = SATCAT.Catalogue[::int(len(SATCAT.Catalogue)/10)]

    print("Propagating space objects in parallel...")

    iterable = [(space_object, jd_start, jd_stop, step_size, output_freq, integrator_type) for space_object in SATCAT.Catalogue]
    with Pool(processes=cpu_count()) as pool:
        with tqdm(total=len(iterable)) as pbar:
            results = []
            for result in pool.imap_unordered(propagate_space_object, iterable):
                results.append(result)
                pbar.update()

    SATCAT.Catalogue = results
    print("Exporting results...")
    dump_pickle(f'src/data/results/propagated_catalogs/{scenario_name}.pickle', SATCAT)
    print(f"Simulation complete. Results saved to: {get_path(f'src/data/results/propagated_catalogs/{scenario_name}.pickle')}")

if __name__ == '__main__':
    #list all the json files in src/data/specify_simulations
    sims = os.listdir(get_path('src/data/specify_simulations'))
    for sim in sims:
        if sim.endswith('.json'):
            print(f"Running simulation: {sim}")
            settings = json.load(open(get_path(f'src/data/specify_simulations/{sim}'), 'r'))
            check_json_file(settings)#check if the json file is filled out correctly
            run_parallel_sim(settings)
            print(f"Simulation {sim} complete")