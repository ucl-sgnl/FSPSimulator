import datetime
import sys
import os
import json
import pickle
from tqdm import tqdm
from multiprocessing import Pool, cpu_count

from utils.SpaceCatalogue import SpaceCatalogue
from utils.Conversions import utc_to_jd

def get_path(*args):
    return os.path.join(os.getcwd(), *args)

def load_pickle(file_path):
    with open(get_path(file_path), 'rb') as f:
        return pickle.load(f)

def dump_pickle(file_path, data):
    with open(get_path(file_path), 'wb') as f:
        pickle.dump(data, f)

def propagate_satellite(args):

    satellite, jd_start, jd_stop, step_size = args
    # Execute the prop_catobject method
    satellite.prop_catobject(jd_start=jd_start, jd_stop=jd_stop, step_size=step_size, propagator="RK45")

    return satellite

def run_parallel_sim(settings):
    # This function will define the simulation settings and run the simulation based on the SpaceCatalogue class instantiation.
    # The SpaceCatalogue class will create a catalogue(list) of SpaceObjects based on the settings provided.
    # This list of SpaceObjects is just metadata (empty ephemerides) until we propagate the SpaceObjects using the prop_catobject method.

    SATCAT = SpaceCatalogue(settings["sim_object_type"], settings["sim_object_catalogue"], settings["repull_catalogues"])
    jd_start = float(utc_to_jd(settings["sim_start_date"])[0])
    jd_stop = float(utc_to_jd(settings["sim_end_date"])[0])
    step_size = int(settings["integrator_step_size"]) # in seconds
    scenario_name = settings["scenario_name"]

    print("Number of Satellites in catalogue specified: ", len(SATCAT.Catalogue))
    print(f"Propagating Satellites and saving state vector every {settings['output_frequency']} days...")

    decayed_before_start = 0
    for satellite in SATCAT.Catalogue:
        if satellite.decay_date < datetime.datetime.strptime(settings["sim_start_date"], '%Y-%m-%d'):
            SATCAT.Catalogue.remove(satellite)
            decayed_before_start += 1
    print("# sats decayed before sim start date: ", decayed_before_start)

    #TODO: testing
    #slice SATCAT.Catalogue to retain only the first and last 100 satellites for testing (first 100 are from JSR/SpaceTrack, last 100 are from FSP predictions)
    SATCAT.Catalogue = SATCAT.Catalogue[:100]

    print("Propagating satellites in parallel...")

    iterable = [(satellite, jd_start, jd_stop, step_size) for satellite in SATCAT.Catalogue]
    with Pool(processes=cpu_count()) as pool:
        with tqdm(total=len(iterable)) as pbar:
            results = []
            for result in pool.imap_unordered(propagate_satellite, iterable):
                results.append(result)
                pbar.update()

    SATCAT.Catalogue = results
    print("Exporting results...")
    dump_pickle(f'src/data/catalogue/prop_{scenario_name}.pickle', SATCAT)

    print(f"Output: {get_path(f'src/data/catalogue/{scenario_name}.pickle')}")
    print(f"Number of Satellites: {len(SATCAT.Catalogue)}")
    print("Simulation Complete")

if __name__ == '__main__':
    ##### Profiling #####
    settings = json.load(open(get_path('src/data/prediction_csv/sim_settings.json'), 'r'))
    run_parallel_sim(settings)