import datetime
import sys
import os
import json
import pickle
from tqdm import tqdm
from tqdm.contrib.concurrent import process_map 
from multiprocessing import Pool, cpu_count
import cProfile

from utils.SpaceCatalogue import SpaceCatalogue
from utils.LaunchModel import Prediction2SpaceObjects
from utils.Conversions import utc_to_jd

def get_path(*args):
    return os.path.join(os.getcwd(), *args)

def load_pickle(file_path):
    with open(get_path(file_path), 'rb') as f:
        return pickle.load(f)

def dump_pickle(file_path, data):
    with open(get_path(file_path), 'wb') as f:
        pickle.dump(data, f)

def propagate_satellite_profiling(args):
    satellite, jd_start, jd_stop = args

    # Create a cProfile.Profile object and enable it to start profiling
    profiler = cProfile.Profile()
    profiler.enable()

    # Execute the prop_catobject method
    satellite.prop_catobject(jd_start=jd_start, jd_stop=jd_stop, step_size=20, propagator="RK45")

    # Disable the profiler
    profiler.disable()

    # Create the directory if it doesn't exist
    directory = 'src/tests/profiling_results/'
    os.makedirs(directory, exist_ok=True)

    # Create the file path
    file_path = os.path.join(directory, f'profiling_{satellite.rso_name}.txt')

    # Save the stats of the profiler into a unique file for each satellite.
    with open(file_path, 'w') as f:
        sys.stdout = f  # redirect output to file
        profiler.print_stats(sort='time')
        sys.stdout = sys.__stdout__  # reset output to normal

    return satellite

def run_parallel_sim(settings):
    SATCAT = SpaceCatalogue(settings["sim_object_type"], settings["sim_object_catalogue"], settings["repull_catalogues"])
    jd_start = float(utc_to_jd(settings["sim_start_date"])[0])
    jd_stop = float(utc_to_jd(settings["sim_end_date"])[0])
    scenario_name = settings["scenario_name"]

    if settings["repull_catalogues"] and os.path.exists(get_path('src/data/catalogue/All_catalogue_latest.csv')):
        if settings["sim_object_type"] == "all":
            SATCAT.CreateCatalogueAll()
        elif settings["sim_object_type"] != "debris":  
            SATCAT.CreateCatalogueActive()

        SATCAT.Catalogue2SpaceObjects()
        
        dump_pickle('src/data/catalogue/SATCAT_before_prop.pickle', SATCAT)

    else:
        SATCAT.Catalogue = load_pickle('src/data/catalogue/SATCAT_before_prop.pickle')

    if settings["environment"] == "development" and os.path.exists(get_path(f'src/data/catalogue/SATCAT_before_prop_{scenario_name}.pickle')):
        SATCAT = load_pickle(f'src/data/catalogue/SATCAT_before_prop_{scenario_name}.pickle')
    else:
        if settings["scenario_name"] != "baseline":
            print("Creating Launch Model...")
            launch_file_object = Prediction2SpaceObjects('src/data/prediction_csv/FSP_Predictions.csv', 'src/data/prediction_csv/sim_settings.json')
            SATCAT.Catalogue.extend(launch_file_object)

    print("Number of Satellites: ", len(SATCAT.Catalogue))
    print(f"Propagating Satellites and saving state every {settings['output_frequency']} days...")

    decayed_before_start = 0
    for satellite in SATCAT.Catalogue:
        if satellite.decay_date < datetime.datetime.strptime(settings["sim_start_date"], '%Y-%m-%d'):
            SATCAT.Catalogue.remove(satellite)
            decayed_before_start += 1
    print("# sats decayed before sim start date: ", decayed_before_start)

    # Propagate satellites in parallel
    print("Propagating satellites in parallel...")

    #slice SATCAT.Catalogue to retain only 100 satellites for testing
    SATCAT.Catalogue = SATCAT.Catalogue[:100]

    chunksize = len(SATCAT.Catalogue) // cpu_count()
    SATCAT.Catalogue = process_map(propagate_satellite_profiling, [(satellite, jd_start, jd_stop) for satellite in SATCAT.Catalogue], max_workers=cpu_count(), chunksize=chunksize)

    print("Exporting results...")
    dump_pickle(f'src/data/catalogue/prop_{scenario_name}.pickle', SATCAT)

    print(f"Output: {get_path(f'src/data/catalogue/{scenario_name}.pickle')}")
    print(f"Number of Satellites: {len(SATCAT.Catalogue)}")
    print("Simulation Complete")

if __name__ == '__main__':
    ##### Profiling #####
    with open('src/tests/profiling_results/profiling_main.txt', 'w') as f:
        sys.stdout = f  # redirect output to file
        policy = json.load(open(get_path('src/data/prediction_csv/sim_settings.json'), 'r'))
        profiler = cProfile.Profile()
        profiler.enable()
        run_parallel_sim(policy)
        profiler.disable()
        profiler.print_stats(sort='time')
    sys.stdout = sys.__stdout__  # reset output to normal