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

def run_sim(settings):
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
            launch_file_object = Prediction2SpaceObjects('src/data/prediction_csv/_full.csv', 'src/data/prediction_csv/sim_settings.json')
            SATCAT.Catalogue.extend(launch_file_object)

    print("Number of Satellites: ", len(SATCAT.Catalogue))
    # Propagate
    output_frequency = int(settings["output_frequency"])
    timestep = output_frequency * 86400
    print(f"Propagating Satellites and saving state every {output_frequency} days...")

    # Filter satellites based on decay_date
    decayed_before_start = 0
    for satellite in SATCAT.Catalogue:  
        # print(satellite.rso_name)
        if satellite.decay_date < datetime.datetime.strptime(settings["sim_start_date"], '%Y-%m-%d'): # if we know that the decay date is before the start of the simulation, we can remove it from the catalogue
            SATCAT.Catalogue.remove(satellite)
            decayed_before_start += 1
    print("# sats decayed before sim start date: ", decayed_before_start)

    for satellite in tqdm(SATCAT.Catalogue): #tqdm is a progress bar
        satellite.prop_catobject(jd_start = jd_start, jd_stop = jd_stop, step_size = 10, propagator="RK45") # propagate each satellite using sgp4, and output a state vector every timestep
    # Export
    print("Exporting results...")
    dump_pickle(f'src/data/catalogue/prop_{scenario_name}.pickle', SATCAT)

    print(f"Output: {get_path(f'src/data/catalogue/{scenario_name}.pickle')}")
    print(f"Number of Satellites: {len(SATCAT.Catalogue)}")
    print("Simulation Complete")

def propagate_satellite(args):

    satellite, jd_start, jd_stop = args
    # Execute the prop_catobject method
    satellite.prop_catobject(jd_start=jd_start, jd_stop=jd_stop, step_size=20, propagator="RK45")

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
            launch_file_object = Prediction2SpaceObjects('src/data/prediction_csv/FSP_Predictions_full.csv', 'src/data/prediction_csv/sim_settings.json')
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

    #slice SATCAT.Catalogue to retain only 1000 satellites for testing
    SATCAT.Catalogue = SATCAT.Catalogue[:100]

    chunksize = len(SATCAT.Catalogue) // cpu_count()
    SATCAT.Catalogue = process_map(propagate_satellite, [(satellite, jd_start, jd_stop) for satellite in SATCAT.Catalogue], max_workers=cpu_count(), chunksize=chunksize)

    print("Exporting results...")
    dump_pickle(f'src/data/catalogue/prop_{scenario_name}.pickle', SATCAT)

    print(f"Output: {get_path(f'src/data/catalogue/{scenario_name}.pickle')}")
    print(f"Number of Satellites: {len(SATCAT.Catalogue)}")
    print("Simulation Complete")

if __name__ == '__main__':
    ##### Profiling #####
    policy = json.load(open(get_path('src/data/prediction_csv/sim_settings.json'), 'r'))
    run_parallel_sim(policy)