import datetime
import sys
import os
import json
import pickle
from tqdm import tqdm
from multiprocessing import Pool, cpu_count, Lock
from concurrent.futures import ThreadPoolExecutor
from utils.SpaceCatalogue import SpaceCatalogue, check_json_file
from utils.Conversions import utc_to_jd, initialize_orekit

def get_path(*args):
    return os.path.join(os.getcwd(), *args)

def load_pickle(file_path):
    with open(get_path(file_path), 'rb') as f:
        return pickle.load(f)

def dump_pickle(file_path, data):
    with open(get_path(file_path), 'wb') as f:
        pickle.dump(data, f)

def propagate_space_object(args):
    space_object, jd_start, jd_stop, step_size, output_freq, integrator_type, force_model, long_term_sgp4 = args
    # Execute the prop_catobject method on the space object
    print(f"Propagating {space_object.rso_name}...")
    try:
        if long_term_sgp4:
            initialize_orekit()
        space_object.prop_catobject(jd_start=jd_start, jd_stop=jd_stop, step_size=step_size, output_freq=output_freq, integrator_type=integrator_type, force_model=force_model, long_term_sgp4=long_term_sgp4)
    except Exception as e:
        print(f"An error occurred while propagating {space_object.rso_name}: {e}")
        import traceback
        traceback.print_exc()

    return space_object

from concurrent.futures import ProcessPoolExecutor

def run_parallel_sim(settings):
    SATCAT = SpaceCatalogue(settings=settings)
    jd_start = float(utc_to_jd(settings["sim_start_date"])[0])
    jd_stop = float(utc_to_jd(settings["sim_end_date"])[0])
    step_size = int(settings["integrator_step_size"])
    output_freq = int(settings["output_frequency"])
    scenario_name = str(settings["scenario_name"])
    integrator_type = str(settings["integrator_type"])
    sgp4_long_term = bool(settings["sgp4_long_term"])
    force_model = settings["force_model"]

    # Initialize for the main process, but also will be initialized within each child process.
    initialize_orekit()

    print("Number of space_object in catalogue specified:", len(SATCAT.Catalogue))
    print(f"Propagating SpaceObjects and saving state vectors every {settings['output_frequency']} seconds...")

    decayed_before_start = 0
    for space_object in SATCAT.Catalogue:
        if space_object.decay_date < datetime.datetime.strptime(settings["sim_start_date"], '%Y-%m-%d'):
            SATCAT.Catalogue.remove(space_object)
            decayed_before_start += 1
    print("# sats decayed before sim start date:", decayed_before_start)

    # Slice SATCAT.Catalogue to select every 1000th space object (for testing)
    SATCAT.Catalogue = SATCAT.Catalogue[::1000]

    print("Propagating space objects in parallel...")

    # Create iterable
    iterable = [(space_object, jd_start, jd_stop, step_size, output_freq, integrator_type, force_model, sgp4_long_term) for space_object in SATCAT.Catalogue]

    results = []
    with tqdm(total=len(iterable)) as pbar:
        def callback(result):
            results.append(result)
            pbar.update()

        # Use ProcessPoolExecutor for multiprocessing
        with ProcessPoolExecutor(max_workers=os.cpu_count) as executor:
            # Since the map function blocks until all results are available, we can simply use it here.
            for result in executor.map(propagate_space_object, iterable):
                callback(result)

    SATCAT.Catalogue = results
    print("Exporting results...")
    save_path = os.path.dirname(f'src/data/results/propagated_catalogs/{scenario_name}.pickle')
    os.makedirs(save_path, exist_ok=True)
    dump_pickle(f'src/data/results/propagated_catalogs/{scenario_name}.pickle', SATCAT)

    print(f"Simulation complete. Results saved to: {get_path(f'src/data/results/propagated_catalogs/{scenario_name}.pickle')}")
    return

def run_sim(settings):
    SATCAT = SpaceCatalogue(settings=settings)
    jd_start = float(utc_to_jd(settings["sim_start_date"])[0])
    jd_stop = float(utc_to_jd(settings["sim_end_date"])[0])
    step_size = int(settings["integrator_step_size"])
    output_freq = int(settings["output_frequency"])
    scenario_name = str(settings["scenario_name"])
    integrator_type = str(settings["integrator_type"])
    sgp4_long_term = bool(settings["sgp4_long_term"])
    force_model = settings["force_model"]
    initialize_orekit()

    print("Number of space objects in catalogue specified:", len(SATCAT.Catalogue))

    # Filter out space objects that decayed before the simulation start date
    SATCAT.Catalogue = [space_object for space_object in SATCAT.Catalogue if space_object.decay_date >= datetime.datetime.strptime(settings["sim_start_date"], '%Y-%m-%d')]

    print("# sats decayed before sim start date:", len(SATCAT.Catalogue))

    # Slice SATCAT.Catalogue to select every 1000th space object (for testing)
    SATCAT.Catalogue = SATCAT.Catalogue[::1000]

    print("Propagating space objects...")

    results = []
    
    for space_object in tqdm(SATCAT.Catalogue):
        lock = Lock() # Create a lock object
        with lock:
            result = propagate_space_object((space_object, jd_start, jd_stop, step_size, output_freq, integrator_type, force_model, sgp4_long_term))
            results.append(result)

    SATCAT.Catalogue = results
    print("Exporting results...")
    save_path = os.path.dirname(f'src/data/results/propagated_catalogs/{scenario_name}.pickle')
    os.makedirs(save_path, exist_ok=True)
    dump_pickle(f'src/data/results/propagated_catalogs/{scenario_name}.pickle', SATCAT)

    print(f"Simulation complete. Results saved to: {get_path(f'src/data/results/propagated_catalogs/{scenario_name}.pickle')}")

if __name__ == '__main__':
    #list all the json files in src/data/specify_simulations
    sims = os.listdir(get_path('src/data/specify_simulations/'))
    for sim in sims:
        if sim.endswith('.json'):
            print(f"Running simulation: {sim}")
            settings = json.load(open(get_path(f'src/data/specify_simulations/{sim}'), 'r'))
            check_json_file(settings)#check if the json file is filled out correctly
            run_parallel_sim(settings)
            print(f"Simulation {sim} complete")