import datetime
import sys
import os
import json
import pickle
from tqdm import tqdm
from multiprocessing import Pool, cpu_count
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
        space_object.prop_catobject(jd_start=jd_start, jd_stop=jd_stop, step_size=step_size, output_freq=output_freq, integrator_type=integrator_type, force_model=force_model, long_term_sgp4=long_term_sgp4)
    except Exception as e:
        print(f"An error occurred while propagating {space_object.rso_name}: {e}")
        import traceback
        traceback.print_exc()

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

    batch = 1
    batch_size = 100
    while SATCAT.Catalogue:
        current_batch = []

        # Propagate space objects for the current batch
        for _ in range(batch_size):
            if not SATCAT.Catalogue:
                break

            space_object = SATCAT.Catalogue.pop(0)
            propagate_space_object((space_object, jd_start, jd_stop, step_size, output_freq, integrator_type, force_model, sgp4_long_term))
            current_batch.append(space_object)

        # Save the current batch
        print(f"Saving batch {batch}...")
        save_path = os.path.join(f'src/data/results/propagated_catalogs/', f"{scenario_name}_batch_{batch}.pickle")
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        dump_pickle(save_path, current_batch)
        batch += 1

    print(f"Simulation complete. Results saved to batches in: {get_path(f'src/data/results/propagated_catalogs/')}")


if __name__ == '__main__':
    #list all the json files in src/data/specify_simulations
    sims = os.listdir(get_path('src/data/specify_simulations/'))
    for sim in sims:
        if sim.endswith('.json'):
            print(f"Running simulation: {sim}")
            settings = json.load(open(get_path(f'src/data/specify_simulations/{sim}'), 'r'))
            check_json_file(settings)#check if the json file is filled out correctly
            run_sim(settings)
            print(f"Simulation {sim} complete")

###### MULTIPROCESSING ######
# def run_parallel_sim(settings):
#     # This function will define the simulation settings and run the simulation based on the SpaceCatalogue class instantiation.
#     # The SpaceCatalogue class will create a catalogue(list) of SpaceObjects based on the settings provided.
#     # This list of SpaceObjects is just metadata (empty ephemerides) until we propagate the SpaceObjects using the prop_catobject method.

#     SATCAT = SpaceCatalogue(settings = settings)
#     jd_start = float(utc_to_jd(settings["sim_start_date"])[0])
#     jd_stop = float(utc_to_jd(settings["sim_end_date"])[0])
#     step_size = int(settings["integrator_step_size"]) # in seconds
#     output_freq = int(settings["output_frequency"]) # in seconds
#     scenario_name = str(settings["scenario_name"]) # this will be used to name the output pickle file
#     integrator_type = str(settings["integrator_type"])  # must be one of "RK45", "RK23", "DOP853", "Radau", "BDF", "LSODA"
#     sgp4_long_term = bool(settings["sgp4_long_term"]) # boolean, if true, use sgp4 long term propagator
#     force_model = settings["force_model"] # list of strings, can be  # "all" or "grav_mono" and/or "j2"  and/or "sun_grav" and/or "moon_grav" and/or "drag_aero" and/or "srp" 

#     print("Number of space_object in catalogue specified: ", len(SATCAT.Catalogue))
#     print(f"Propagating SpaceObjects and saving state vectors every {settings['output_frequency']} seconds...")

#     decayed_before_start = 0
#     for space_object in SATCAT.Catalogue:
#         if space_object.decay_date < datetime.datetime.strptime(settings["sim_start_date"], '%Y-%m-%d'):
#             SATCAT.Catalogue.remove(space_object)
#             decayed_before_start += 1
#     print("# sats decayed before sim start date: ", decayed_before_start)

#     #TODO: testing
#     #slice SATCAT.Catalogue to select everth 1000th space object
#     SATCAT.Catalogue = SATCAT.Catalogue[::1000]

#     print("Propagating space objects in parallel...")

#     iterable = [(space_object, jd_start, jd_stop, step_size, output_freq, integrator_type, force_model, sgp4_long_term) for space_object in SATCAT.Catalogue]
#     with Pool(initializer=initialize_orekit, processes=cpu_count()) as pool:
#         with tqdm(total=len(iterable)) as pbar:
#             results = []
#             for result in pool.imap_unordered(propagate_space_object, iterable):
#                 results.append(result)
#                 pbar.update()

#     SATCAT.Catalogue = results
#     print("Exporting results...")
#     save_path = os.path.dirname(f'src/data/results/propagated_catalogs/{scenario_name}.pickle')
#     os.makedirs(save_path, exist_ok=True)
#     dump_pickle(f'src/data/results/propagated_catalogs/{scenario_name}.pickle', SATCAT)

#     print(f"Simulation complete. Results saved to: {get_path(f'src/data/results/propagated_catalogs/{scenario_name}.pickle')}")

###### MULTITHREADING ######

# from concurrent.futures import ThreadPoolExecutor

# def run_parallel_sim(settings):
#     SATCAT = SpaceCatalogue(settings=settings)
#     jd_start = float(utc_to_jd(settings["sim_start_date"])[0])
#     jd_stop = float(utc_to_jd(settings["sim_end_date"])[0])
#     step_size = int(settings["integrator_step_size"])
#     output_freq = int(settings["output_frequency"])
#     scenario_name = str(settings["scenario_name"])
#     integrator_type = str(settings["integrator_type"])
#     sgp4_long_term = bool(settings["sgp4_long_term"])
#     force_model = settings["force_model"]
#     initialize_orekit()

#     print("Number of space_object in catalogue specified: ", len(SATCAT.Catalogue))
#     print(f"Propagating SpaceObjects and saving state vectors every {settings['output_frequency']} seconds...")

#     decayed_before_start = 0
#     for space_object in SATCAT.Catalogue:
#         if space_object.decay_date < datetime.datetime.strptime(settings["sim_start_date"], '%Y-%m-%d'):
#             SATCAT.Catalogue.remove(space_object)
#             decayed_before_start += 1
#     print("# sats decayed before sim start date: ", decayed_before_start)

#     # Slice SATCAT.Catalogue to select every 1000th space object (for testing)
#     SATCAT.Catalogue = SATCAT.Catalogue[::100]

#     print("Propagating space objects in parallel...")

#     # Function to update the progress bar
#     def update_result(result):
#         results.append(result)
#         pbar.update()

#     # Create iterable
#     iterable = [(space_object, jd_start, jd_stop, step_size, output_freq, integrator_type, force_model, sgp4_long_term) for space_object in SATCAT.Catalogue]

#     with tqdm(total=len(iterable)) as pbar:
#         results = []
#         # Use ThreadPoolExecutor for multithreading
#         with ThreadPoolExecutor(max_workers=cpu_count()) as executor:
#             # Submit tasks
#             futures = [executor.submit(propagate_space_object, args) for args in iterable]
#             # Register callback to update progress bar
#             for future in futures:
#                 future.add_done_callback(lambda f: update_result(f.result()))

#     SATCAT.Catalogue = results
#     print("Exporting results...")
#     save_path = os.path.dirname(f'src/data/results/propagated_catalogs/{scenario_name}.pickle')
#     os.makedirs(save_path, exist_ok=True)
#     dump_pickle(f'src/data/results/propagated_catalogs/{scenario_name}.pickle', SATCAT)

#     print(f"Simulation complete. Results saved to: {get_path(f'src/data/results/propagated_catalogs/{scenario_name}.pickle')}")
