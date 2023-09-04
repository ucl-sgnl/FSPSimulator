import os
import json
import pickle
from tqdm import tqdm
from fspsim.utils.SpaceCatalogue import SpaceCatalogue, check_json_file
from fspsim.utils.Conversions import utc_to_jd

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
    #print(f"Propagating {space_object.rso_name}...")
    try:
        space_object.prop_catobject(jd_start=jd_start, jd_stop=jd_stop, step_size=step_size, output_freq=output_freq, integrator_type=integrator_type, force_model=force_model, long_term_sgp4=long_term_sgp4)
    except Exception as e:
        print(f"An error occurred while propagating {space_object.rso_name}: {e}")

    return

def run_sim(settings):
    SATCAT = SpaceCatalogue(settings=settings)
    total_objects = len(SATCAT.Catalogue)
    
    jd_start = float(utc_to_jd(settings["sim_start_date"]))
    jd_stop = float(utc_to_jd(settings["sim_end_date"]))
    step_size = int(settings["integrator_step_size"])
    output_freq = int(settings["output_frequency"])
    scenario_name = str(settings["scenario_name"])
    integrator_type = str(settings["integrator_type"])
    sgp4_long_term = bool(settings["sgp4_long_term"])
    force_model = settings["force_model"]

    batch = 1
    batch_size = 100000000

    SATCAT.Catalogue = SATCAT.Catalogue[::1000] # Slice for testing   

    # remove satellites that have an altitude higher than 2000, or are HEO
    for sat in SATCAT.Catalogue:
        if (sat.apogee + sat.perigee) / 2 > 2000 or sat.eccentricity > 0.2:
            SATCAT.Catalogue.remove(sat)

    # Create a progress bar
    pbar = tqdm(total=len(SATCAT.Catalogue), desc="Propagating")

    while SATCAT.Catalogue:
        current_batch = []

        # Propagate space objects for the current batch
        for _ in range(batch_size):
            if not SATCAT.Catalogue:
                break

            space_object = SATCAT.Catalogue.pop(0)
            propagate_space_object((space_object, jd_start, jd_stop, step_size, output_freq, integrator_type, force_model, sgp4_long_term))
            current_batch.append(space_object)
            pbar.update(1)  # Update the progress bar

        # Save the current batch
        print(f"Saving batch {batch}...")
        save_path = os.path.join(f'src/fspsim/data/results/propagated_catalogs/', f"{scenario_name}_batch_{batch}.pickle")
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        dump_pickle(save_path, current_batch)
        batch += 1

    pbar.close()  # Close the progress bar
    print(f"Simulation complete. Results saved to batches in: {get_path(f'src/fspsim/data/results/propagated_catalogs/')}")

if __name__ == '__main__':
    #list all the json files in src/fspsim/data/specify_simulations
    sims = os.listdir(get_path('src/fspsim/data/specify_simulations/'))
    for sim in sims:
        if sim.endswith('.json'):
            print(f"Running simulation: {sim}")
            settings = json.load(open(get_path(f'src/fspsim/data/specify_simulations/{sim}'), 'r'))
            check_json_file(settings)#check if the json file is filled out correctly
            run_sim(settings)
            print(f"Simulation {sim} complete")