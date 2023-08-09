import os
import json
import pickle
import datetime
from tqdm import tqdm
from utils.SpaceCatalogue import SpaceCatalogue, check_json_file
from utils.Conversions import utc_to_jd, initialize_orekit


BUFFER_SIZE = 100


def get_path(*args):
    return os.path.join(os.getcwd(), *args)


def load_pickle(file_path):
    with open(get_path(file_path), 'rb') as f:
        data = []
        while True:
            try:
                data.append(pickle.load(f))
            except EOFError:
                break
        return data


def append_pickle(file_path, data_buffer):
    with open(get_path(file_path), 'ab') as f:
        for item in data_buffer:
            pickle.dump(item, f)


def propagate_space_object(args):
    space_object, jd_start, jd_stop, step_size, output_freq, integrator_type, force_model, long_term_sgp4 = args
    print(f"Propagating {space_object.rso_name}...")
    try:
        space_object.prop_catobject(
            jd_start=jd_start,
            jd_stop=jd_stop,
            step_size=step_size,
            output_freq=output_freq,
            integrator_type=integrator_type,
            force_model=force_model,
            long_term_sgp4=long_term_sgp4
        )
    except Exception as e:
        print(f"Error propagating {space_object.rso_name}: {e}")
        import traceback
        traceback.print_exc()

    return space_object


def run_sim(settings):
    if settings["sgp4_long_term"]:
        initialize_orekit()
    SATCAT = SpaceCatalogue(settings=settings)
    jd_start = float(utc_to_jd(settings["sim_start_date"])[0])
    jd_stop = float(utc_to_jd(settings["sim_end_date"])[0])

    SATCAT.Catalogue = [space_object for space_object in SATCAT.Catalogue if space_object.decay_date >= datetime.datetime.strptime(settings["sim_start_date"], '%Y-%m-%d')]


    print(f"Propagating {len(SATCAT.Catalogue)} space objects...")

    results_buffer = []
    save_path = f'src/data/results/propagated_catalogs/{settings["scenario_name"]}.pickle'
    os.makedirs(os.path.dirname(save_path), exist_ok=True)

    for space_object in tqdm(SATCAT.Catalogue):
        result = propagate_space_object(
            (
                space_object,
                jd_start,
                jd_stop,
                settings["integrator_step_size"],
                settings["output_frequency"],
                settings["integrator_type"],
                settings["force_model"],
                settings["sgp4_long_term"]
            )
        )
        results_buffer.append(result)

        if len(results_buffer) >= BUFFER_SIZE:
            append_pickle(save_path, results_buffer)
            results_buffer.clear()

    # Save remaining results in the buffer
    if results_buffer:
        append_pickle(save_path, results_buffer)

    print(f"Simulation complete. Results saved to: {get_path(save_path)}")


if __name__ == '__main__':

    simulation_files = [sim for sim in os.listdir(get_path('src/data/specify_simulations/')) if sim.endswith('.json')]
    
    for sim in simulation_files:
        print(f"Running simulation: {sim}")
        with open(get_path(f'src/data/specify_simulations/{sim}'), 'r') as file:
            settings = json.load(file)

        check_json_file(settings)
        run_sim(settings)
        print(f"Simulation {sim} complete")


# TODO: PARALLEL DOES NOT WORK
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

#     # Initialize for the main process, but also will be initialized within each child process.
#     initialize_orekit()

#     print("Number of space_object in catalogue specified:", len(SATCAT.Catalogue))
#     print(f"Propagating SpaceObjects and saving state vectors every {settings['output_frequency']} seconds...")

#     decayed_before_start = 0
#     for space_object in SATCAT.Catalogue:
#         if space_object.decay_date < datetime.datetime.strptime(settings["sim_start_date"], '%Y-%m-%d'):
#             SATCAT.Catalogue.remove(space_object)
#             decayed_before_start += 1
#     print("# sats decayed before sim start date:", decayed_before_start)

#     # Slice SATCAT.Catalogue to select every 1000th space object (for testing)
#     SATCAT.Catalogue = SATCAT.Catalogue[::1000]

#     print("Propagating space objects in parallel...")

#     # Create iterable
#     iterable = [(space_object, jd_start, jd_stop, step_size, output_freq, integrator_type, force_model, sgp4_long_term) for space_object in SATCAT.Catalogue]

#     results = []
#     mp.set_start_method('spawn') #This ensures a clean start of each new process
#     with tqdm(total=len(iterable)) as pbar:
#         def callback(result):
#             results.append(result)
#             pbar.update()

#         # Use ProcessPoolExecutor for multiprocessing
#         with ProcessPoolExecutor(max_workers=os.cpu_count()) as executor:
#             # Since the map function blocks until all results are available, we can simply use it here.
#             for result in executor.map(propagate_space_object, iterable):
#                 callback(result)

#     SATCAT.Catalogue = results
#     print("Exporting results...")
#     save_path = os.path.dirname(f'src/data/results/propagated_catalogs/{scenario_name}.pickle')
#     os.makedirs(save_path, exist_ok=True)
#     dump_pickle(f'src/data/results/propagated_catalogs/{scenario_name}.pickle', SATCAT)

#     print(f"Simulation complete. Results saved to: {get_path(f'src/data/results/propagated_catalogs/{scenario_name}.pickle')}")
#     return
