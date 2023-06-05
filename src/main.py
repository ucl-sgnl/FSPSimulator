import datetime
import sys
import os
import json
import pickle
from dateutil.parser import parse
import dateutil.relativedelta as rd

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))) 

from src.utils.SpaceCatalogue import SpaceCatalogue
from src.utils.LaunchModel import Prediction2SpaceObjects
from src.utils.coords import utc_to_jd

def run_simple_sim(settings):
    # initialise the catalogue, repull if required
    catalogue = SpaceCatalogue(settings["sim_object_type"], settings["sim_object_catalogue"], settings["repull_catalogues"])
    jd_start = utc_to_jd(settings["sim_start_date"])
    jd_stop = utc_to_jd(settings["sim_end_date"])
    step_size = 60*60*24*365 # 1 year
    policy_name = settings["scenario_name"]

    # pull the catalogue if required
    if settings["repull_catalogues"] and os.path.exists(os.path.join(os.getcwd(), f'src/data/catalogue/All_catalogue_latest.csv')):
        if settings["sim_object_type"] != "all":
            catalogue.CreateCatalogueActive()
        if settings["sim_object_type"] == "debris":
            return # Currently does not exist
        else:
            catalogue.CreateCatalogueAll()

        # create list of space objects from the merged catalogue
        catalogue.Catalogue2SpaceObjects()

        with open(os.path.join(os.getcwd(), f'src/data/catalogue/SATCAT_before_prop.pickle'), 'wb') as f:
            SATCAT_before_prop = catalogue.ReturnCatalogue()
            pickle.dump(SATCAT_before_prop, f)
    else: 
        # set the current catalogue to a previous version
        with open(os.path.join(os.getcwd(), f'src/data/catalogue/SATCAT_before_prop.pickle'), 'rb') as f:
            catalogue.SetCatalogue(pickle.load(f))

    # Only run if results don't already exist
    if settings["environment"] == "development" and os.path.exists(os.path.join(os.getcwd(), f'src/data/catalogue/SATCAT_before_prop_{policy_name}.pickle')):
        with open(os.path.join(os.getcwd(), f'src/data/catalogue/SATCAT_before_prop_{policy_name}.pickle'), 'rb') as f:
            SATCAT_before_prop = pickle.load(f)
    if settings["scenario_name"] == "baseline":
        # don't create a launch model if baseline
        SATCAT_before_prop = catalogue.ReturnCatalogue()
    else:
        # Launch Files
        print("Creating Launch Model...")
        in_file = 'src/data/prediction_csv/FSP_Predictions.csv'
        settings_path = 'src/data/prediction_csv/sim_settings.json'

        # apply policy to launch file and load into UCL Catalogue
        launch_file_object = Prediction2SpaceObjects(in_file, settings_path)
        SATCAT_before_prop = catalogue.ReturnCatalogue()
        SATCAT_before_prop = SATCAT_before_prop + launch_file_object

    # Propagate
    timestep = int(settings["metric_timestep"])*24*60*60
    print("Propagating Satellites...")

    # remove all satellites that have decayed before sim_start_date since we wont be able to propagate them anyway
    decayed_before_start = 0
    for satellite in SATCAT_before_prop:  
        if satellite.decay_date < datetime.datetime.strptime(settings["sim_start_date"], '%Y-%m-%d'):
            SATCAT_before_prop.remove(satellite)
            decayed_before_start += 1
    print("# sats decayed before sim start date: ", decayed_before_start)

    for satellite in SATCAT_before_prop:
        # print("propagating satellite: ", satellite.rso_name)
        satellite.prop_catobjects(jd_start[0], jd_stop[0], timestep) # convert days to seconds

    # Export
    print("Exporting results...")
    with open(os.path.join(os.getcwd(), f'src/data/catalogue/prop_{policy_name}.pickle'), 'wb') as f:
        pickle.dump(SATCAT_before_prop, f)

    print("Output: " + os.path.join(os.getcwd(), f'src/data/catalogue/prop_{policy_name}.pickle'))
    print("Number of Satellites: " + str(len(SATCAT_before_prop)))
    print("Simulation Complete")

def run_sim_timestep(settings, to_csv=False):
    # initialise the catalogue, repull if required
    catalogue = SpaceCatalogue(settings["sim_object_type"], settings["sim_object_catalogue"], settings["repull_catalogues"])
    jd_start = utc_to_jd(settings["sim_start_date"])

    #TODO: explain the format and use of output_frequency in the README or something
    # Convert output_frequency to a relativedelta object

    # if output_frequency is exists use it: 
    step_size = int(settings["output_frequency"][:-1])  # extract the number
    step_unit = settings["output_frequency"][-1]  # extract the unit (d, m, y)

    if step_unit == "d" or step_unit == "D":
        step_delta = rd.relativedelta(days=step_size)
    elif step_unit == "m" or step_unit == "M":
        step_delta = rd.relativedelta(months=step_size)
    elif step_unit == "y" or step_unit == "Y":
        step_delta = rd.relativedelta(years=step_size)

    current_end_date = parse(settings["sim_start_date"])
    final_end_date = parse(settings["sim_end_date"])

    while current_end_date <= final_end_date:
        jd_stop = utc_to_jd(current_end_date.strftime('%Y-%m-%d'))
        step_size = settings["output_frequency"]
        policy_name = settings["scenario_name"]

        # pull the catalogue if required
        if settings["repull_catalogues"] and os.path.exists(os.path.join(os.getcwd(), f'src/data/catalogue/All_catalogue_latest.csv')):
            if settings["sim_object_type"] != "all":
                catalogue.CreateCatalogueActive() 
            if settings["sim_object_type"] == "debris":
                return  # Currently does not exist
            else:
                catalogue.CreateCatalogueAll()  # this will start from space-track

            # create list of space objects from the merged catalogue
            catalogue.Catalogue2SpaceObjects()

            with open(os.path.join(os.getcwd(), f'src/data/catalogue/SATCAT_before_prop.pickle'), 'wb') as f:
                SATCAT_before_prop = catalogue.ReturnCatalogue()
                pickle.dump(SATCAT_before_prop, f)
        else:
            # set the current catalogue to a previous version
            with open(os.path.join(os.getcwd(), f'src/data/catalogue/SATCAT_before_prop.pickle'), 'rb') as f:
                catalogue.SetCatalogue(pickle.load(f))

        # Only run if results don't already exist
        if settings["environment"] == "development" and os.path.exists(os.path.join(os.getcwd(), f'src/data/catalogue/SATCAT_before_prop_{policy_name}.pickle')):
            with open(os.path.join(os.getcwd(), f'src/data/catalogue/SATCAT_before_prop_{policy_name}.pickle'), 'rb') as f:
                SATCAT_before_prop = pickle.load(f)
        if settings["scenario_name"] == "baseline":
            # don't create a launch model if baseline
            SATCAT_before_prop = catalogue.ReturnCatalogue()
        else:
            # Launch Files
            print("Creating Launch Model...")
            in_file = 'src/data/prediction_csv/FSP_Predictions.csv'
            settings_path = 'src/data/prediction_csv/sim_settings.json'

            # apply policy to launch file and load into UCL Catalogue
            launch_file_object = Prediction2SpaceObjects(in_file, settings_path)
            SATCAT_before_prop = catalogue.ReturnCatalogue()
            SATCAT_before_prop = SATCAT_before_prop + launch_file_object

        # Propagate
        timestep = int(settings["metric_timestep"])*24*60*60
        print("Propagating Satellites...")

        # remove all satellites that have decayed
        # TODO: do we really want to remove the decayed satellites?
        # would we not want to keep them for the metadata and just give them no ephemeris?
        for satellite in SATCAT_before_prop:
            try:
                if satellite.decay_date < datetime.datetime.today():
                    SATCAT_before_prop.remove(satellite)
            except TypeError:  # catch when the satellite decay date is None
                continue

        for satellite in SATCAT_before_prop:
            satellite.prop_catobjects(jd_start[0], jd_stop[0], timestep)  # convert days to seconds

        # Export
        if to_csv:
            output_filename = os.path.join(os.getcwd(), f'src/data/catalogue/prop_{policy_name}_{step_unit}{step_size}_{current_end_date.strftime("%Y%m%d")}.pickle')
            with open(output_filename, 'wb') as f:
                pickle.dump(SATCAT_before_prop, f)

            print("Output: " + output_filename)

        print("Number of Satellites: " + str(len(SATCAT_before_prop)))

        # update current_end_date
        current_end_date += step_delta

    print("Simulation Complete")

if __name__ == '__main__':
    with open(os.path.join(os.getcwd(), 'src/data/prediction_csv/sim_settings.json'), 'r') as f:
        json_data = f.read()

    # This just runs the sim to the end
    # Parse the JSON string into a Python object
    policy = json.loads(json_data)
    run_simple_sim(policy)

    #This runs the sim in steps of the specified size
    # Parse the JSON string into a Python object
    # settings = json.loads(json_data)
    # run_sim_timestep(settings, to_csv=False)
