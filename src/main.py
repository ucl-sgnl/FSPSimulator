# if your debugger or python console is not working, try adding the following lines to the top of the file
import datetime
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))) 

import json
import pickle
from src.utils.SpaceObject import SpaceObject
from src.utils.SpaceCatalogue import SpaceCatalogue
from src.utils.LaunchModel import Prediction2SpaceObjects
from src.utils.coords import utc_to_jd

def run_simulation(policy):
    # initialise the catalogue, repull if required
    catalogue = SpaceCatalogue(policy["sim_object_type"], policy["sim_object_catalogue"], policy["repull_catalogues"])
    jd_start = utc_to_jd(policy["sim_start_date"])
    jd_stop = utc_to_jd(policy["sim_end_date"])
    step_size = 60*60*24*365 # 1 year
    policy_name = policy["scenario_name"]

    # pull the catalogue if required
    if policy["repull_catalogues"] and os.path.exists(os.path.join(os.getcwd(), f'src/data/catalogue/All_catalogue_latest.csv')):
        if policy["sim_object_type"] != "all":
            catalogue.CreateCatalogueActive() # this will start with celestrak
        if policy["sim_object_type"] == "debris":
            return # Currently does not exist
        else:
            catalogue.CreateCatalogueAll() # this will start from space-track
             
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
    if policy["environment"] == "development" and os.path.exists(os.path.join(os.getcwd(), f'src/data/catalogue/SATCAT_before_prop_{policy_name}.pickle')):
        with open(os.path.join(os.getcwd(), f'src/data/catalogue/SATCAT_before_prop_{policy_name}.pickle'), 'rb') as f:
            SATCAT_before_prop = pickle.load(f)
    if policy["scenario_name"] == "baseline":
        # don't create a launch model if baseline
        SATCAT_before_prop = catalogue.ReturnCatalogue()
    else:
        # Launch Files
        print("Creating Launch Model...")
        in_file = 'src/data/prediction_csv/FSP_Predictions.csv'
        policy_path = 'src/data/prediction_csv/policy_fsptest.json'

        # apply policy to launch file and load into UCL Catalogue
        launch_file_object = Prediction2SpaceObjects(in_file, policy_path)
        SATCAT_before_prop = catalogue.ReturnCatalogue()
        SATCAT_before_prop = SATCAT_before_prop + launch_file_object

    # Propagate
    timestep = int(policy["metric_timestep"])*24*60*60
    print("Propagating Satellites...")

    # remove all satellites that have decayed
    for satellite in SATCAT_before_prop:
        try:    
            if satellite.decay_date < datetime.datetime.today():
                SATCAT_before_prop.remove(satellite)
        except TypeError: # catch when the satellite decay date is None
            continue 

    for satellite in SATCAT_before_prop:
        print(satellite.__dict__)
        satellite.prop_catobjects(jd_start[0], jd_stop[0], timestep) # convert days to seconds
    
    # Export
    with open(os.path.join(os.getcwd(), f'src/data/catalogue/SATCAT_after_prop_{policy_name}.pickle'), 'wb') as f:
        pickle.dump(SATCAT_before_prop, f)

    print("Output: " + os.path.join(os.getcwd(), f'src/data/catalogue/SATCAT_after_prop_{policy_name}.pickle'))
    print("Number of Satellites: " + str(len(SATCAT_before_prop)))
    print("Simulation Complete")

if __name__ == '__main__':
    with open(os.path.join(os.getcwd(), 'src/data/prediction_csv/policy_fsptest.json'), 'r') as f:
        json_data = f.read()

    # Parse the JSON string into a Python object
    policy = json.loads(json_data)
    run_simulation(policy)