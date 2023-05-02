# if your debugger or python console is not working, try adding the following lines to the top of the file
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))) 
import json
import pickle
from src.utils.SpaceObject import SpaceObject
from src.utils.SpaceCatalogue import SpaceCatalogue
from src.utils.LaunchModel import Prediction2SpaceObjects
from src.utils.coords import utc_to_jd
import timeit

def run_simulation(policy):
    # First load in the most recent version of JSR catalogue, this will occur in initialisation of the SpaceCatalogue Class
    catalogue = SpaceCatalogue(policy["sim_object_type"], policy["sim_object_catalogue"])
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


    # check if the simulation has already been run
    if policy["environment"] == "development" and os.path.exists(os.path.join(os.getcwd(), f'src/data/catalogue/SATCAT_before_prop_{policy_name}.pickle')):
        # as an updated satellite object already exists, load it in
        with open(os.path.join(os.getcwd(), f'src/data/catalogue/SATCAT_before_prop_{policy_name}.pickle'), 'rb') as f:
            SATCAT_before_prop = pickle.load(f)
    else:
        # Load in up to date future launch files
        print("Creating Launch Model...")
        in_file = 'src/data/prediction_csv/04_04_23_fsp.csv'
        policy = 'src/data/prediction_csv/policy_fsptest.json'
        # applies policy to the launch model and then creates space objects
        launch_file_object = Prediction2SpaceObjects(in_file, policy)

        # Load it into the UCL Space Catalogue
        SATCAT_before_prop = catalogue.ReturnCatalogue()
        SATCAT_before_prop = SATCAT_before_prop + launch_file_object

    # Run the simulation
    print("Propogating Satellites...")
    for satellite in SATCAT_before_prop:
        satellite.sgp4_prop_catobjects(jd_start[0], jd_stop[0], step_size)

    # save the final Satellite Catalogue in the data folder
    with open(os.path.join(os.getcwd(), f'src/data/catalogue/SATCAT_after_prop_{policy_name}.pickle'), 'wb') as f:
        pickle.dump(SATCAT_before_prop, f)

if __name__ == '__main__':
    with open(os.path.join(os.getcwd(), 'src/data/prediction_csv/policy_fsptest.json'), 'r') as f:
        json_data = f.read()

    # Parse the JSON string into a Python object
    policy = json.loads(json_data)
    run_simulation(policy)