# if your debugger or python console is not working, try adding the following lines to the top of the file
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))) 
import json
import pickle
from src.utils.SpaceObject import SpaceObject
from src.utils.SpaceCatalogue import SpaceCatalogue
from src.utils.LaunchModel import Prediction2SpaceObjects

def run_simulation(policy):
    # First load in the most recent version of JSR catalogue, this will occur in initialisation of the SpaceCatalogue Class
    catalogue = SpaceCatalogue(policy["sim_object_type"], policy["sim_object_catalogue"])

    if policy["environment"] == "development" and os.path.exists(os.path.join(os.getcwd(), 'src/data/catalogue/SATCAT_before_prop.pickle')):
        # as an updated satellite object already exists, load it in
        with open(os.path.join(os.getcwd(), 'src/data/catalogue/SATCAT_before_prop.pickle'), 'rb') as f:
            # Load the data from the pickle file
            SATCAT_before_prop = pickle.load(f)
        # Run the simulation
        for sat in SATCAT_before_prop:
            sat.sgp4_prop(policy["sim_end_date"])
    else:
        # create a merged catalogue based on the user definition
        if policy["sim_object_type"] != "all":
            catalogue.CreateCatalogueActive() # this will start with celestrak
        else:
            catalogue.CreateCatalogueAll() # this will start from space-track

        # create list of space objects from the merged catalogue
        catalogue.Catalogue2SpaceObjects()

        # Load in up to date future launch files
        in_file = 'src/data/prediction_csv/04_04_23_fsp.csv'
        policy = 'src/data/prediction_csv/policy_fsptest.json'
        # applies policy to the launch model and then creates space objects
        launch_file_object = Prediction2SpaceObjects(in_file, policy)

        # Load it into the UCL Space Catalogue
        SATCAT_before_prop = catalogue.ReturnCatalogue()
        SATCAT_before_prop = SATCAT_before_prop + launch_file_object

        # save the catalogue to a pickle file in the src/data folder
        with open(os.path.join(os.getcwd(), 'src/data/catalogue/SATCAT_before_prop.pickle'), 'wb') as f:
            pickle.dump(SATCAT_before_prop, f)
        
        # Run the simulation
        for satellite in SATCAT_before_prop:
            satellite.sgp4_prop(policy["sim_end_date"])
        
    # Create graphs

if __name__ == '__main__':
    with open(os.path.join(os.getcwd(), 'src/data/prediction_csv/policy_fsptest.json'), 'r') as f:
        json_data = f.read()

    # Parse the JSON string into a Python object
    policy = json.loads(json_data)
    run_simulation(policy)