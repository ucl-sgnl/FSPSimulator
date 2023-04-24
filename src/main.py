from src.utils.SpaceObject import SpaceObject
from src.utils.SpaceCatalogue import SpaceCatalogue
from src.utils.LaunchModel import Prediction2SpaceObjects
import json

# def run_simulation(policy):
def run_simulation():
    # id = policy.get('id')
    # print(f'Running simulation for policy {id}')

    # First load in the most recent version of JSR catalogue
    catalogue = SpaceCatalogue()

    # pull down and merge JSR and Celestrak catalogues
    catalogue.CreateMergedSatelliteCatalogue()
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
    
    # read json file
    with open('src/data/prediction_csv/policy_fsptest.json') as f:
        policy = json.load(f)
    
    # Run the simulation
    for satellite in SATCAT_before_prop:
        satellite.sgp4_prop(policy["sim_end_date"])
    # Create graphs

if __name__ == '__main__':
    run_simulation()