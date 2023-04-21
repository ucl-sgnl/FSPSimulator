from src.utils.SpaceObject import SpaceObject
from src.utils.SpaceCatalogue import SpaceCatalogue
from src.utils.LaunchModel import Prediction2SpaceObjects

# def run_simulation(policy):
def run_simulation():
    # id = policy.get('id')
    # print(f'Running simulation for policy {id}')

    # First load in the most recent version of JSR catalogue
    catalogue = SpaceCatalogue()
    catalogue.create_merged_space_catalogue()
    catalogue.Catalogue2SpaceObjects()

    # Load in up to date launch files
    in_file = 'src/data/prediction_csv/04_04_23_fsp.csv'
    policy = 'src/data/prediction_csv/policy_fsptest.json'
    launch_file_object = Prediction2SpaceObjects(in_file, policy)

    # Load it into the UCL Space Catalogue
    currentcat_objects = catalogue.Catalogue2SpaceObjects("src/data/prediction_csv/policy_fsptest.json")

    # Run the simulation

if __name__ == '__main__':
    run_simulation()