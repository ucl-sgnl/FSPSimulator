from utils.SpaceObject import Satellite
from utils.SpaceCatalogue import SpaceCatalogue
import pandas as pd

# def run_simulation(policy):
def run_simulation():
    # id = policy.get('id')
    # print(f'Running simulation for policy {id}')

    # First load in the most recent version of JSR catalogue
    space_catalogue = SpaceCatalogue()

    # Load it into the UCL Space Catalogue
    currentcat = pd.read_csv("src\data\external\currentcat.tsv", sep='\t')
    space_catalogue.AddSatellitesTSV(currentcat)


if __name__ == '__main__':
    run_simulation()