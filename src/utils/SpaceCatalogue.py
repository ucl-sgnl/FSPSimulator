# from SpaceObject import Satellite
from UpdateCatalogue import update_catalogue_jsr, update_catalogue_celestrak, create_merged_catalogue
import pandas as pd
import os

class SpaceCatalogue:
    def __init__(self):
        self.PullCatalogue()
        self.Satellites = []
        self.Catalogue = None

    def PullCatalogue(self):
        # pull in the latest data based from data/external/catalogue.csv
        # update the jsr catalogue if it is older than 1 day
        update_catalogue_jsr()
        update_catalogue_celestrak()
        create_merged_space_catalogue()

    def create_merged_space_catalogue():
    # then pull the most recent file form data/external/currentcat.tsv
        jsr_cat = pd.read_csv(os.path.join(os.getcwd(), 'src/data/external/currentcat.tsv'), sep='\t')
        celestrak_cat = pd.read_csv(os.path.join(os.getcwd(), 'src/data/external/celestrak_active.txt'), sep='\n', header=None)

        # convert the celestrak 3LE into one string
        with open(os.path.join(os.getcwd(), 'src/data/external/celestrak_active.txt'), 'r') as f:
            three_line_elements = [line.strip() for line in f.readlines()]

        tles = []
        for i in range(0, len(three_line_elements), 3):
            tles.append(three_line_elements[i:i+3])

        # read in 3LE correctly

        # merge the two dataframes

if __name__ == '__main__':
    catalogue = SpaceCatalogue()

            