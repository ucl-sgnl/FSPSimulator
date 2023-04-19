# from SpaceObject import Satellite
from UpdateCatalogue import update_catalogue_jsr, update_catalogue_celestrak
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

        # then pull the most recent file form data/external/currentcat.tsv
        jsr_cat = pd.read_csv(os.path.join(os.getcwd(), 'src/data/external/currentcat.tsv'), sep='\t')
        print(jsr_cat.head())

       

if __name__ == '__main__':
    catalogue = SpaceCatalogue()

            