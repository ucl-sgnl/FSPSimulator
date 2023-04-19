# from SpaceObject import Satellite
from UpdateCatalogue import update_catalogue_jsr, update_catalogue_celestrak
from coords import tle_parse
import pandas as pd
import os

class SpaceCatalogue:
    def __init__(self):
        self.PullCatalogue()
        self.Satellites = []
        self.Catalogue = None

    def create_merged_space_catalogue(self):
    # then pull the most recent file form data/external/currentcat.tsv
        jsr_cat = pd.read_csv(os.path.join(os.getcwd(), 'src/data/external/currentcat.tsv'), sep='\t')
        celestrak_cat = pd.read_csv(os.path.join(os.getcwd(), 'src/data/external/celestrak_active.txt'), sep='\n', header=None)

        # convert the celestrak 3LE into one string
        with open(os.path.join(os.getcwd(), 'src/data/external/celestrak_active.txt'), 'r') as f:
            three_line_elements = f.readlines()

        tles = []
        for i in range(0, len(three_line_elements), 3):
            tles.append(three_line_elements[i:i+3])

        tles_parsed= []
        # read in 3LE correctly
        for tle in tles:
            test = ''.join(tle)
            tles_parsed.append(tle_parse(test))


        # first convert the list of tles to a dataframe
        tles_dict = pd.DataFrame(tles_parsed)

        # then merge the two dataframes
        self.Catalogue = pd.merge(jsr_cat, tles_dict, left_on='Satcat', right_on='satellite catalog number')

        self.Catalogue.to_csv(os.path.join(os.getcwd(), 'src/data/external/output_test.tsv'), index=False)
        


    # create a dataframe from the 3LE
    def PullCatalogue(self):
        update_catalogue_jsr()
        update_catalogue_celestrak()
        self.create_merged_space_catalogue()

if __name__ == '__main__':
    catalogue = SpaceCatalogue()

            