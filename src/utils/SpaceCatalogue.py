from utils.SpaceObject import Satellite

class SpaceCatalogue:
    def __init__(self):
        self.PullCatalogue()
        self.Satellites = []

    def PullCatalogue(self):
        # pull in the latest data based from data/external/catalogue.csv
        pass
    
    # Adds JSR satellites to the catalogue
    def AddSatellitesTSV(self, catalogue):
        catalogue = catalogue.drop(index=0) # drop the header file
        for index, row in catalogue.iterrows():
            print(row["Period"])
            