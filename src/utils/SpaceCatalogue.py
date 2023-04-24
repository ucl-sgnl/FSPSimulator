# from SpaceObject import Satellite
from src.utils.UpdateCatalogue import update_catalogue_jsr, update_catalogue_celestrak
from src.utils.SpaceObject import SpaceObject
from src.utils.Coords import tle_parse
import pandas as pd
import os

class SpaceCatalogue:
    def __init__(self):
        self.PullCatalogue()
        self.Satellites = []
        self.Catalogue = []
        self.CurrentCatalogue = []

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
        tles_dict.set_index('satellite catalog number', inplace=True)

        # then merge the two dataframes
        self.CurrentCatalogue = tles_dict.merge(jsr_cat, right_on='Satcat', left_on='satellite catalog number')

        print(self.CurrentCatalogue)

        self.CurrentCatalogue.to_csv(os.path.join(os.getcwd(), 'src/data/external/active_jsr_celestrak_latest.csv'), index=False)
        
    # create a dataframe from the 3LE
    def PullCatalogue(self):
        # update_catalogue_jsr()
        # update_catalogue_celestrak()
        pass # use this for testing as it takes ages

    def Catalogue2SpaceObjects(self):
        # for each satellite in the current catalogue, create a Space Object Instance
        for index, row in self.CurrentCatalogue.iterrows():
            # line number,classification,International Designator(launch year),International Designator (launch num),International Designator (piece of launch),
            # epoch year,epoch day,first time derivative of mean motion(ballisitc coefficient),second time derivative of mean motion(delta-dot),bstar drag term,
            # ephemeris type,element number,checksum,inclination,right ascension of the ascending node,eccentricity,argument of perigee,mean anomaly,mean motion,
            # revolution number at epoch,

            #JCAT,DeepCat,Satcat,Piece,Active,Type,Name,LDate,Parent,Owner,State,SDate,ExpandedStatus,DDate,ODate,Period,Perigee,PF,
            # Apogee,AF,Inc,IF,OpOrbit     
                self.Catalogue.append(SpaceObject(object_type=row['Type'], 
                                                    payload_operational_status=row['Active'], ## this will need to be updated to celestrak's active status
                                                    application="Unknown", 
                                                    operator=row['Owner'], 
                                                    # mass=mass, 
                                                    # maneuverable=maneuverable, 
                                                    # spin_stabilized=spin_stabilized, 
                                                    # radar_cross_section=radar_cs, 
                                                    # characteristic_area=area, 
                                                    # characteristic_length=length, 
                                                    # propulsion_type=propulsion, 
                                                    # sma=row[], 
                                                    eccentricity=row['eccentricity'], inc=row['inclination'], argp=row['argument of perigee'], 
                                                    raan=row['right ascension of the ascending node'], 
                                                    # tran=np.deg2rad(TRAN_n), 
                                                    # launch_site=launch_site, 
                                                    launch_date=row['LDate'], 
                                                    decay_date=row['DDate'], 
                                                    cospar_id=row['Piece'],
                                                    rso_name=row['Name'],
                                                    perigee_altitude=row['Perigee'],
                                                    apogee_altitude=row['Apogee'])
                                                    )
        return self.Catalogue
    
    def CleanCatalogue():
        """
            As there are discrepencies in each of the catalogues, this function will clean it. Either by removing satellties that have decayed
            Or by adding in extra information if non-existant, specifically the mass of the satellite
        """

    def ReturnCatalogue(self):
        return self.Catalogue

if __name__ == '__main__':
    catalogue = SpaceCatalogue() # for testing purposes only
    catalogue.create_merged_space_catalogue()
    catalogue.Catalogue2SpaceObjects()

