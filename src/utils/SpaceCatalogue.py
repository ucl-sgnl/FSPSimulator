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
        self.CurrentCatalogue = None

    # Pulls the most recent version of the merged catalogue
    def ReturnCatalogue(self):
        return self.Catalogue
    
    # Merges the JSR and Celestrak catalogues and creates a list of SpaceObjects
    def CreateMergedSatelliteCatalogue(self):
        # then pull the most recent file form data/external/currentcat.tsv
        jsr_cat = pd.read_csv(os.path.join(os.getcwd(), 'src/data/external/currentcat.tsv'), sep='\t')
        jsr_cat_extra_info = pd.read_csv(os.path.join(os.getcwd(), 'src/data/external/satcat.tsv'), sep='\t')
        with open(os.path.join(os.getcwd(), 'src/data/external/celestrak_active.txt'), 'r') as f:
            three_line_elements = f.readlines()

        # merge the each three line element into one string
        tles = []
        for i in range(0, len(three_line_elements), 3):
            tles.append(three_line_elements[i:i+3])

        # read in 3LE correctly
        tles_parsed= []
        for tle in tles:
            test = ''.join(tle)
            tles_parsed.append(tle_parse(test))

        # first convert the list of tles to a dataframe
        tles_dict = pd.DataFrame(tles_parsed)

        # then merge the two dataframes, firstly removing the first row of the jsr_cat
        jsr_cat = jsr_cat.drop(0)
        self.CurrentCatalogueDF = tles_dict.merge(jsr_cat, right_on='Satcat', left_on='satellite catalog number')
        self.CurrentCatalogueDF = self.CurrentCatalogueDF.merge(jsr_cat_extra_info, right_on='#JCAT', left_on='#JCAT')

        # clean the catalogue in a format that can be used by the SpaceObject class
        # Get list of non-duplicate column names
        cols_to_keep = self.CurrentCatalogueDF.filter(regex='^(?!.*_y)').columns
        self.CurrentCatalogueDF = self.CurrentCatalogueDF.loc[:, cols_to_keep]
        new_cols = {col: col.replace('_x', '') for col in self.CurrentCatalogueDF.columns}
        self.CurrentCatalogueDF = self.CurrentCatalogueDF.rename(columns=new_cols)
        
        # remove all the satellites that have already been launched, or dont have a perigee or apogee
        self.CurrentCatalogueDF = self.CurrentCatalogueDF[self.CurrentCatalogueDF['DDate'] == '-']
        self.CurrentCatalogueDF = self.CurrentCatalogueDF[self.CurrentCatalogueDF['Apogee'] != '-']
        self.CurrentCatalogueDF = self.CurrentCatalogueDF[self.CurrentCatalogueDF['Perigee'] != '-']     

        # Export for audit purposes
        self.CurrentCatalogueDF.to_csv(os.path.join(os.getcwd(), 'src/data/external/active_jsr_celestrak_latest.csv'), index=False)
    
    # Updates the most recent versions of celestrak and JSR's data
    def PullCatalogue(self):
        # update_catalogue_jsr()
        # update_catalogue_celestrak()
        pass # use this for testing as it takes ages

    # Creates a list of SpaceObjects from the merged catalogue
    def Catalogue2SpaceObjects(self):
        # headings are as follows:
            #     Index(['line number', 'satellite catalog number', 'classification',
            #    'International Designator(launch year)',
            #    'International Designator (launch num)',
            #    'International Designator (piece of launch)', 'epoch year', 'epoch day',
            #    'first time derivative of mean motion(ballisitc coefficient)',
            #    'second time derivative of mean motion(delta-dot)', 'bstar drag term',
            #    'ephemeris type', 'element number', 'checksum', 'inclination',
            #    'right ascension of the ascending node', 'eccentricity',
            #    'argument of perigee', 'mean anomaly', 'mean motion',
            #    'revolution number at epoch', '#JCAT', 'DeepCat', 'Satcat', 'Piece',
            #    'Active', 'Type', 'Name', 'LDate', 'Parent', 'Owner', 'State', 'SDate',
            #    'ExpandedStatus', 'DDate', 'ODate', 'Period', 'Perigee', 'PF', 'Apogee',
            #    'AF', 'Inc', 'IF', 'OpOrbit', 'PLName', 'Primary', 'Status', 'Dest',
            #    'Manufacturer', 'Bus', 'Motor', 'Mass', 'MassFlag', 'DryMass',
            #    'DryFlag', 'TotMass', 'TotFlag', 'Length', 'LFlag', 'Diameter', 'DFlag',
            #    'Span', 'SpanFlag', 'Shape', 'OQUAL', 'AltNames'],
            #   dtype='object')   
        for index, row in self.CurrentCatalogueDF.iterrows():
                self.Catalogue.append(SpaceObject(  object_type=row['Type'], 
                                                    payload_operational_status=row['Active'], ## this will need to be updated to celestrak's active status
                                                    application="Unknown", 
                                                    operator=row['Owner'], 
                                                    mass=row['Mass'], 
                                                    # type = row['Type'],
                                                    # maneuverable=maneuverable, 
                                                    # spin_stabilized=spin_stabilized, 
                                                    # radar_cross_section=radar_cs, 
                                                    # characteristic_area=area, 
                                                    # characteristic_length=length, 
                                                    # propulsion_type=propulsion, 
                                                    # sma=row[], 
                                                    eccentricity=row['eccentricity'], 
                                                    inc=row['inclination'], 
                                                    argp=row['argument of perigee'], 
                                                    raan=row['right ascension of the ascending node'], 
                                                    # launch_site=launch_site,
                                                    source = row['Owner'],
                                                    launch_date=row['LDate'], 
                                                    decay_date=row['DDate'], 
                                                    cospar_id=row['Piece'],
                                                    rso_name=row['Name'],
                                                    perigee_altitude=row['Perigee'],
                                                    apogee_altitude=row['Apogee']
                                                ))
        return self.Catalogue

if __name__ == '__main__':
    catalogue = SpaceCatalogue()
    catalogue.CreateMergedSatelliteCatalogue()
    catalogue.Catalogue2SpaceObjects()

