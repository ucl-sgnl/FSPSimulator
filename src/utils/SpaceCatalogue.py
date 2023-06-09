import os
import json
import requests
import wget
import pickle
import pandas as pd
import datetime
from utils.SpaceObject import SpaceObject
from utils.Conversions import tle_parse
from utils.LaunchModel import Prediction2SpaceObjects

import json

def check_json_file(json):
    # Define the valid keys and their types
    expected_keys = {
        "scenario_name": str,
        "monthly_ton_capacity": str,
        "launch_start_date": str,
        "remove_operators": str,
        "sim_start_date": str,
        "sim_end_date": str,
        "output_frequency": int,
        "integrator_step_size": int,
        "integrator_type": str,
        "sim_object_type": str,
        "sim_object_catalogue": str,
        "environment": str,
        "repull_catalogues": bool,
        "satellite_predictions_csv": str
    }

    # Define the valid values for keys with a limited set of valid options
    valid_values = {
        "integrator_type": ["RK45", "RK23", "DOP853", "Radau", "BDF", "LSODA"],
        "sim_object_type": ["all", "debris", "active"],
        "sim_object_catalogue": ["jsr", "spacetrack", "both"],
    }

    for key, expected_type in expected_keys.items():
        if key not in json:
            print(f"Error: Missing key '{key}' in JSON file.")
            return False
        if not isinstance(json[key], expected_type):
            print(f"Error: Incorrect type for key '{key}'. Expected {expected_type.__name__}, but got {type(json[key]).__name__}.")
            return False
        if key in valid_values and json[key] not in valid_values[key]:
            print(f"Error: Invalid value for key '{key}'. Expected one of {valid_values[key]}, but got {json[key]}.")
            return False

    print("JSON file is valid.")
    return True

def get_path(*args):
    return os.path.join(os.getcwd(), *args)

def load_pickle(file_path):
    with open(get_path(file_path), 'rb') as f:
        loaded_data = pickle.load(f)
        return loaded_data
    
def dump_pickle(file_path, data):
    with open(get_path(file_path), 'wb') as f:
        pickle.dump(data, f)

class SpaceCatalogue:
    def __init__(self, settings):
        satellite_predictions_csv = str("src/data/prediction_csv/") + settings["satellite_predictions_csv"]
        sim_object_type = settings["sim_object_type"] # this can be "active", "all", or "debris"
        sim_object_catalogue = settings["sim_object_catalogue"] # this can be "jsr", "spacetrack", or "both"
        repull_catalogues = settings["repull_catalogues"]
        self.Satellites = []
        self.Catalogue = []
        self.CurrentCatalogue = None
        self.sim_object_type = sim_object_type # this can be "active", "all", or "debris"
        self.sim_object_catalogue = sim_object_catalogue # this can be "jsr", "spacetrack", or "both"
        # raise exception if invalid sim_object_type or sim_object_catalogue is specified
        if self.sim_object_type not in ["active", "all", "debris"]:
            raise Exception("Invalid sim_object_type specified, must be 'active', 'all', or 'debris'")
        if self.sim_object_catalogue not in ["jsr", "spacetrack", "both"]:
            raise Exception("Invalid sim_object_catalogue specified, must be 'jsr', 'spacetrack', or 'both'")
        self.repull_catalogues = bool(repull_catalogues)

        # If we are repulling the catalogues call the appropriate function depending on the sim_object_catalogue
        if self.repull_catalogues == True:
            if self.sim_object_catalogue == "jsr":
                print("re-downloading JSR Cat")
                self.PullCatalogueJSR()
            elif self.sim_object_catalogue == "spacetrack":
                print("re-downloading SpaceTrack Cat")
                self.PullCatalogueSpaceTrack()
            elif self.sim_object_catalogue == "both":
                print("re-downloading JSR and SpaceTrack Cat")
                self.PullCatalogueJSR()
                self.PullCatalogueSpaceTrack()

            print("saving new catalogue data to pickle")
            dump_pickle('src/data/catalogue/SATCAT_before_prop.pickle', self)

        # If we are not repulling the catalogues just use the existing local catalogues
        else:
            print("using existing local catalogues")
            loaded_data = self.load_from_file('src/data/catalogue/SATCAT_before_prop.pickle')
            self.__dict__ = loaded_data.__dict__

        # Now use the catalogues to create a list of SpaceObjects which will be appended to the Catalogue list
        if self.sim_object_type == "active":
            self.CreateCatalogueActive()
        elif self.sim_object_type == "all":
            self.CreateCatalogueAll()
        elif self.sim_object_type == "debris":
            raise Exception("Debris not yet implemented")
        
        self.Catalogue2SpaceObjects()
        
        # Now add the predictions to the Catalogue attribute of the SpaceCatalogue instance by making a list of SpaceObjects using Prediction2SpaceObjects
        predicted_space_objects = Prediction2SpaceObjects(satellite_predictions_csv = satellite_predictions_csv, simsettings=settings)
        self.Catalogue.extend(predicted_space_objects)
        return None
    
    @classmethod
    def load_from_file(cls, file_path):
        with open(get_path(file_path), 'rb') as f:
            loaded_data = pickle.load(f)
            if isinstance(loaded_data, cls):
                return loaded_data
            else:
                raise TypeError(f"When loading catalogue from file, expected a {cls.__name__} instance, but got {type(loaded_data)}")
            
    def CreateCatalogueActive(self):
        """
        This function will merge the JSR and SpaceTrack catalogues
        
        ### Exports
        - List of merged space objects to 'src/data/external/active_jsr_spacetrack.csv'
        """
        # This merges the JSR and Spacetrack active catalogues
        jsr_cat = pd.read_csv(os.path.join(os.getcwd(), 'src/data/external/currentcat.tsv'), sep='\t')
        jsr_cat_extra_info = pd.read_csv(os.path.join(os.getcwd(), 'src/data/external/satcat.tsv'), sep='\t')
        with open(os.path.join(os.getcwd(), 'src/data/external/spacetrack_active.txt'), 'r') as f:
            three_line_elements = f.readlines()

        # merge the each three line element into one string
        tles = []
        for i in range(0, len(three_line_elements), 3):
            tles.append(three_line_elements[i:i+3])

        # read in 3LE correctly, this will also have the entire TLE string as an extra parameter
        tles_parsed= []
        for tle in tles:
            temp = ''.join(tle)
            tles_parsed.append(tle_parse(temp))

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
        # Open the file in 'w' mode (write mode, overwrite if exists)
        with open(os.path.join(os.getcwd(), 'src/data/catalogue/active_jsr_spacetrack_latest.csv'), 'w') as f:
            # Write the header and data to the file
            self.CurrentCatalogueDF.to_csv(f, index=False)

    def CreateCatalogueAll(self):
        """
        Will use Space-track as a base and merge JSR for the active satellites that we have information on. 

        ### Exports
        - Space Catalogue of all tracked objects by Space-track
            - 'src/data/catalogue/All_catalogue_latest.txt'
        """ 
        # Space Track's catalogue is a json
        spacetrack = pd.read_json(os.path.join(os.getcwd(), 'src/data/external/spacetrack_all.json'))
        print("Number of satellites in spacetrack catalogue: ", len(spacetrack))
        jsr_cat_extra_info = pd.read_csv(os.path.join(os.getcwd(), 'src/data/external/satcat.tsv'), sep='\t')
        print("Number of satellites in jsr catalogue: ", len(jsr_cat_extra_info))
        # merge the two dataframes, keeping all items in spacetrack dataframe
        self.CurrentCatalogueDF = spacetrack.merge(jsr_cat_extra_info, right_on='Piece', left_on='OBJECT_ID', how='left')
        print("Number of satellites in catalogues after merge: ", len(self.CurrentCatalogueDF))
        # Count the number of rows that didn't match
        unmatched_rows = self.CurrentCatalogueDF['Piece'].isnull().sum()
        print("Number of Objects in spacetrack catalogue that didn't match jsr catalogue: ", unmatched_rows)
        # Drop the rows that didn't match
        self.CurrentCatalogueDF = self.CurrentCatalogueDF.dropna(subset=['Piece'])
        print("Number of satellites in catalogue after dropping unmatching rows: ", len(self.CurrentCatalogueDF))

        # convert the BStar to scientific notation
        # self.CurrentCatalogueDF['BSTAR'] = self.CurrentCatalogueDF['BSTAR'].map('{:.7e}'.format)
        
        # drop satellites that have a perigee or apogee > 100,000 km
        self.CurrentCatalogueDF = self.CurrentCatalogueDF[self.CurrentCatalogueDF['PERIAPSIS'] < 100000]
        self.CurrentCatalogueDF = self.CurrentCatalogueDF[self.CurrentCatalogueDF['APOAPSIS'] < 100000]
        # make any pergigees or apogees that are negative or 0 set to 1.69420 -> this will keep them in the catalog but also allow us to find them out later
        self.CurrentCatalogueDF['PERIAPSIS'] = self.CurrentCatalogueDF['PERIAPSIS'].mask(self.CurrentCatalogueDF['PERIAPSIS'] <= 0, 1.69420)
        self.CurrentCatalogueDF['APOAPSIS'] = self.CurrentCatalogueDF['APOAPSIS'].mask(self.CurrentCatalogueDF['APOAPSIS'] <= 0, 1.69420)

        # find all the rows that have a decay date (if not they contain None)
        rows_with_decay_date = self.CurrentCatalogueDF[self.CurrentCatalogueDF['DECAY_DATE'].notnull()]
        # now check that these decay dates are before the launch date
        # first we need to convert the dates to datetime objects
        rows_with_decay_date['LAUNCH_DATE'] = pd.to_datetime(rows_with_decay_date['LAUNCH_DATE'])
        rows_with_decay_date['DECAY_DATE'] = pd.to_datetime(rows_with_decay_date['DECAY_DATE'])
        # now we can check if the decay date is before the launch date
        rows_with_bad_decay_date = rows_with_decay_date[rows_with_decay_date['DECAY_DATE'] < rows_with_decay_date['LAUNCH_DATE']]
        # now we can drop these rows from the dataframe
        self.CurrentCatalogueDF = self.CurrentCatalogueDF.drop(rows_with_bad_decay_date.index)

        # Where there is any conflicting data, we will use Space-track's data. 
        # Drop unnecessary columns that exist in Space-track data
        unused_cols = ['Name', # OBJECT_NAME is the same
                           'LDate', # LAUNCH_DATE is the same
                           'DDate', # DECAY_DATE is the same
                           'Perigee', # PERIAPSIS is the same
                            'Apogee', # APOAPSIS is the same
                            'Inc', # INCLINATION is the same
                            'Type', # OBJECT_TYPE is the same
                            'TotMass', # we are just going to use MASS for simplicity
                            'MassFlag', # we are just going to use MASS for simplicity
                            'DryMass', # we are just going to use MASS for simplicity
                            'CCSDS_OMM_VERS', # not needed
                            'COMMENT', # not needed
                            'OBJECT_ID', # not needed we have OBJECT_NAME and NORAD_CAT_ID
                            'CENTER_NAME', # not needed
                            'REF_FRAME', # not needed
                            'TIME_SYSTEM', # not needed
                            'EPHEMERIS_TYPE', # not sure what this is
                            'RCS_SIZE', # not needed
                            'FILE', # not needed
                            'GP_ID', # not needed
                            '#JCAT', # not needed this is just the JCAT internal ID
                            'Satcat', # not needed as we merge on OBJECT_ID. Mostly empty anyway.
                            'Piece', # not needed
                            'Type', # not needed as we have OBJECT_TYPE
                            'Name', # not needed as we have OBJECT_NAME
                            'PLName', # not needed
                            'Parent', # not needed
                            'SDate', # not needed
                            'Primary', # not needed
                            'Dest', # not needed
                            'Manufacturer', # not needed
                            'Bus', # not needed
                            'Motor', # not needed
                            'DryFlag', # not needed
                            'TotFlag', # not needed
                            'TotMass', # not needed
                            'LFlag', # not needed
                            'DFlag', # not needed
                            'Span', # not needed
                            'SpanFlag', # not needed
                            'Shape', # not needed
                            'IF', # not needed
                            'OQUAL', # not needed
                            'AltNames'] # not needed
        self.CurrentCatalogueDF = self.CurrentCatalogueDF.drop(columns=unused_cols)

        # now we will drop rows that do not have certain data as they are required for the simulation
        # drop a row if it doenst have : Launch date, decay date, perigee, apogee, eccentricity, inclination, right ascension of the ascending node, argument of perigee
        self.CurrentCatalogueDF = self.CurrentCatalogueDF.dropna(subset=['LAUNCH_DATE', 'PERIAPSIS', 'APOAPSIS', 'ECCENTRICITY', 'INCLINATION', 'RA_OF_ASC_NODE', 'ARG_OF_PERICENTER'])
        print("Number of satellites in catalogue after dropping rows that were missing data required for simulation: ", len(self.CurrentCatalogueDF))
        #TODO: ~10,000 rows are dropped here. Need to investigate why this is happening.  
        # Export the cleaned catalogue for sanity checking
        self.CurrentCatalogueDF.to_csv(os.path.join(os.getcwd(), 'src/data/catalogue/All_catalogue_latest.csv'), index=False)
        
    def Catalogue2SpaceObjects(self):
        """
        This function will convert the current catalogue into a list of SpaceObjects

        ### Returns
        - List of SpaceObjects
            - View the SpaceObject class for more information
        """
        if self.sim_object_type == "active":
            for index, row in self.CurrentCatalogueDF.iloc[1:].iterrows():
                    self.Catalogue.append(SpaceObject(  object_type=row['Type'], 
                                                        payload_operational_status=row['Active'],
                                                        application="Unknown", 
                                                        operator=row['Owner'], 
                                                        mass=row['Mass'], 
                                                        eccentricity=row['eccentricity'], 
                                                        inc=row['inclination'], 
                                                        argp=row['argument of perigee'], 
                                                        raan=row['right ascension of the ascending node'], 
                                                        source = row['Owner'],
                                                        launch_date=row['LDate'], 
                                                        decay_date=row['DDate'], 
                                                        cospar_id=row['Piece'],
                                                        rso_name=row['Name'],
                                                        perigee=row['Perigee'],
                                                        apogee=row['Apogee'],
                                                        tle=row['TLE']
                                                    ))
        elif self.sim_object_type == "all":
            # the dataframe we are looking at is the merged JSR and SpaceTrack catalogue
            # it contains the columns from both catalogues
            print("building space objects from JSR/Space-track merged catalogue")
            for _, row in self.CurrentCatalogueDF.iloc[1:].iterrows():
                tle = row["TLE_LINE1"] + "\n" + row["TLE_LINE2"]
                self.Catalogue.append(SpaceObject(object_type=row['OBJECT_TYPE'], 
                                                        mass=row['Mass'], 
                                                        launch_site=row["SITE"],
                                                        sma=row['SEMIMAJOR_AXIS'], 
                                                        eccentricity=row['ECCENTRICITY'], 
                                                        inc=row['INCLINATION'], 
                                                        argp=row['ARG_OF_PERICENTER'], 
                                                        raan=row['RA_OF_ASC_NODE'], 
                                                        source = row['COUNTRY_CODE'],
                                                        launch_date=row['LAUNCH_DATE'], 
                                                        decay_date=row['DECAY_DATE'], 
                                                        rso_name=row['OBJECT_NAME'],
                                                        perigee=row['PERIAPSIS'],
                                                        apogee=row['APOAPSIS'],
                                                        tle=tle,
                                                        epoch=row['EPOCH']
                                                    ))
        else:
            raise Exception("Invalid sim_object_type specified")
        return self.Catalogue

    def DownloadJSRCatalogueIfNewer(self, local_path, url):
        """Download a file from a URL if it is newer than the local file."""
        if os.path.exists(local_path):
            local_last_modified = datetime.datetime.fromtimestamp(os.path.getmtime(local_path))

            response = requests.head(url)
            remote_last_modified = datetime.datetime.strptime(response.headers['Last-Modified'], '%a, %d %b %Y %H:%M:%S %Z')

            if remote_last_modified <= local_last_modified:
                return

        wget.download(url, local_path)

    def PullCatalogueJSR(self):
        """Update the JSR catalogue files."""
        cwd = os.getcwd()
        external_dir = os.path.join(cwd, 'src/data/external/')
        tsv_cat_path = external_dir + 'currentcat.tsv'
        payload_cat_path = external_dir + 'payloadcat.tsv'

        urls = {tsv_cat_path: 'http://planet4589.org/space/gcat/tsv/derived/currentcat.tsv',
            payload_cat_path: 'http://planet4589.org/space/gcat/tsv/cat/psatcat.tsv'}

        for path, url in urls.items():
            self.DownloadJSRCatalogueIfNewer(path, url)

    def PullCatalogueSpaceTrack(self):
        #TODO: I Don't think this only pulls active satellites. It pulls all satellites?
        """
        Pull down the latest Active Satellites from SpaceTrack
        
        ### Exports
        - Text file of the latest active satellites
            - This is saved as src/data/external/spacetrack_active.txt 
        """
        cwd = os.getcwd()
        external_dir = os.path.join(cwd, 'src/data/external/')
        spacetrack_path = external_dir + f'spacetrack_active.txt'
        url = 'https://www.space-track.org/basicspacedata/query/class/gp/EPOCH/%3Enow-30/orderby/NORAD_CAT_ID,EPOCH/format/3le' # all satellites
        r = requests.get(url, stream=True)
        if r.status_code != 200: #TODO: does this need to be 200?
            # check if the file exists already, if so continue
            if os.path.exists(spacetrack_path):
                print("Failed to get file but a version of the file already exists locally, continuing")
                return
            else:
                raise Exception("Failed to download file")
        with open(spacetrack_path, "wb") as f:
            for chunk in r.iter_content(chunk_size=1024):
                f.write(chunk)

    # def PullCatalogueSpaceTrackAll(self):

    # ----------- REWRITE -----------
    # ----------- THIS DOES NOT PULL ALL THE SPACE TRACK CATALOG. IT PULLS ONLY STARLINKS. -----------

    #     """ Pull down the entire SpaceTrack catalogue """
    #     # location of file
    #     cwd = os.getcwd()
    #     external_dir = os.path.join(cwd, 'src/data/external/')
    #     spacetrack_path = external_dir + f'spacetrack_all.json'

    #     # this should change
    #     if os.path.exists(spacetrack_path):
    #         print("File already exists, not pulling an updated version of spacetrack's entire catalogue")
    #         return

    #     # Will require spacetrak login
    #     load_dotenv()
    #     username = os.getenv('SPACETRACK_USERNAME')
    #     password = os.getenv('SPACETRACK_PASSWORD')
    #     output = os.getenv('SPACETRACK_OUTPUT')
    #     siteCred = {'identity': username, 'password': password}

    #     # required urls
    #     uriBase                = "https://www.space-track.org"
    #     requestLogin           = "/ajaxauth/login"
    #     requestCmdAction       = "/basicspacedata/query" 
    #     requestFindStarlinks   = "/class/tle_latest/NORAD_CAT_ID/>40000/ORDINAL/1/OBJECT_NAME/STARLINK~~/format/json/orderby/NORAD_CAT_ID%20asc"
    #     requestOMMStarlink1    = "/class/omm/NORAD_CAT_ID/"
    #     requestOMMStarlink2    = "/orderby/EPOCH%20asc/format/json"
    #     entire_catalogue = '/class/gp/orderby/NORAD_CAT_ID asc/emptyresult/show'


    #     # use requests package to drive the RESTful session with space-track.org
    #     with requests.Session() as session:
    #         # run the session in a with block to force session to close if we exit
    #         # need to log in first. note that we get a 200 to say the web site got the data, not that we are logged in
    #         resp = session.post(uriBase + requestLogin, data = siteCred)
    #         if resp.status_code != 200:
    #             raise Exception(resp, "POST fail on login, please check credentials in .env file")

    #         # this query picks up all objects from the catalog. Note - a 401 failure shows you have bad credentials 
    #         resp = session.get(uriBase + requestCmdAction + entire_catalogue)
    #         if resp.status_code != 200:
    #             print(resp)
    #             raise Exception(resp, "GET fail on request for Starlink satellites")

    #         # use the json package to break the json formatted response text into a Python structure (a list of dictionaries)
    #         retData = json.loads(resp.text)
    #         satCount = len(retData)
    #         print(satCount, "satellites found from spacetrack")

    #         # convert json to string
    #         retDataStr = json.dumps(retData)

    #         # Open a file and write the JSON string to it
    #         with open(spacetrack_path, 'w') as f:
    #             f.write(retDataStr)
    #     session.close()

if __name__ == '__main__':
    with open(os.path.join(os.getcwd(), 'src/data/prediction_csv/sim_settings.json'), 'r') as f:
        json_data = f.read()

    # Parse the JSON string into a Python object
    policy = json.loads(json_data)
    catalogue = SpaceCatalogue(policy["sim_object_type"], policy["sim_object_catalogue"])
    catalogue.CreateCatalogueAll()
    # catalogue.Catalogue2SpaceObjects()

