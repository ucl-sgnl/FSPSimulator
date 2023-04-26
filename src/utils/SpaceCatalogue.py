import os
import json
import requests
import wget
import numpy as np
import pandas as pd
from dotenv import load_dotenv
import datetime
from src.utils.SpaceObject import SpaceObject
from src.utils.coords import tle_parse

class SpaceCatalogue:
    def __init__(self, sim_object_type, sim_object_catalogue):
        self.PullAllCataloguesIfNewer()
        self.Satellites = []
        self.Catalogue = []
        self.CurrentCatalogue = None
        self.sim_object_type = sim_object_type # this will tell you whether to which part of the catalogue to use, active, inactive or all
        self.sim_object_catalogue = sim_object_catalogue

    def ReturnCatalogue(self):
        """
        This will return the current space catalogue. 

        Warning: Due to the policy configuration, it may not have merged JSR and Celestrak.

        ### Returns
        - List of SpaceObjects
        """
        return self.Catalogue

    def CreateCatalogueActive(self):
        """
        This function will merge the JSR and Celestrak catalogues and create a list of SpaceObjects 
        
        ### Exports
        - List of merged space objects to 'src/data/external/active_jsr_celestrak.csv'
        """
        # firstly define which combination of catalogues they would like to use
        jsr_cat = pd.read_csv(os.path.join(os.getcwd(), 'src/data/external/currentcat.tsv'), sep='\t')
        jsr_cat_extra_info = pd.read_csv(os.path.join(os.getcwd(), 'src/data/external/satcat.tsv'), sep='\t')
        with open(os.path.join(os.getcwd(), 'src/data/external/celestrak_active.txt'), 'r') as f:
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
        with open(os.path.join(os.getcwd(), 'src/data/catalogue/active_jsr_celestrak_latest.csv'), 'w') as f:
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
        spacetrack = pd.read_json(os.path.join(os.getcwd(), 'src/data/external/celestrak_all.json'))
        jsr_cat_extra_info = pd.read_csv(os.path.join(os.getcwd(), 'src/data/external/satcat.tsv'), sep='\t')

        # merge the two dataframes, keeping all items in spacetrack dataframe
        self.CurrentCatalogueDF = spacetrack.merge(jsr_cat_extra_info, right_on='Piece', left_on='OBJECT_ID', how='left')

        # also export the catalogue for audit purposes
        self.CurrentCatalogueDF.to_csv(os.path.join(os.getcwd(), 'src/data/catalogue/All_catalogue_latest.csv'), index=False)
        
    def Catalogue2SpaceObjects(self):
        """
        This function will convert the current catalogue into a list of SpaceObjects

        ### Returns
        - List of SpaceObjects
            - View the SpaceObject class for more information
        """
        # headings available:
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
        if self.sim_object_type == "active":
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
                                                        apogee_altitude=row['Apogee'],
                                                        tle=row['TLE']
                                                    ))
        elif self.sim_object_type == "all":
            for index, row in self.CurrentCatalogueDF.iterrows():
                tle = row["TLE_LINE0"] + "\n" + row["TLE_LINE1"] + "\n" + row["TLE_LINE2"]
                self.Catalogue.append(SpaceObject( object_type=row['OBJECT_TYPE'], 
                                                        #payload_operational_status=row['Active'], ## this will need to be updated to celestrak's active status
                                                        application="Unknown", 
                                                        #operator=row['Owner'], 
                                                        mass=row['Mass'], 
                                                        launch_site=row["SITE"],
                                                        # type = row['Type'],
                                                        # maneuverable=maneuverable, 
                                                        # spin_stabilized=spin_stabilized, 
                                                        # radar_cross_section=radar_cs, 
                                                        # characteristic_area=area, 
                                                        # characteristic_length=length, 
                                                        # propulsion_type=propulsion, 
                                                        sma=row['SEMIMAJOR_AXIS'], 
                                                        eccentricity=row['ECCENTRICITY'], 
                                                        inc=row['INCLINATION'], 
                                                        argp=row['ARG_OF_PERICENTER'], 
                                                        raan=row['RA_OF_ASC_NODE'], 
                                                        # launch_site=launch_site,
                                                        source = row['COUNTRY_CODE'],
                                                        launch_date=row['LAUNCH_DATE'], 
                                                        decay_date=row['DECAY_DATE'], 
                                                        cospar_id=row['OBJECT_ID'],
                                                        rso_name=row['OBJECT_NAME'],
                                                        perigee_altitude=row['PERIAPSIS'],
                                                        apogee_altitude=row['APOAPSIS'],
                                                        tle=tle,
                                                        epoch=row['EPOCH']
                                                    ))
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

        urls = {tsv_cat_path: 'https://planet4589.org/space/gcat/tsv/derived/currentcat.tsv',
            payload_cat_path: 'https://planet4589.org/space/gcat/tsv/cat/psatcat.tsv'}

        for path, url in urls.items():
            self.DownloadJSRCatalogueIfNewer(path, url)

    def PullCatalogueCelestrakActive(self):
        """
        Pull down the latest Active Satellites from Celestrak
        
        ### Exports
        - Text file of the latest active satellites
            - This is saved as src/data/external/celestrak_active.txt 
        """
        cwd = os.getcwd()
        external_dir = os.path.join(cwd, 'src/data/external/')
        celestrak_path = external_dir + f'celestrak_active.txt'
        # url = 'http://celestrak.org/NORAD/elements/gp.php?GROUP=active&FORMAT=tle' # just active satellites
        url = 'https://www.space-track.org/basicspacedata/query/class/gp/EPOCH/%3Enow-30/orderby/NORAD_CAT_ID,EPOCH/format/3le' # all satellites
        r = requests.get(url, stream=True)
        if r.status_code != 200:
            # check if the file exists already, if so continue
            if os.path.exists(celestrak_path):
                print("File already exists, continuing")
                return
            else:
                raise Exception("Failed to download file")
        with open(celestrak_path, "wb") as f:
            for chunk in r.iter_content(chunk_size=1024):
                f.write(chunk)

    def PullCatalogueSpaceTrackAll(self):
        """ Pull down the entire Celestrak catalogue """
        # location of file
        cwd = os.getcwd()
        external_dir = os.path.join(cwd, 'src/data/external/')
        celestrak_path = external_dir + f'celestrak_all.json'

        # this should change
        if os.path.exists(celestrak_path):
            print("File already exists, not pulling an updated version of celestrak's entire catalogue")
            return

        # Will require spacetrak login
        load_dotenv()
        username = os.getenv('SPACETRACK_USERNAME')
        password = os.getenv('SPACETRACK_PASSWORD')
        output = os.getenv('SPACETRACK_OUTPUT')
        siteCred = {'identity': username, 'password': password}

        # required urls
        uriBase                = "https://www.space-track.org"
        requestLogin           = "/ajaxauth/login"
        requestCmdAction       = "/basicspacedata/query" 
        requestFindStarlinks   = "/class/tle_latest/NORAD_CAT_ID/>40000/ORDINAL/1/OBJECT_NAME/STARLINK~~/format/json/orderby/NORAD_CAT_ID%20asc"
        requestOMMStarlink1    = "/class/omm/NORAD_CAT_ID/"
        requestOMMStarlink2    = "/orderby/EPOCH%20asc/format/json"
        entire_catalogue = '/class/gp/orderby/NORAD_CAT_ID asc/emptyresult/show'


        # use requests package to drive the RESTful session with space-track.org
        with requests.Session() as session:
            # run the session in a with block to force session to close if we exit
            # need to log in first. note that we get a 200 to say the web site got the data, not that we are logged in
            resp = session.post(uriBase + requestLogin, data = siteCred)
            if resp.status_code != 200:
                raise Exception(resp, "POST fail on login, please check credentials in .env file")

            # this query picks up all objects from the catalog. Note - a 401 failure shows you have bad credentials 
            resp = session.get(uriBase + requestCmdAction + entire_catalogue)
            if resp.status_code != 200:
                print(resp)
                raise Exception(resp, "GET fail on request for Starlink satellites")

            # use the json package to break the json formatted response text into a Python structure (a list of dictionaries)
            retData = json.loads(resp.text)
            satCount = len(retData)
            print(satCount, "satellites found from celestrak")

            # convert json to string
            retDataStr = json.dumps(retData)

            # Open a file and write the JSON string to it
            with open(celestrak_path, 'w') as f:
                f.write(retDataStr)
        session.close()

    def PullAllCataloguesIfNewer(self):
        """
        Checks for each catalogue that is available (currently, JSR and Celestrak) whether a newer version is available.

        If a newer version is available, it will pull it down and save it to the external directory.

        During testing, some of these may be commented out due to the time it takes to pull the data.
        """
        self.PullCatalogueJSR()
        self.PullCatalogueCelestrakActive()
        self.PullCatalogueSpaceTrackAll()
        # pass

if __name__ == '__main__':
    with open(os.path.join(os.getcwd(), 'src/data/prediction_csv/policy_fsptest.json'), 'r') as f:
        json_data = f.read()

    # Parse the JSON string into a Python object
    policy = json.loads(json_data)
    catalogue = SpaceCatalogue(policy["sim_object_type"], policy["sim_object_catalogue"])
    catalogue.CreateCatalogueAll()
    # catalogue.Catalogue2SpaceObjects()

