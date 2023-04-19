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
            temp_satellite = Satellite()
            temp_satellite.cospar_id = row['JCAT']
            temp_satellite.rso_name = row['Name']
            temp_satellite.rso_type = row['Type']
            temp_satellite.source = row['Operator']
            temp_satellite.decay_date = row['DDate']
            temp_satellite.orbital_period = row['Period']
            temp_satellite.perigee_altitude = row['Perigee']
            temp_satellite.apogee_altitude = row['Apogee']
            temp_satellite.inc = row['Inc']
            temp_satellite.ArgPerigee = row['ArgPerigee']
            temp_satellite.MeanAnomaly = row['MeanAnomaly']
            temp_satellite.MeanMotion = row['MeanMotion']
            temp_satellite.Revolution = row['Revolution']
            temp_satellite.LaunchDate = row['LaunchDate']
            temp_satellite.LaunchSite = row['LaunchSite']
            temp_satellite.LaunchVehicle = row['LaunchVehicle']
       


            