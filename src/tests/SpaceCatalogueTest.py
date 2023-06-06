import unittest
import pandas as pd
import os.path
import json
from src.utils.SpaceCatalogue import SpaceCatalogue

class TestMyClass(unittest.TestCase):
    
    def setUp(self):
        with open(os.path.join(os.getcwd(), 'src/data/prediction_csv/sim_settings.json'), 'r') as f:
            json_data = f.read()
        policy = json.loads(json_data)
        self.SpaceCat = SpaceCatalogue(policy["sim_object_type"], policy["sim_object_catalogue"])
        
    def test_ReturnCatalogue(self):
        # Test if the function returns a list of SpaceObjects
        catalogue = self.SpaceCat.ReturnCatalogue()
        self.assertIsInstance(catalogue, list)
        
        # Test that the list is empty if no catalogue has been created
        self.assertEqual(len(catalogue), 0)

    
    def test_CreateCatalogueActive(self):
        # Test if the function exports a file
        self.SpaceCat.CreateCatalogueActive()
        self.assertTrue(os.path.exists(os.path.join(os.getcwd(), 'src/data/catalogue/active_jsr_spacetrack_latest.csv')))
       
        # Test if the file contains data
        df = pd.read_csv(os.path.join(os.getcwd(), 'src/data/catalogue/active_jsr_spacetrack_latest.csv'))
        self.assertGreater(len(df), 0)
        # Additional test cases:
        # Test if the exported file has the expected columns
        # Test if the data in the exported file is as expected
    
    def test_CreateCatalogueAll(self):
        # Test if the function exports a file
        self.SpaceCat.CreateCatalogueAll()
        self.assertTrue(os.path.exists(os.path.join(os.getcwd(), 'src/data/catalogue/All_catalogue_latest.csv')))
        # Test if the file contains data
        df = pd.read_csv(os.path.join(os.getcwd(), 'src/data/catalogue/All_catalogue_latest.csv'))
        self.assertGreater(len(df), 0)
        # Additional test cases:
        # Test if the exported file has the expected columns
        # Test if the data in the exported file is as expected
