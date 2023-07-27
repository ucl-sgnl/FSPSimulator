import json
import os
import unittest
from src.api import api

class TestAPI(unittest.TestCase):
    def setUp(self):
        api.testing = True
        self.client = api.test_client()

    def test_hello(self):
        response = self.client.get('/')
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.get_data(as_text=True), 'Hello, welcome to the UCL Future Space Populations API!')

    def test_post_json(self):
        data = {
            'id': 'test_id',
            'parameter1': 1,
            'parameter2': 2
        }
        response = self.client.post('/newsim', json=data)
        self.assertEqual(response.status_code, 201)
        self.assertEqual(json.loads(response.get_data(as_text=True)), {'success': True})
        filename = os.path.join(os.getcwd(), 'src', 'data', 'policy', 'test_id.json')
        with open(filename) as f:
            self.assertEqual(json.load(f), data)

    def test_get_file(self):
        data = {
            'id': 'test_id',
            'parameter1': 1,
            'parameter2': 2
        }
        filename = os.path.join(os.getcwd(), 'src', 'data', 'policy', 'test_id.json')
        with open(filename, 'w') as f:
            json.dump(data, f)
        response = self.client.get('/getfile?type=policy&id=test_id')
        self.assertEqual(response.status_code, 200)
        self.assertEqual(json.loads(response.get_data(as_text=True)), data)

    def test_update_external_cat(self):
        response = self.client.get('/updatecatalogue?type=celestrak')
        self.assertEqual(response.status_code, 201)
        self.assertEqual(json.loads(response.get_data(as_text=True)), {'success': True})

        # Add more tests for update_catalogue_jsr and invalid type parameters

if __name__ == '__main__':
    unittest.main()
