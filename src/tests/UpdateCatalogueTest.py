import os
import unittest
from unittest.mock import patch, mock_open
import datetime
import requests
import wget
import pandas as pd
from src.utils.UpdateCatalogue import read_csv, jsr_download_if_newer, update_catalogue_jsr

class TestYourModule(unittest.TestCase):
    def test_read_csv(self):
        data = 'column1,column2\n1,2\n3,4'
        with patch('src.utils.UpdateCatalogue.pd.read_csv', return_value=pd.DataFrame([[1, 2], [3, 4]])) as mocked_read_csv:
            with patch('builtins.open', new=mock_open(read_data=data)):
                file_path = 'dummy.csv'
                df = read_csv(file_path)
                mocked_read_csv.assert_called_once_with(file_path, sep=',')
                self.assertEqual(len(df), 2)
                self.assertEqual(df.at[0, 'column1'], 1)

    def test_jsr_download_if_newer_local_newer(self):
        local_path = 'dummy_path'
        url = 'https://dummy_url'
        with patch('src.utils.UpdateCatalogue.os.path.exists', return_value=True), \
             patch('src.utils.UpdateCatalogue.os.path.getmtime', return_value=datetime.datetime.now().timestamp()), \
             patch('src.utils.UpdateCatalogue.requests.head', return_value=requests.Response()), \
             patch('src.utils.UpdateCatalogue.datetime.datetime.strptime', return_value=datetime.datetime.now() - datetime.timedelta(days=1)), \
             patch('src.utils.UpdateCatalogue.wget.download') as mocked_wget:
            
            jsr_download_if_newer(local_path, url)
            mocked_wget.assert_not_called()

    def test_jsr_download_if_newer_local_older(self):
        local_path = 'dummy_path'
        url = 'https://dummy_url'
        with patch('src.utils.UpdateCatalogue.os.path.exists', return_value=True), \
             patch('src.utils.UpdateCatalogue.os.path.getmtime', return_value=(datetime.datetime.now() - datetime.timedelta(days=2)).timestamp()), \
             patch('src.utils.UpdateCatalogue.requests.head', return_value=requests.Response()), \
             patch('src.utils.UpdateCatalogue.datetime.datetime.strptime', return_value=datetime.datetime.now() - datetime.timedelta(days=1)), \
             patch('src.utils.UpdateCatalogue.wget.download') as mocked_wget:
            
            jsr_download_if_newer(local_path, url)
            mocked_wget.assert_called_once_with(url, local_path)

    def test_update_catalogue_jsr(self):
        with patch('src.utils.UpdateCatalogue.os.getcwd', return_value='dummy_cwd'), \
             patch('src.utils.UpdateCatalogue.os.path.join', return_value='dummy_cwd/src/data/external/'), \
             patch('src.utils.UpdateCatalogue.jsr_download_if_newer') as mocked_jsr_download:
            
            update_catalogue_jsr()
            self.assertEqual(mocked_jsr_download.call_count, 2)

if __name__ == '__main__':
    unittest.main()
