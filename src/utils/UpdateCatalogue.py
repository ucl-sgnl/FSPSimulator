import os
import datetime
import requests
import wget
import numpy as np
import math
import pandas as pd

cwd_path = os.getcwd()

def read_csv(file_path):
    df = pd.read_csv(file_path, sep=',')
    return df

def jsr_download_if_newer(local_path, url):
    # function to download a file from a url if the local file is older than the remote file
    # avoids having to re-download the file unnecessarily
    if os.path.exists(local_path):
        local_last_modified = datetime.datetime.fromtimestamp(os.path.getmtime(local_path))

        response = requests.head(url)
        remote_last_modified = datetime.datetime.strptime(response.headers['Last-Modified'], '%a, %d %b %Y %H:%M:%S %Z')

        if remote_last_modified <= local_last_modified:
            print(f"{local_path} is already up to date.")
            return

    wget.download(url, local_path)
    print(f"{local_path} has been updated.")

def update_catalogue_jsr():
    # set paths
    cwd = os.getcwd()
    external_dir = os.path.join(cwd, 'src/data/external/')
    tsv_cat_path = external_dir + 'currentcat.tsv'
    payload_cat_path = external_dir + 'payloadcat.tsv'

    # current cat is the general catalog
    # payload cat is the payload specific catalog and contains detailed info about payloads (e.g. when they are active)
    urls = {tsv_cat_path: 'https://planet4589.org/space/gcat/tsv/derived/currentcat.tsv',
        payload_cat_path: 'https://planet4589.org/space/gcat/tsv/cat/psatcat.tsv'}

    for path, url in urls.items():
        jsr_download_if_newer(path, url)

def update_catalogue_celestrak():
    print('Updating celestrak files')

if __name__ == '__main__':
    update_catalogue_jsr()
    # update_catalogue_celestrak()
