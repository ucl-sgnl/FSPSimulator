def update_catalogue_celestrak():
    print('Updating celestrak files')

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
    """Download a file from a URL if it is newer than the local file."""
    if os.path.exists(local_path):
        local_last_modified = datetime.datetime.fromtimestamp(os.path.getmtime(local_path))

        response = requests.head(url)
        remote_last_modified = datetime.datetime.strptime(response.headers['Last-Modified'], '%a, %d %b %Y %H:%M:%S %Z')

        if remote_last_modified <= local_last_modified:
            return

    wget.download(url, local_path)

def update_catalogue_jsr():
    """Update the JSR catalogue files."""
    cwd = os.getcwd()
    external_dir = os.path.join(cwd, 'src/data/external/')
    tsv_cat_path = external_dir + 'currentcat.tsv'
    payload_cat_path = external_dir + 'payloadcat.tsv'

    urls = {tsv_cat_path: 'http://planet4589.org/space/gcat/tsv/derived/currentcat.tsv',
        payload_cat_path: 'http://planet4589.org/space/gcat/tsv/cat/psatcat.tsv'}

    for path, url in urls.items():
        jsr_download_if_newer(path, url)

def update_catalogue_celestrak():
    """Pull down the latest Active Satellites from Celestrak"""
    cwd = os.getcwd()
    external_dir = os.path.join(cwd, 'src/data/external/')
    celestrak_path = external_dir + f'celestrak_active.txt'
    url = 'http://celestrak.org/NORAD/elements/gp.php?GROUP=active&FORMAT=tle'
    r = requests.get(url, stream=True)
    if r.status_code != 200:
        raise Exception("Failed to download file")
    with open(celestrak_path, "wb") as f:
        for chunk in r.iter_content(chunk_size=1024):
            f.write(chunk)
    
if __name__ == '__main__':
    update_catalogue_jsr()
    # update_catalogue_celestrak()
