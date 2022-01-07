
import os
from multiprocessing import Pool

import dustmaps.sfd
import dustmaps.bayestar
import numpy as np
import pandas as pd
import progressbar
import requests
from lightkurve import search_lightcurvefile
from lightkurve import LightCurveFileCollection

LC_DATA_PATH = 'lc_data/KIC-{}.csv'
FLARE_DATA_PATH = '../data_files/filtered_flares.csv'

def download_light_curve(KIC_ID):
    """
    Downloads and saves the light curve for the given KIC ID
    Args:
        KIC_ID (int): Kepler Input Catalogue ID.
    """

    KIC = "KIC {}".format(KIC_ID)
    if not os.path.isfile(LC_DATA_PATH.format(KIC_ID)):
        lcf = search_lightcurvefile(KIC, mission = 'Kepler').download_all()
        new_lcf =  LightCurveFileCollection([x for x in lcf if x.targetid == int(KIC_ID)])
        lc = new_lcf.SAP_FLUX.stitch()
        lc.to_csv(LC_DATA_PATH.format(KIC_ID))

def download_dust_maps():
    """
    Fetches the sfd and bayestar dust maps
    """
    dustmaps.sfd.fetch()
    dustmaps.bayestar.fetch()


def download_kic():
    dir_path = '../data_files'
    if not os.path.exists(dir_path):
        raise ValueError('{} is not found, run from the project root or create a folder'.format(dir_path))
    file_path = os.path.join(dir_path, 'kepler_kic_v10.csv.gz')
    url = 'https://archive.stsci.edu/pub/kepler/catalogs/kepler_kic_v10.csv.gz'
    chunk_size = 8192
    with requests.get(url, stream=True) as resp:
        resp.raise_for_status()
        size = int(resp.headers['content-length'])
        if os.path.exists(file_path) and os.stat(file_path).st_size == size:
            return
        with open(file_path, 'wb') as fh, progressbar.ProgressBar(max_value=size) as bar:
            for i, chunk in enumerate(resp.iter_content(chunk_size=chunk_size)):
                fh.write(chunk)
                bar.update(i * chunk_size)

def download_all():
    # 1. Downloading light curves and storing them
    print('Downloading light curves')
    df = pd.read_csv(FLARE_DATA_PATH)
    kic_id_array = np.unique(df['KIC'])
    index = 0

    with progressbar.ProgressBar(max_value = len(kic_id_array)) as bar:
            for kic_id in kic_id_array:
                download_light_curve(kic_id)
                index += 1
                bar.update(index)

    # 2. Downloading dust maps and storing them
    print('Downloading dust maps')
    download_dust_maps()

    # 3. Downloading KIC file if it doesn't exist
    print('Downloading KIC')
    download_kic()

    print('Data downloaded!')