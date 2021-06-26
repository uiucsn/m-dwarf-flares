from multiprocessing import Pool

import dustmaps.sfd
import dustmaps.bayestar
import lightkurve as lk
import numpy as np
import pandas as pd
import progressbar
from lightkurve import search_lightcurvefile
from lightkurve import LightCurveFileCollection

LC_DATA_PATH = 'lc_data/KIC-{}.csv'
FLARE_DATA_PATH = 'data_files/filtered_flares.csv'

def download_light_curve(KIC_ID):
    """
    Downloads and saves the light curve for the given KIC ID

    Args:
        KIC_ID (int): Kepler Input Catalogue ID.
    """
    KIC = "KIC {}".format(KIC_ID)
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


def main():
    # 1. Downloading light curves and storing them
    print('Downloading light curves')
    df = pd.read_csv(FLARE_DATA_PATH)
    kic_id_array = np.unique(df['KIC'])
    index = 0

    with progressbar.ProgressBar(max_value = len(kic_id_array)) as bar:
        with Pool(processes=16) as pool:
            for index, _ in enumerate(pool.imap(download_light_curve, kic_id_array)):
                index += 1
                bar.update(index)

    # 2. Downloading dust maps and storing them
    print('Downloading dust maps')
    download_dust_maps()

    print('Data downloaded!')


if __name__ == '__main__':
    main()
