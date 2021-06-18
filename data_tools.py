from astropy.io import fits
import astropy
import numpy as np
import pandas as pd
from astropy.table import QTable
from pyvo.dal import TAPService
import lightkurve as lk
import matplotlib.pyplot as plt
from lc_tools import load_light_curve, get_flare_lc_from_time, get_normalized_lc


def save_effective_kepler_temps():
    """
    Function to save the effective temp for the the 541 kic objects from the kepler input catalogue
    """
    mag = pd.read_csv('data_files/mag.csv')
    kepler_cataloug = pd.read_csv('data_files/kepler_kic_v10.csv.gz')
    
    df = pd.DataFrame()
    k = kepler_cataloug.loc[kepler_cataloug['kic_kepler_id'].isin(mag['KIC ID'].astype(int))]

    k.to_csv('data_files/eff_temp.csv')


def remove_incomplete_entries_from_flare_data():
    """
    Removes flares for which Gaia distance data does not exist
    """
    dist = pd.read_csv('data_files/dist_new.csv')
    flare_data = astropy.io.ascii.read('data_files/apjaa8ea2t3_mrt.txt', quotechar="\s")
    df = flare_data.to_pandas()

    filtered = df.loc[df['KIC'].isin(dist['KIC ID'].astype(int))]
    filtered.to_csv('data_files/filtered_flares.csv')

def fetch_Gaia_Data():
    """
    Function to query Gaia Data
    """

    FLARE_DATA_PATH = 'data_files/apjaa8ea2t3_mrt.txt'

    flare_data = astropy.io.ascii.read(FLARE_DATA_PATH, quotechar="\s")
    kic = np.array(flare_data['KIC'])
    list = np.unique(kic)

    table = QTable(names = ('KIC ID' ,'d', 'r_med_geo',  'r_lo_geo' , 'r_hi_geo', 'r_med_photogeo', 'r_lo_photogeo', 'r_hi_photogeo', 'phot_g_mean_mag'),
    dtype= ('int32', 'float64', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32'))
    for KIC_ID in list:

        KIC = "KIC {}".format(KIC_ID)
        lcf = lk.search_lightcurvefile(KIC_ID, mission = 'Kepler')

        ra_deg = lcf.ra
        dec_deg = lcf.dec

        search_radius_arcsec = 1
        string_query = f'''SELECT distance(ra, dec, {ra_deg[0]}, {dec_deg[0]}) as d, r_med_geo, r_lo_geo, r_hi_geo, r_med_photogeo,r_lo_photogeo, r_hi_photogeo, phot_g_mean_mag FROM gedr3dist.main JOIN gaia.edr3lite USING (source_id) WHERE distance(ra, dec, {ra_deg[0]}, {dec_deg[0]}) < {search_radius_arcsec / 3600.0}'''
        
        tap = TAPService('https://dc.zah.uni-heidelberg.de/tap')
        response = tap.search(string_query)
        
        r = response.to_table()

        if len(r) > 1:
            index = np.where(r['r_med_geo'] == np.amin(r['r_med_geo']))
            table.add_row((KIC_ID, r[index]['d'], r[index]['r_med_geo'], r[index]['r_lo_geo'], r[index]['r_hi_geo'], r[index]['r_med_photogeo'], r[index]['r_lo_photogeo'], r[index]['r_hi_photogeo'], r[index]['phot_g_mean_mag']))
        elif len(r) == 0:
            continue
        else:
            table.add_row((KIC_ID, r[0]['d'], r[0]['r_med_geo'], r[0]['r_lo_geo'], r[0]['r_hi_geo'], r[0]['r_med_photogeo'], r[0]['r_lo_photogeo'], r[0]['r_hi_photogeo'], r[0]['phot_g_mean_mag']))
        table.write('data_files/dist_new.csv', format = 'ascii.csv')

        print(table)

 
def save_flare_flux_amps():
    """
    Function to save flare's relative flux ampslitude in the filtered flares file
    """
    df = pd.read_csv('data_files/filtered_flares.csv')
    print(df)
    amps = []
    for i in range(len(df)):
        lc = load_light_curve(df['KIC'][i])
        flare_lc = get_flare_lc_from_time(lc, df['St-BKJD'][i], df['End-BKJD'][i])
        normalized_flare = get_normalized_lc(flare_lc)
        amps.append(np.amax(normalized_flare.flux))
        print((i / len(df)) * 100, '%')
    
    df['flux_amp'] = amps
    df.to_csv('data_files/filtered_flares.csv')

def filter_out_flares_with_nans():
    """
    Removes all flares that contain nan values in their flux arrays
    """
    df = pd.read_csv('data_files/filtered_flares.csv')
    count = 0
    indices = []
    for i in range(len(df)):
        lc = load_light_curve(df['KIC'][i])
        flare_lc = get_flare_lc_from_time(lc, df['St-BKJD'][i], df['End-BKJD'][i])
        print((i / len(df)) * 100, '% Number of bad flares', count, end='\r')
        if np.isnan(flare_lc.flux).any():
            print('bad flare', i)
            indices.append(i)
            count += 1
    df.drop(df.index[indices], inplace=True)
    print(len(df))
    df.to_csv('data_files/filtered_flares.csv')
