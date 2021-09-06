import os.path
from functools import lru_cache

import lightkurve as lk
import light_curve as lc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import astropy
from lightkurve import search_lightcurvefile
from lightkurve import LightCurveFileCollection

LC_DATA_PATH = 'lc_data/KIC-{}.csv'
FLARE_DATA_PATH = 'data_files/apjaa8ea2t3_mrt.txt'
FLARE_INSTANCES_PATH = 'flare_instances/KIC-{}/'
FLARE_INSTANCE = 'flare_instances/KIC-{kic}/{start}-{end}.csv'

def extract_features_in_lsst_passbands(lc_dict):
    for passband in lc_dict.keys():
        print(passband)
        get_lc_features(lc_dict[passband])

def get_lc_features(lcf):

    t = lcf.time
    m = lcf.flux
    sigma = None

    # Removing nan values
    index = np.logical_not(np.isnan(m))
    m = m[index]
    t = t[index] 

    # Moment statistics
    std = lc.StandardDeviation()
    skew = lc.Skew()
    kurtosis = lc.Kurtosis()

    # Shape based features
    etaE = lc.EtaE()
    stetsonK = lc.StetsonK()


    # Statistical features
    amplitude = lc.Amplitude()
    beyond_1_std = lc.BeyondNStd(1)
    beyond_2_std = lc.BeyondNStd(2)
    beyond_3_std = lc.BeyondNStd(3)
    magnitude_percentile_ratio = lc.MagnitudePercentageRatio(0.4, 0.05)
    cusum = lc.Cusum()
    
    extr = lc.Extractor(std, skew, kurtosis, etaE, stetsonK, amplitude, beyond_1_std, beyond_2_std, beyond_3_std, magnitude_percentile_ratio, cusum)

    print(f'Half-amplitude is {amplitude(t, m, sigma).item()}')
    print(f'Fraction of observations beyond 2 std from mean is {beyond_2_std(t, m, sigma).item()}')
    print(f'Eta E fit is {etaE(t, m, sigma).item()}')
    print(f'All features: {extr(t, m, sigma)}')

if __name__=="__main__":
    get_features()

