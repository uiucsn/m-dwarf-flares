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

labels = ['StdDeviation', 'Skew', 'Kurtosis', 'EtaE', 'StetsonK', 'Amplitude', 'Beyond1Std', 'Beyond2Std', 'Beyond3Std', 'MagnitudePercentageRatio', 'Cusum']
tables = {
    'u': pd.DataFrame(columns=labels),
    'g': pd.DataFrame(columns=labels),
    'r': pd.DataFrame(columns=labels),
    'i': pd.DataFrame(columns=labels),
    'z': pd.DataFrame(columns=labels),
    'y': pd.DataFrame(columns=labels),
    'kep': pd.DataFrame(columns=labels),
}

def extract_features_in_lsst_passbands(lc_dict, index):
    for passband in lc_dict.keys():
        features = get_lc_features(lc_dict[passband])
        tables[passband].loc[index] = features
        tables[passband].to_csv('simulation_features/{}.csv'.format(passband))


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
    return extr(t, m, sigma)

def plot_all_feature_distributions():
    for passband in ['u','g','r','i','z','y','kep']:
        tables[passband] = pd.read_csv('simulation_features/{}.csv'.format(passband))
    
    for column in tables['u'].columns:
        all_data = tables['u'][column] + tables['g'][column] + tables['r'][column] + tables['i'][column] + tables['z'][column] + tables['y'][column]

        # max_amp = max(all_data)
        # min_amp = min(all_data)
        # bin_width = (max_amp - min_amp) / 20
        # bins = np.arange(min_amp, max_amp + bin_width, bin_width)

        fig, ax = plt.subplots(2, 3, constrained_layout=True,
                        sharex=True, sharey=True)

        fig.suptitle('Distribution of {0} in LSST passbands. x axis = {0}'.format(column))

        ax[0][0].hist(tables['u'][column], histtype='step', color='m', label = 'u band') 
        ax[0][1].hist(tables['g'][column], histtype='step', color='g', label = 'g band') 
        ax[0][2].hist(tables['r'][column], histtype='step', color='r', label = 'r band') 
        ax[1][0].hist(tables['i'][column], histtype='step', color='c', label = 'i band') 
        ax[1][1].hist(tables['z'][column], histtype='step', color='b', label = 'z band') 
        ax[1][2].hist(tables['y'][column], histtype='step', color='y', label = 'y band') 

        ax[0][0].legend(loc='upper right')
        ax[0][1].legend(loc='upper right')
        ax[0][2].legend(loc='upper right')
        ax[1][0].legend(loc='upper right')
        ax[1][1].legend(loc='upper right')
        ax[1][2].legend(loc='upper right')

        plt.xlabel(column)
        plt.ylabel("Number of m dwarfs")
        
        
        plt.savefig("feature_distribution/{}_distribution.pdf".format(column))

if __name__ == "__main__":
    plot_all_feature_distributions()
