import astropy
from astropy.io import ascii
from astropy.modeling import models
from astropy import units as u
import lightkurve as lk
from lightkurve import search_lightcurvefile
from lightkurve import LightCurveFileCollection
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import os.path
import pandas as pd
import csv
import numpy as np
import math
from lc_tools import *

def get_spectra_data(lc, KIC_ID, temperature):
    # Reading flare data file.
    flare_data = astropy.io.ascii.read(FLARE_DATA_PATH, quotechar="\s")
    bb_model = models.BlackBody(temperature=temperature*u.K) 

    wavelenghts = {
        'kepler_min': 4183.66,
        'kepler_max': 9050.23,
        'lsst_u_min': 3060.00,
        'lsst_u_max': 4080.00,
        'lsst_g_min': 3869.99,
        'lsst_g_max': 5665.01,
        'lsst_r_min': 5375.02,
        'lsst_r_max': 7054.98,
        'lsst_i_min': 6764.99,
        'lsst_i_max': 8325.01,
        'lsst_z_min': 8034.99,
        'lsst_z_max': 9380.01,
        'lsst_y_min': 9090.18,
        'lsst_y_max': 10999.99,

    }
    #Plotting points identified as flares.
    for flare in flare_data:
        if flare['KIC'] != int(KIC_ID):
            continue
            
        # Obtaining start and end indices for flare instances
        flare_lc = get_flare_lc(lc,flare)

        kepler_flux_mean = compute_mean_flux(wavelenghts['kepler_min'], wavelenghts['kepler_max'], bb_model)

        lsst_flux_u_mean = compute_mean_flux(wavelenghts['lsst_u_min'], wavelenghts['lsst_u_max'], bb_model)
        lsst_flux_g_mean = compute_mean_flux(wavelenghts['lsst_g_min'], wavelenghts['lsst_g_max'], bb_model)
        lsst_flux_r_mean = compute_mean_flux(wavelenghts['lsst_r_min'], wavelenghts['lsst_r_max'], bb_model)
        lsst_flux_i_mean = compute_mean_flux(wavelenghts['lsst_i_min'], wavelenghts['lsst_i_max'], bb_model)
        lsst_flux_z_mean = compute_mean_flux(wavelenghts['lsst_z_min'], wavelenghts['lsst_z_max'], bb_model)
        lsst_flux_y_mean = compute_mean_flux(wavelenghts['lsst_y_min'], wavelenghts['lsst_y_max'], bb_model)

        percentage_change_u = (lsst_flux_u_mean - kepler_flux_mean) / kepler_flux_mean
        percentage_change_g = (lsst_flux_g_mean - kepler_flux_mean) / kepler_flux_mean
        percentage_change_r = (lsst_flux_r_mean - kepler_flux_mean) / kepler_flux_mean
        percentage_change_i = (lsst_flux_i_mean - kepler_flux_mean) / kepler_flux_mean
        percentage_change_z = (lsst_flux_z_mean - kepler_flux_mean) / kepler_flux_mean
        percentage_change_y = (lsst_flux_y_mean - kepler_flux_mean) / kepler_flux_mean
        
        u_band = [(1 + percentage_change_u) * flux for flux in flare_lc.flux]
        g_band = [(1 + percentage_change_g) * flux for flux in flare_lc.flux]
        r_band = [(1 + percentage_change_r) * flux for flux in flare_lc.flux]
        i_band = [(1 + percentage_change_i) * flux for flux in flare_lc.flux]
        z_band = [(1 + percentage_change_z) * flux for flux in flare_lc.flux]
        y_band = [(1 + percentage_change_y) * flux for flux in flare_lc.flux]

        plt.plot(flare_lc.time, flare_lc.flux, color = 'black', label = 'Kepler', linestyle='-', marker='.')
        plt.plot(flare_lc.time, u_band, color = 'purple', label = 'Lsst u band', linestyle='-', marker='.')
        plt.plot(flare_lc.time, g_band, color = 'green', label = 'Lsst g band', linestyle='-', marker='.')
        plt.plot(flare_lc.time, r_band, color = 'orange', label = 'Lsst r band', linestyle='-', marker='.')
        plt.plot(flare_lc.time, i_band, color = 'red', label = 'Lsst i band', linestyle='-', marker='.')
        plt.plot(flare_lc.time, z_band, color = 'grey', label = 'Lsst z band', linestyle='-', marker='.')
        plt.plot(flare_lc.time, y_band, color = 'brown', label = 'Lsst y band', linestyle='-', marker='.')
        
        plt.legend()
        plt.show()

def compute_mean_flux(min_wav_len, max_wav_len, bb_model):
    
    wavlengths = np.arange(min_wav_len,	max_wav_len) * u.AA
    kepler_flux = bb_model(wavlengths)
    kepler_flux_mean = np.mean(kepler_flux)
    return kepler_flux_mean