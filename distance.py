import time
import concurrent.futures
import matplotlib.pyplot as plt
import pickle
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
from spectra_tools import *

dist_data = astropy.io.ascii.read('dist.csv')
mag_data = astropy.io.ascii.read('mag.csv')

def get_luminosity_with_magnitude(KIC):

    for row in dist_data:
        if row['KIC ID'] == int(KIC):
            if row['r_med_photogeo'] != 0:
                distance = row['r_med_photogeo']
            else:
                distance = row['r_med_geo']
    
    for row in mag_data:
        if row['KIC ID'] == int(KIC):
            ab_magnitude = row['kic_kepmag']


    flux = (10 ** (-ab_magnitude / 2.5)) * 3630.7805 * astropy.units.Jy
    luminosity = flux * 4 * math.pi * (distance * astropy.units.pc) ** 2
    return luminosity

def get_spectra_data_with_luminosity(flare_lc, KIC_ID, temp, luminosity):
    # Reading flare data file.
    T = temp * u.K
    lmbd, t = get_kepler_transmission()
    intensityKepler = simps(x=lmbd, y=bb(lmbd * u.AA, T).value * t) / simps(x=lmbd, y=t)

    intensity_u = compute_band_intensity("u", temp)
    intensity_g = compute_band_intensity("g", temp)
    intensity_r = compute_band_intensity("r", temp)
    intensity_i = compute_band_intensity("i", temp)
    intensity_z = compute_band_intensity("z", temp)
    intensity_y = compute_band_intensity("y", temp)

    relative_change_u = (intensity_u - intensityKepler)/intensityKepler
    relative_change_g = (intensity_g - intensityKepler)/intensityKepler
    relative_change_r = (intensity_r - intensityKepler)/intensityKepler
    relative_change_i = (intensity_i - intensityKepler)/intensityKepler
    relative_change_z = (intensity_z - intensityKepler)/intensityKepler
    relative_change_y = (intensity_y - intensityKepler)/intensityKepler

    u_band = [(1+relative_change_u) * flux for flux in flare_lc.flux]
    g_band = [(1+relative_change_g) * flux for flux in flare_lc.flux]
    r_band = [(1+relative_change_r) * flux for flux in flare_lc.flux]
    i_band = [(1+relative_change_i) * flux for flux in flare_lc.flux]
    z_band = [(1+relative_change_z) * flux for flux in flare_lc.flux]
    y_band = [(1+relative_change_y) * flux for flux in flare_lc.flux]

    fig1, (ax1, ax2) = plt.subplots(2,1,figsize=(15, 12))
    #fig1.suptitle('KIC {id} at {temp} K: Flare from {start} to {end}'.format(start = str(flare['St-BKJD']), end = str(flare['End-BKJD']), id = str(KIC_ID), temp = str(temp)))

    ax1.plot(flare_lc.time, flare_lc.flux * luminosity, color = 'black', label = 'Kepler K band', linestyle='-', marker='.')
    ax1.plot(flare_lc.time, u_band * luminosity, color = 'purple', label = 'Lsst u band', linestyle='-', marker='.')
    ax1.plot(flare_lc.time, g_band * luminosity, color = 'green', label = 'Lsst g band', linestyle='-', marker='.')
    ax1.plot(flare_lc.time, r_band * luminosity, color = 'orange', label = 'Lsst r band', linestyle='-', marker='.')
    ax1.plot(flare_lc.time, i_band * luminosity, color = 'red', label = 'Lsst i band', linestyle='-', marker='.')
    ax1.plot(flare_lc.time, z_band * luminosity, color = 'grey', label = 'Lsst z band', linestyle='-', marker='.')
    ax1.plot(flare_lc.time, y_band * luminosity, color = 'brown', label = 'Lsst y band', linestyle='-', marker='.')
    
    ax1.set_ylabel("Luminosity in Watts")
    ax1.set_xlabel("Time in BKJD days")

    ax1.legend()

    ax2.plot(flare_lc.time, flare_lc.flux * luminosity,color = 'black', label = 'Kepler K band', linestyle='-', marker='.')
    ax2.legend()
    ax2.set_ylabel("Luminosity in Watts")
    ax2.set_xlabel("Time in BKJD days")

    plt.show()