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
from scipy.integrate import simps
from astropy.modeling.blackbody import blackbody_lambda as bb

passband_props = dict(
    u=dict(eff=3751.36, min=3205.54, max=4081.24, width=473.19),
    g=dict(eff=4741.64, min=3873.01, max=5665.07, width=1253.26),
    r=dict(eff=6173.23, min=5375.74, max=7054.95, width=1206.92),
    i=dict(eff=7501.62, min=6765.00, max=8324.70, width=1174.76),
    z=dict(eff=8679.19, min=8034.98, max=9374.71, width=997.51),
    y=dict(eff=9711.53, min=9088.94, max=10897.23, width=871.76),
)

bands = tuple(passband_props)

def get_spectra_data(lc, KIC_ID, temp):
    # Reading flare data file.
    flare_data = astropy.io.ascii.read(FLARE_DATA_PATH, quotechar="\s")

    #Plotting points identified as flares.
    for flare in flare_data:
        if flare['KIC'] != int(KIC_ID):
            continue
            
        # Obtaining start and end indices for flare instances
        flare_lc = get_flare_lc(lc,flare)

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
        fig1.suptitle('KIC {id} at {temp} K: Flare from {start} to {end}'.format(start = str(flare['St-BKJD']), end = str(flare['End-BKJD']), id = str(KIC_ID), temp = str(temp)))

        ax1.plot(flare_lc.time, flare_lc.flux, color = 'black', label = 'Kepler K band', linestyle='-', marker='.')
        ax1.plot(flare_lc.time, u_band, color = 'purple', label = 'Lsst u band', linestyle='-', marker='.')
        ax1.plot(flare_lc.time, g_band, color = 'green', label = 'Lsst g band', linestyle='-', marker='.')
        ax1.plot(flare_lc.time, r_band, color = 'orange', label = 'Lsst r band', linestyle='-', marker='.')
        ax1.plot(flare_lc.time, i_band, color = 'red', label = 'Lsst i band', linestyle='-', marker='.')
        ax1.plot(flare_lc.time, z_band, color = 'grey', label = 'Lsst z band', linestyle='-', marker='.')
        ax1.plot(flare_lc.time, y_band, color = 'brown', label = 'Lsst y band', linestyle='-', marker='.')
        
        ax1.set_ylabel("Relative Flux")
        ax1.set_xlabel("Time")

        ax1.legend()

        ax2.plot(flare_lc.time, flare_lc.flux,color = 'black', label = 'Kepler K band', linestyle='-', marker='.')
        ax2.legend()
        ax2.set_ylabel("Relative Flux")
        ax2.set_xlabel("Time")

        plt.show()

def get_spectra_data(flare_lc, KIC_ID, temp):
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

    ax1.plot(flare_lc.time, flare_lc.flux, color = 'black', label = 'Kepler K band', linestyle='-', marker='.')
    ax1.plot(flare_lc.time, u_band, color = 'purple', label = 'Lsst u band', linestyle='-', marker='.')
    ax1.plot(flare_lc.time, g_band, color = 'green', label = 'Lsst g band', linestyle='-', marker='.')
    ax1.plot(flare_lc.time, r_band, color = 'orange', label = 'Lsst r band', linestyle='-', marker='.')
    ax1.plot(flare_lc.time, i_band, color = 'red', label = 'Lsst i band', linestyle='-', marker='.')
    ax1.plot(flare_lc.time, z_band, color = 'grey', label = 'Lsst z band', linestyle='-', marker='.')
    ax1.plot(flare_lc.time, y_band, color = 'brown', label = 'Lsst y band', linestyle='-', marker='.')
    
    ax1.set_ylabel("Relative Flux")
    ax1.set_xlabel("Time in BKJD days")

    ax1.legend()

    ax2.plot(flare_lc.time, flare_lc.flux,color = 'black', label = 'Kepler K band', linestyle='-', marker='.')
    ax2.legend()
    ax2.set_ylabel("Relative Flux")
    ax2.set_xlabel("Time in BKJD days")

    plt.show()

def compute_band_intensity(band, temp):
    T = temp * u.K
    lmbd, t = get_transmission(band)
    intensity = simps(x=lmbd, y=bb(lmbd * u.AA, T).value * t) / simps(x=lmbd, y=t)
    return intensity

def get_transmission(band):
    lmbd, t = np.genfromtxt(f'filters/LSST_LSST.{band}.dat', unpack=True)
    return lmbd, t

def get_kepler_transmission():
    lmbd, t = np.genfromtxt('filters/Kepler_Kepler.K.dat', unpack=True)
    return lmbd, t
