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
from distance import *
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

def get_baseline_luminosity_in_lsst_passband(flare_lc, KIC_ID, temp, luminosity):
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

    u_band = (1+relative_change_u) * luminosity.value
    g_band = (1+relative_change_g) * luminosity.value 
    r_band = (1+relative_change_r) * luminosity.value
    i_band = (1+relative_change_i) * luminosity.value
    z_band = (1+relative_change_z) * luminosity.value
    y_band = (1+relative_change_y) * luminosity.value
    

    dict = {
        'u': u_band,
        'g': g_band,
        'r': r_band,
        'i': i_band,
        'z': z_band,
        'y': y_band,
        'kep': luminosity.value,
    }
    return dict

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

def get_flare_luminosities_in_lsst_passbands(flare_lc, KIC_ID, temp, luminosity):
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

    kep_band = [flux * luminosity.value for flux in flare_lc.flux]
    u_band = [(1+relative_change_u) * flux * luminosity.value for flux in flare_lc.flux]
    g_band = [(1+relative_change_g) * flux * luminosity.value for flux in flare_lc.flux]
    r_band = [(1+relative_change_r) * flux * luminosity.value for flux in flare_lc.flux]
    i_band = [(1+relative_change_i) * flux * luminosity.value for flux in flare_lc.flux]
    z_band = [(1+relative_change_z) * flux * luminosity.value for flux in flare_lc.flux]
    y_band = [(1+relative_change_y) * flux * luminosity.value for flux in flare_lc.flux]

    fig1, (ax1, ax2) = plt.subplots(2,1,figsize=(15, 12))
    #fig1.suptitle('KIC {id} at {temp} K: Flare from {start} to {end}'.format(start = str(flare['St-BKJD']), end = str(flare['End-BKJD']), id = str(KIC_ID), temp = str(temp)))

    ax1.plot(flare_lc.time, kep_band, color = 'black', label = 'Kepler K band', linestyle='-', marker='.')
    ax1.plot(flare_lc.time, u_band, color = 'purple', label = 'Lsst u band', linestyle='-', marker='.')
    ax1.plot(flare_lc.time, g_band, color = 'green', label = 'Lsst g band', linestyle='-', marker='.')
    ax1.plot(flare_lc.time, r_band, color = 'orange', label = 'Lsst r band', linestyle='-', marker='.')
    ax1.plot(flare_lc.time, i_band, color = 'red', label = 'Lsst i band', linestyle='-', marker='.')
    ax1.plot(flare_lc.time, z_band, color = 'grey', label = 'Lsst z band', linestyle='-', marker='.')
    ax1.plot(flare_lc.time, y_band, color = 'brown', label = 'Lsst y band', linestyle='-', marker='.')
    
    ax1.set_ylabel("Luminosity in Watts / Hz")
    ax1.set_xlabel("Time in BKJD days")

    ax1.legend()

    ax2.plot(flare_lc.time, kep_band,color = 'black', label = 'Kepler K band', linestyle='-', marker='.')
    ax2.legend()
    ax2.set_ylabel("Luminosity in Watts / Hz")
    ax2.set_xlabel("Time in BKJD days")

    plt.show()

    dict = {
        'u': lk.LightCurve(time = flare_lc.time, flux = u_band, flux_err = flare_lc.flux_err),
        'g': lk.LightCurve(time = flare_lc.time, flux = g_band, flux_err = flare_lc.flux_err),
        'r': lk.LightCurve(time = flare_lc.time, flux = r_band, flux_err = flare_lc.flux_err),
        'i': lk.LightCurve(time = flare_lc.time, flux = i_band, flux_err = flare_lc.flux_err),
        'z': lk.LightCurve(time = flare_lc.time, flux = z_band, flux_err = flare_lc.flux_err),
        'y': lk.LightCurve(time = flare_lc.time, flux = y_band, flux_err = flare_lc.flux_err),
        'kep': lk.LightCurve(time = flare_lc.time, flux = kep_band, flux_err = flare_lc.flux_err),
    }
    return dict

def fit_flare_on_base(flare, base):

    # Reading flare data file.

    u_band = [(luminosity + base['u']) for luminosity in flare['u'].flux]
    g_band = [(luminosity + base['g']) for luminosity in flare['g'].flux]
    r_band = [(luminosity + base['r']) for luminosity in flare['r'].flux]
    i_band = [(luminosity + base['i']) for luminosity in flare['i'].flux]
    z_band = [(luminosity + base['z']) for luminosity in flare['z'].flux]
    y_band = [(luminosity + base['y']) for luminosity in flare['y'].flux]
    kep_band = y_band = [(luminosity + base['kep']) for luminosity in flare['kep'].flux]

    fig1, (ax1, ax2) = plt.subplots(2,1,figsize=(15, 12))
    #fig1.suptitle('KIC {id} at {temp} K: Flare from {start} to {end}'.format(start = str(flare['St-BKJD']), end = str(flare['End-BKJD']), id = str(KIC_ID), temp = str(temp)))

    ax1.plot(flare['kep'].time, kep_band, color = 'black', label = 'Kepler K band', linestyle='-', marker='.')
    ax1.plot(flare['u'].time, u_band, color = 'purple', label = 'Lsst u band', linestyle='-', marker='.')
    ax1.plot(flare['g'].time, g_band, color = 'green', label = 'Lsst g band', linestyle='-', marker='.')
    ax1.plot(flare['r'].time, r_band, color = 'orange', label = 'Lsst r band', linestyle='-', marker='.')
    ax1.plot(flare['i'].time, i_band, color = 'red', label = 'Lsst i band', linestyle='-', marker='.')
    ax1.plot(flare['z'].time, z_band, color = 'grey', label = 'Lsst z band', linestyle='-', marker='.')
    ax1.plot(flare['y'].time, y_band, color = 'brown', label = 'Lsst y band', linestyle='-', marker='.')
    
    ax1.set_ylabel("Luminosity in Watts / Hz")
    ax1.set_xlabel("Time in BKJD days")

    ax1.legend()


    ax2.plot(flare['kep'].time, kep_band,color = 'black', label = 'Kepler K band', linestyle='-', marker='.')
    ax2.legend()
    ax2.set_ylabel("Luminosity in Watts / Hz")
    ax2.set_xlabel("Time in BKJD days")

    plt.show()

    dict = {
        'u': lk.LightCurve(time = flare['u'].time, flux = u_band),
        'g': lk.LightCurve(time = flare['g'].time, flux = g_band),
        'r': lk.LightCurve(time = flare['r'].time, flux = r_band),
        'i': lk.LightCurve(time = flare['i'].time, flux = i_band),
        'z': lk.LightCurve(time = flare['z'].time, flux = z_band),
        'y': lk.LightCurve(time = flare['y'].time, flux = y_band),
        'kep': lk.LightCurve(time = flare['kep'].time, flux = kep_band),
    }
    return dict
    
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
