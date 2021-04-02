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
#from spectra_tools import get_kepler_transmission, compute_band_intensity

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

def get_fluxes_in_lsst_passbands(lsst_luminosities, distance):

    u_band = [(luminosity / (4 * math.pi * (distance ** 2))) for luminosity in lsst_luminosities['u'].flux]
    g_band = [(luminosity / (4 * math.pi * (distance ** 2)))for luminosity in lsst_luminosities['g'].flux]
    r_band = [(luminosity / (4 * math.pi * (distance ** 2))) for luminosity in lsst_luminosities['r'].flux]
    i_band = [(luminosity / (4 * math.pi * (distance ** 2))) for luminosity in lsst_luminosities['i'].flux]
    z_band = [(luminosity / (4 * math.pi * (distance ** 2)))for luminosity in lsst_luminosities['z'].flux]
    y_band = [(luminosity / (4 * math.pi * (distance ** 2)))for luminosity in lsst_luminosities['y'].flux]
    kep_band = y_band = [(luminosity / (4 * math.pi * (distance ** 2))) for luminosity in lsst_luminosities['kep'].flux]

    dict = {
        'u': lk.LightCurve(time = lsst_luminosities['u'].time, flux = u_band),
        'g': lk.LightCurve(time = lsst_luminosities['g'].time, flux = g_band),
        'r': lk.LightCurve(time = lsst_luminosities['r'].time, flux = r_band),
        'i': lk.LightCurve(time = lsst_luminosities['i'].time, flux = i_band),
        'z': lk.LightCurve(time = lsst_luminosities['z'].time, flux = z_band),
        'y': lk.LightCurve(time = lsst_luminosities['y'].time, flux = y_band),
        'kep': lk.LightCurve(time = lsst_luminosities['kep'].time, flux = kep_band),
    }
    
    return dict
    