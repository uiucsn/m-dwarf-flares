import astropy
from astropy import units as u
import lightkurve as lk
import math
import numpy as np

dist_data = astropy.io.ascii.read('data_files/dist_new.csv')
mag_data = astropy.io.ascii.read('data_files/mag.csv')

def get_luminosity_with_magnitude(KIC):
    """
    Computes the luminosity of a KIC object based on distance data from Gaia
    DR3 and Kepler band magnitude.

    Args:
        KIC (Kepler): Kepler ID for which luminosity needs to be found.

    Returns:
        luminosity: Estimated luminosity value based on Kepler and Gaia Data.
    """

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

def get_mags_in_lsst_passbands(model_luminosities, distance):
    """
    Converts a dictionary of luminosity time sequences in lsst and kepler passbands to their
    magnitude values based on the distance.

    Args:
        model_luminosities (dictionary of light cruves): A dictionary of luminosity time sequences
        in all 6 lsst passbands and the Kepler passbands.

        distance (parsec): Distance values in parsec which is used to compute the magnitude.

    Returns:
        magnitudes: A dictionary of magnitude time sequences in lsst and kepler passbands.
    """

    u_band = [((lum * u.Watt * u.s)/ (4 * math.pi * (distance * u.pc)**2)).si.to(u.ABmag).value for lum in model_luminosities['u'].flux]
    g_band = [((lum * u.Watt * u.s) / (4 * math.pi * (distance * u.pc)**2)).si.to(u.ABmag).value for lum in model_luminosities['g'].flux]
    r_band = [((lum * u.Watt * u.s) / (4 * math.pi * (distance * u.pc)**2)).si.to(u.ABmag).value for lum in model_luminosities['r'].flux]
    i_band = [((lum * u.Watt * u.s) / (4 * math.pi * (distance * u.pc)**2)).si.to(u.ABmag).value for lum in model_luminosities['i'].flux]
    z_band = [((lum * u.Watt * u.s) / (4 * math.pi * (distance * u.pc)**2)).si.to(u.ABmag).value for lum in model_luminosities['z'].flux]
    y_band = [((lum * u.Watt * u.s) / (4 * math.pi * (distance * u.pc)**2)).si.to(u.ABmag).value for lum in model_luminosities['y'].flux]
    kep_band = [((lum * u.Watt * u.s) / (4 * math.pi * (distance * u.pc)**2)).si.to(u.ABmag).value for lum in model_luminosities['kep'].flux]

    dict = {
        'u': lk.LightCurve(time = model_luminosities['u'].time, flux = u_band),
        'g': lk.LightCurve(time = model_luminosities['g'].time, flux = g_band),
        'r': lk.LightCurve(time = model_luminosities['r'].time, flux = r_band),
        'i': lk.LightCurve(time = model_luminosities['i'].time, flux = i_band),
        'z': lk.LightCurve(time = model_luminosities['z'].time, flux = z_band),
        'y': lk.LightCurve(time = model_luminosities['y'].time, flux = y_band),
        'kep': lk.LightCurve(time = model_luminosities['kep'].time, flux = kep_band),
    }

    return dict
