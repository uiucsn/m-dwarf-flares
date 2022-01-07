import astropy
from astropy import units as u
import lightkurve as lk
import math
import numpy as np
import os

dist_data = astropy.io.ascii.read(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data_files','dist_new.csv'))
mag_data = astropy.io.ascii.read(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data_files','mag.csv'))

def get_stellar_luminosity(KIC):
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
            if row['r_med_photogeo'] != 0 and not np.isnan(row['r_med_photogeo']):
                distance = row['r_med_photogeo']
            else:
                distance = row['r_med_geo']
    
    for row in mag_data:
        if row['KIC ID'] == int(KIC):
            ab_magnitude = row['kic_kepmag'] 

    flux = (ab_magnitude * u.ABmag).to(u.Jy)
    luminosity = flux * 4 * math.pi * (distance * u.pc) ** 2
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

    dict = {
        band: lk.LightCurve(
            time=model_lum.time,
            flux=((model_lum.flux * u.Watt * u.s)/ (4 * math.pi * (distance)**2)).si.to(u.ABmag).value,
        )
        for band, model_lum in model_luminosities.items()
    }

    return dict
