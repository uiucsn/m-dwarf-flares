from functools import lru_cache

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from dustmaps.bayestar import BayestarQuery
from dustmaps.sfd import SFDQuery


BAYESTAR = BayestarQuery(version='bayestar2019')
SFD = SFDQuery()


def get_extinction_after_symmetric_interpolation(coordinates):
    """
    Interpolates the Bayestar Dust Data for the entire sky by inverting the galactic coordinates of 
    points where real data is not avalible.

    Args:
        coordinates (Astropy sky coords): coordinates for which extinction needs to be calculated.

    Returns:
        Bayestar 19 reddening values: The Reddeining values from the Bayestar 19 map.
    """

    gal = coordinates.galactic
    mirrored_coordinates = SkyCoord(l = - 1 * gal.l, b = -1 * gal.b, distance = gal.distance, frame = 'galactic').icrs

    reddening = BAYESTAR(coordinates, mode='median')
    reddening_mirrored = BAYESTAR(mirrored_coordinates, mode='median')

    reddening_interpolated = np.where(np.isnan(reddening), reddening_mirrored, reddening)

    return reddening_interpolated

def get_extinction_after_symmetric_interpolation_with_sfd_factor(coordinates):
    """
    Interpolates the Bayestar Dust Data for the entire sky by inverting the galactic coordinates of 
    points where real data is not avalible and applying a correction factor from SFD 2D data.

    Args:
        coordinates (Astropy sky coords): coordinates for which extinction needs to be calculated.

    Returns:
        Bayestar 19 reddening values: The Reddeining values from the Bayestar 19 map.
    """

    gal = coordinates.galactic
    mirrored_coordinates = SkyCoord(l = - 1 * gal.l, b = -1 * gal.b, distance = gal.distance, frame = 'galactic').icrs

    reddening_bayes = BAYESTAR(coordinates, mode='median')
    reddening_bayes_mirrored = BAYESTAR(mirrored_coordinates, mode='median')
    reddening_sfd_norm = SFD(coordinates)
    reddening_sfd_flip = SFD(mirrored_coordinates)

    # This replaces any undefined proportionality factors (i.e 1 / 0) by 1
    sfd_factor = np.divide(reddening_sfd_norm, reddening_sfd_flip, out = np.ones_like(reddening_sfd_norm), where = reddening_sfd_flip != 0)
    interpolated_reddening_3D = sfd_factor * reddening_bayes_mirrored
    reddening_interpolated = np.where(np.isnan(reddening_bayes), interpolated_reddening_3D, reddening_bayes)

    return reddening_interpolated

@lru_cache()
def get_LSST_extinction():
    df = pd.read_csv('data_files/LSST_rv_3_1.csv')
    rv = pd.Series(df['value'].values,index=df['lsst_passband']).to_dict()
    return rv

def get_extinction_in_lsst_passbands(coordinates):
    """
    Computes extinction in lsst passbands after applying appropriate correction factors to 
    results from the get_extinction_after_symmetric_interpolation_with_sfd_factor() function

    Args:
        coordinates (Astropy sky coords): coordinates for which extinction needs to be calculated in lsst passbands.

    Returns:
        dictionary: Extinction values in the lsst ugrizy passbands for the given coordinates.
    """

    rv = get_LSST_extinction()

    bayestar_coefficient = 0.884
    bayestar_reddening = get_extinction_after_symmetric_interpolation_with_sfd_factor(coordinates)
    e_bv = bayestar_coefficient * bayestar_reddening

    lsst_ext = {
        'u': rv['u'] * e_bv,
        'g': rv['g'] * e_bv,
        'r': rv['r'] * e_bv,
        'i': rv['i'] * e_bv,
        'z': rv['z'] * e_bv,
        'y': rv['y'] * e_bv,
    }
    return lsst_ext

def apply_extinction_to_lsst_mags(magnitudes, extinction_values):
    """
    Applies the extinction value to the lsst model magnitudes

    Args:
        magnitudes (dictionary): A dictionary of lightcurves in lsst ugrizy and kep passbands
        extinction_values (dictionary): A dictionary of extinction values in lsst ugrizy passbands

    Returns:
        dictionary: A dictionary of lightcurves in lsst ugrizy and kep passbands
    """
    lsst_mags_with_ext = {
        'u': magnitudes['u'] + extinction_values['u'],
        'g': magnitudes['g'] + extinction_values['g'],
        'r': magnitudes['r'] + extinction_values['r'],
        'i': magnitudes['i'] + extinction_values['i'],
        'z': magnitudes['z'] + extinction_values['z'],
        'y': magnitudes['y'] + extinction_values['y'],
        'kep': magnitudes['kep'],
    }    
    return lsst_mags_with_ext

def plot_extinction_after_symmetric_interpolation():
    """
    Plots a sky map of extinction computed using get_extinction_after_symmetric_interpolation()
    """

    ra, dec = np.meshgrid(np.linspace(-180, 180, 361, endpoint=False), np.linspace(-90, 90, 181, endpoint=False))
    coordinates =  SkyCoord(ra = ra * u.degree, dec = dec * u.degree, distance = 2000 * u.pc)
    reddening = get_extinction_after_symmetric_interpolation(coordinates)

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection="mollweide")
    scatter = ax.scatter(-coordinates.ra.wrap_at(180 * u.deg).radian, coordinates.dec.wrap_at(180 * u.deg).radian,c = reddening, cmap = 'gist_stern', s=3, vmin=0)
    ax.grid(True)
    ax.set_xticklabels(['10h', '8h', '6h', '4h', '2h', '0h', '22h', '20h', '18h', '16h', '14h'])
    plt.colorbar(scatter, label = "Extinction")
    plt.show()

def plot_extinction_after_symmetric_interpolation_with_sfd_factor():
    """
    Plots a sky map of extinction computed using get_extinction_after_symmetric_interpolation_with_sfd_factor()
    """

    ra, dec = np.meshgrid(np.linspace(-180, 180, 361, endpoint=False), np.linspace(-90, 90, 181, endpoint=False))
    coordinates =  SkyCoord(ra = ra * u.degree, dec = dec * u.degree, distance = 2000 * u.pc)
    reddening = get_extinction_after_symmetric_interpolation_with_sfd_factor(coordinates)

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection="mollweide")
    scatter = ax.scatter(-coordinates.ra.wrap_at(180 * u.deg).radian, coordinates.dec.wrap_at(180 * u.deg).radian,c = reddening, cmap = 'gist_stern', s=3, vmin=0)
    ax.grid(True)
    ax.set_xticklabels(['10h', '8h', '6h', '4h', '2h', '0h', '22h', '20h', '18h', '16h', '14h'])
    plt.colorbar(scatter, label = "Extinction")
    plt.show()
