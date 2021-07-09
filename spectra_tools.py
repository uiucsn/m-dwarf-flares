from functools import lru_cache

import lightkurve as lk
import numpy as np
from astropy import units as u
from astropy.modeling.blackbody import blackbody_lambda as bb
from scipy.integrate import simps

passband_props = dict(
    u=dict(eff=3751.36, min=3205.54, max=4081.24, width=473.19),
    g=dict(eff=4741.64, min=3873.01, max=5665.07, width=1253.26),
    r=dict(eff=6173.23, min=5375.74, max=7054.95, width=1206.92),
    i=dict(eff=7501.62, min=6765.00, max=8324.70, width=1174.76),
    z=dict(eff=8679.19, min=8034.98, max=9374.71, width=997.51),
    y=dict(eff=9711.53, min=9088.94, max=10897.23, width=871.76),
)

bands = tuple(passband_props)

def get_baseline_luminosity_in_lsst_passband(flare_lc, KIC_ID, spectrum, luminosity):
    """
    Returns the baseline luminosities for all six LSST passbands and the Kepler band based 
    on distance data from Gaia DR3, magnitude data in the Kepler band, and spectrum data
    modelled after a simple black body.

    Args:
        flare_lc (flare light curve): The light kurve object of the flare.
        KIC_ID (Kepler ID): The Kepler Input Catalogue building.
        temp (float): The temperature of the star in Kelvin
        luminosity (float): The baseline luminosity in the Kepler band.

    Returns:
        Dictionary : A dictionary with baseline luminosity values in all 6 LSST passbands
        and the Kepler passband.
    """
    lmbd, t = get_kepler_transmission()
    intensityKepler = simps(x=lmbd, y=spectrum(lmbd * u.AA).value * t) / simps(x=lmbd, y=t)

    intensity_u = compute_band_intensity("u", spectrum)
    intensity_g = compute_band_intensity("g", spectrum)
    intensity_r = compute_band_intensity("r", spectrum)
    intensity_i = compute_band_intensity("i", spectrum)
    intensity_z = compute_band_intensity("z", spectrum)
    intensity_y = compute_band_intensity("y", spectrum)

    kep_band = luminosity.value
    u_band = luminosity.value * (intensity_u / intensityKepler)
    g_band = luminosity.value * (intensity_g / intensityKepler)
    r_band = luminosity.value * (intensity_r / intensityKepler)
    i_band = luminosity.value * (intensity_i / intensityKepler)
    z_band = luminosity.value * (intensity_z / intensityKepler)
    y_band = luminosity.value * (intensity_y / intensityKepler)
    

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

def get_flare_luminosities_in_lsst_passbands(flare_lc, KIC_ID, spectrum, luminosity):
    """
    Returns the flare luminosities for all six LSST passbands and the Kepler band based 
    on distance data from Gaia DR3, magnitude data in the Kepler band, and spectrum data
    modelled after a simple black body.

    Args:
        flare_lc (flare light curve): The light kurve object of the flare.
        KIC_ID (Kepler ID): The Kepler Input Catalogue building.
        temp (float): The temperature of the star in Kelvin
        luminosity (float): The baseline luminosity in the Kepler band.

    Returns:
        Dictionary : A dictionary with flare luminosity values in all 6 LSST passbands
        and the Kepler passband.
    """

    lmbd, t = get_kepler_transmission()
    intensityKepler = simps(x=lmbd, y=spectrum(lmbd * u.AA).value * t) / simps(x=lmbd, y=t)

    intensity_u = compute_band_intensity("u", spectrum)
    intensity_g = compute_band_intensity("g", spectrum)
    intensity_r = compute_band_intensity("r", spectrum)
    intensity_i = compute_band_intensity("i", spectrum)
    intensity_z = compute_band_intensity("z", spectrum)
    intensity_y = compute_band_intensity("y", spectrum)

    kep_band = flare_lc.flux * luminosity.value
    u_band = flare_lc.flux * luminosity.value * (intensity_u / intensityKepler)
    g_band = flare_lc.flux * luminosity.value * (intensity_g / intensityKepler)
    r_band = flare_lc.flux * luminosity.value * (intensity_r / intensityKepler)
    i_band = flare_lc.flux * luminosity.value * (intensity_i / intensityKepler)
    z_band = flare_lc.flux * luminosity.value * (intensity_z / intensityKepler)
    y_band = flare_lc.flux * luminosity.value * (intensity_y / intensityKepler)

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
    """
    Retuns a dictionary of luminosities in all LSST passbands and the Kepler passband with the 
    the flare luminosities fitted on top of the baseline luminosities. 

    Args:
        flare (Dictionary of light kurves): Flare luminosities in all 6 lsst passbands and kepler passband
        base ([Dictionary of numbers]): Base luminosities in all 6 lsst passbands and the kepler passband

    Returns:
        Dictionary of luminosities: Light kurves objects in all lsst passbands and the kepler passband with 
        the flare fitted on top of the baseline luminosities.
    """
    u_band = flare['u'].flux + base['u']
    g_band = flare['g'].flux + base['g']
    r_band = flare['r'].flux + base['r']
    i_band = flare['i'].flux + base['i']
    z_band = flare['z'].flux + base['z']
    y_band = flare['y'].flux + base['y']
    kep_band = flare['kep'].flux + base['kep']

    # Adding baseline stellar luminosity observations before and after the flare for LCLIB
    u_band = np.concatenate([[base['u']],u_band,[base['u']]])
    g_band = np.concatenate([[base['g']],g_band,[base['g']]])
    r_band = np.concatenate([[base['r']],r_band,[base['r']]])
    i_band = np.concatenate([[base['i']],i_band,[base['i']]])
    z_band = np.concatenate([[base['z']],z_band,[base['z']]])
    y_band = np.concatenate([[base['y']],y_band,[base['y']]])
    kep_band =  np.concatenate([[base['kep']],kep_band,[base['kep']]])

    # Adding times for the new luminosity points
    start_time = np.amin(flare['kep'].time)
    end_time = np.amax(flare['kep'].time)
    delta_t = (end_time - start_time) / len(flare['kep'].time)
    time =  np.concatenate([[start_time - delta_t],flare['kep'].time,[end_time + delta_t]])

    dict = {
        'u': lk.LightCurve(time = time, flux = u_band),
        'g': lk.LightCurve(time = time, flux = g_band),
        'r': lk.LightCurve(time = time, flux = r_band),
        'i': lk.LightCurve(time = time, flux = i_band),
        'z': lk.LightCurve(time = time, flux = z_band),
        'y': lk.LightCurve(time = time, flux = y_band),
        'kep': lk.LightCurve(time = time, flux = kep_band),
    }
    return dict
    
def compute_band_intensity(band, spectrum):
    """
    Computes the intensity in a particular lsst passband at a certain temperature, 
    modelled after a black body.

    Args:
        band (char): Letter corresponding to an lsst passband
        temp (temperature): Temperature of the black body at a given temperature.

    Returns:
        Intensity: Internsity in the given passband at the given temperature.
    """
    lmbd, t = get_transmission(band)
    intensity = simps(x=lmbd, y=spectrum(lmbd * u.AA).value * t) / simps(x=lmbd, y=t)
    return intensity

@lru_cache()
def get_transmission(band, n_grid=200):
    """
    Returns the transmission for a given lsst passband.

    Args:
        band (char): Letter corresponding to an lsst passband.
        n_grid (int): Number of grid points

    Returns:
        [type]: [description]
    """
    lmbd, t = np.genfromtxt(f'filters/LSST_LSST.{band}.dat', unpack=True)
    step = lmbd.size // n_grid or 1
    return lmbd[::step], t[::step]

@lru_cache()
def get_kepler_transmission(n_grid=200):
    """
    Returns the transmission for the kepler passband.

    Args:
        n_grid (int): Number of grid points

    Returns:
        [type]: [description]
    """
    lmbd, t = np.genfromtxt('filters/Kepler_Kepler.K.dat', unpack=True)
    step = lmbd.size // n_grid or 1
    return lmbd[::step], t[::step]
