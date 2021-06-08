import astropy.coordinates as coord
import numpy as np
import pandas as pd
import math

from lc_tools import  load_light_curve, get_flare_lc_from_time, get_normalized_lc, dump_modeled_data_to_LCLIB, add_LCLIB_header
from spectra_tools import get_baseline_luminosity_in_lsst_passband, get_flare_luminosities_in_lsst_passbands, fit_flare_on_base
from distance import get_luminosity_with_magnitude, get_mags_in_lsst_passbands
from ch_vars.spatial_distr import MilkyWayDensityJuric2008 as MWDensity
from plotting_tools import plotGenricSkyMapWithDistances, plotGenricSkyMap, plotGeneric2DHistogram, plotGenericHistogram

FLARE_DATA_PATH = 'data_files/filtered_flares.csv'

TOTAL_KEPLER_M_DWARF_COUNT = 4664
TOTAL_KEPLER_DURATION_IN_DAYS = 1460
TOTAL_FLARE_INSTANCES_COUNT = 103187

LOCAL_RADIUS = 20
NUMBER_OF_LOCAL_M_DWARFS = 1082

KEPLER_MEAN_EFFECTIVE_TEMP_FOR_M_DWARFS = 3743.4117647058824
KEPLER_STD_EFFECTIVE_TEMP_FOR_M_DWARFS = 161.37182827551771

def run_generator(flare_count):

    rng = np.random.default_rng(40)

    print("Sampling coordinates of stars")
    coordinates, distances = get_realistically_distributed_spherical_coordinates(flare_count, rng)

    print("Obtaining reference flare")
    kic_id, start_time, end_time = get_random_flare_events(flare_count, rng)
    
    print("Sampling star temperature")
    star_temp = get_normally_distributed_star_temp(flare_count, rng)
    
    print("Sampling flare temperature")
    flare_temp = get_normally_distributed_flare_temp(flare_count, rng)
    
    print("Begining flare modelling")
    add_LCLIB_header(flare_count)

    for i in range(flare_count):
        generate_model_flare_file(i, coordinates[i], distances[i], kic_id[i], start_time[i], end_time[i], star_temp[i], flare_temp[i])
        

def generate_model_flare_file(index, coordinates, distance, KIC_ID, start_time, end_time, star_temp, flare_temp):
    # Data extraction
    lc = load_light_curve(KIC_ID)
    flare_lc = get_flare_lc_from_time(lc, start_time, end_time)

    new_lc = get_normalized_lc(flare_lc)
    luminosity = get_luminosity_with_magnitude(KIC_ID).si  
    # Data modelling
    flare_luminosities = get_flare_luminosities_in_lsst_passbands(new_lc, KIC_ID, flare_temp, luminosity)
    
    baseline_luminosities = get_baseline_luminosity_in_lsst_passband(new_lc, KIC_ID, star_temp, luminosity)
    model_luminosities = fit_flare_on_base(flare_luminosities, baseline_luminosities)
    model_mags = get_mags_in_lsst_passbands(model_luminosities, distance)

    dump_modeled_data_to_LCLIB(index, coordinates.ra, coordinates.dec, KIC_ID, start_time, end_time, star_temp, flare_temp, distance, model_mags)
   
def get_number_of_expected_flares(radius, duration):
    """
    Computes the number of expected flares for a sphere of given radius during a given time period.
    This is based on data from the Kepler Space telescope and SUPERBLINK survey.

    Args:
        radius (float): Radius of sphere in parsec
        duration (float): Duration in days

    Returns:
        int : An estimate of the number of flares expcted during this time period in the sphere
    """

    flare_density = (NUMBER_OF_LOCAL_M_DWARFS / (4/3 * math.pi * LOCAL_RADIUS**2)) * (TOTAL_FLARE_INSTANCES_COUNT / (TOTAL_KEPLER_M_DWARF_COUNT * TOTAL_KEPLER_DURATION_IN_DAYS))
    volume = 4/3 * math.pi * (radius ** 3)
    flare_count = math.floor(flare_density * volume * duration)

    return flare_count

def get_realistically_distributed_spherical_coordinates(count, rng):
    mw = MWDensity()
    coordinates = mw.sample_eq(count, rng)
    return coordinates, coordinates.distance.value


def get_uniformly_distributed_spherical_coordinates(radius, count, rng, chunk_size=1<<10):
    """
    Returns a collection of uniformally distributed coordinates (RA, Dec) and distances that all fall 
    within a sphere of the parameter radius. 

    Args:
        radius (float): The radius of the sphere in parsec.
        count (int): Number of coordinates to be returned.
        chunk_size ([type], optional): [description]. Defaults to 1<<10.

    Returns:
        tuple: Contains two numpy arrays containing the skycoords and distances in parsec.
    """

    x_ = []
    y_ = []
    z_ = []
    n = 0
    while n < count:
        x, y, z = rng.uniform(-1 * radius,radius, (3, chunk_size))
        idx = x ** 2 + y ** 2 + z ** 2 < radius * radius
        x_.append(x[idx])
        y_.append(y[idx])
        z_.append(z[idx])
        n += np.count_nonzero(idx)
    x = np.concatenate(x_)
    y = np.concatenate(y_)
    z = np.concatenate(z_)
    
    r, dec, ra = coord.cartesian_to_spherical(x[:count], y[:count], z[:count])
    coordinates = coord.SkyCoord(ra=ra, dec=dec)
    return coordinates, r.value



def get_random_flare_events(count, rng):
    """
    Returns a tuple of 3 numpy arrays containing the Kepler Input Catalogue ID, flare start time and flare end time 
    of flares randomly selected from the filtered_flares.csv file.

    Args:
        count (int): Number of flares for which the data needs to be returned.

    Returns:
        tuple: Contains numpy arrays for the KIC ID, flare start time and flare end time respectively.
    """

    KIC = []
    St_time = []
    End_time = []

    df = pd.read_csv(FLARE_DATA_PATH)
    indices = rng.integers(low = 0, high = len(df), size = count)
    
    for index in indices:
        KIC.append(df['KIC'][index])
        St_time.append(df['St-BKJD'][index])
        End_time.append(df['End-BKJD'][index])
    
    return KIC, St_time, End_time


def get_normally_distributed_star_temp(count, rng):
    """
    Returns a numpy array of star temperatures modelled after a normal distribution of effective star temperatures
    based on data from the Kepler Input Catalogue.

    Args:
        count (int): Length of numpy array to be returned

    Returns:
        numpy array: numpy array containing the star temperatures with length = count
    """

    return rng.normal(KEPLER_MEAN_EFFECTIVE_TEMP_FOR_M_DWARFS, KEPLER_STD_EFFECTIVE_TEMP_FOR_M_DWARFS, count)


def get_normally_distributed_flare_temp(count, rng):
    """
    Returns a numpy array of flare temperatures modelled after a normal distribution.

    Args:
        count (int): Length of numpy array to be returned

    Returns:
        numpy array: numpy array containing the flare temperatures with length = count
    """

    return rng.normal(30000, 500, count)

if __name__ == "__main__":
    run_generator(1000)