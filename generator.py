from functools import lru_cache

import astropy.coordinates as coord
import numpy as np
import pandas as pd
import math
import time
import argparse
import sys
import progressbar

from lc_tools import  load_light_curve, get_flare_lc_from_time, get_normalized_lc, dump_modeled_data_to_LCLIB, add_LCLIB_header
from spectra_tools import get_baseline_luminosity_in_lsst_passband, get_flare_luminosities_in_lsst_passbands, fit_flare_on_base
from distance import get_stellar_luminosity, get_mags_in_lsst_passbands
from ch_vars.spatial_distr import MilkyWayDensityJuric2008 as MWDensity
from plotting_tools import save_simulation_plots
from extinction_tools import get_extinction_in_lsst_passbands, apply_extinction_to_lsst_mags

FLARE_DATA_PATH = 'data_files/filtered_flares.csv'

# Do not change these values
KEPLER_MEAN_EFFECTIVE_TEMP_FOR_M_DWARFS = 3743.4117647058824
KEPLER_STD_EFFECTIVE_TEMP_FOR_M_DWARFS = 161.37182827551771

# Minimum Relative Flux Amplitude of the flares. Flares below this relative flux amplitude will be filtered out.
MIN_RELATIVE_FLUX_AMPLITUDE = 0.01
# Maximum magnitude for a flare in the LSST u passband. Flares above this mag value will be filtered out.
PEAK_MAGNITUDE_THRESHOLD = {
    'u': 23.66 + 1,
    'g': 24.69 + 1,
    'r': 24.06 + 1,
    'i': 23.45 + 1,
    'z': 22.54 + 1,
    'y': 21.62 + 1,
}
# Minimum magnitude amplitude of the simulated flare in the u passband. Flares below this mag amplitude will be filtered out.
U_BAND_AMPLITUDE_THRESHOLD = 0.1 

PARAMETER_COUNT_MULTIPLIER = 50

def run_generator(flare_count, file_path, start_index, remove_header, to_plot):
    """
    Runs the generator functions. Samples the respective distributions for the parameters and writes
    simulated flare instances to an LCLIB file. 

    Args:
        flare_count (int): number of flares to be generated
        file_path (int): Path of the file where the results need to be saved.
    """

    number_of_nominal_flares = 0
    number_of_simulated_flares = 0
    nominal_flare_indices = []
    nominal_flare_instance = []
    rng = np.random.default_rng(start_index)
    parameter_count = PARAMETER_COUNT_MULTIPLIER * flare_count # Generating more parameters to avoid reloading of dust map and other files

    with open(file_path, 'w') as output_file:

        # Adding lc lib to header
        if not remove_header:
            add_LCLIB_header(flare_count, output_file)

        with progressbar.ProgressBar(max_value = flare_count) as bar:
            # While loop keeps executing until the number of flares generated matches the flare count
            while (number_of_nominal_flares < flare_count):

                print("--- Generating new parameters for flare modelling ---")

                print("1. Sampling coordinates of stars ...")
                coordinates = get_realistically_distributed_spherical_coordinates(parameter_count, rng)
                galactic_coordinates = coord.SkyCoord(coordinates).galactic
                distances = coordinates.distance
                
                print("2. Computing extinction values in lsst passbands ...")
                extinction_values = get_extinction_in_lsst_passbands(coord.SkyCoord(coordinates))

                print("3. Obtaining reference flares ...")
                kic_id, start_time, end_time = get_random_flare_events(parameter_count, rng, MIN_RELATIVE_FLUX_AMPLITUDE)
                
                print("4. Sampling star temperature ...")
                star_temp = get_normally_distributed_star_temp(parameter_count, rng)
                
                print("5. Sampling flare temperature ...")
                flare_temp = get_normally_distributed_flare_temp(parameter_count, rng)
                
                print("6. Commencing flare modelling ...")
                for i in range(parameter_count):
                    # Breaking out of the loop if the correct number of flares are generated
                    bar.update(number_of_nominal_flares)
                    if (number_of_nominal_flares == flare_count):
                        number_of_simulated_flares = i
                        break
                    extinction = {
                        'u': extinction_values['u'][i],
                        'g': extinction_values['g'][i],
                        'r': extinction_values['r'][i],
                        'i': extinction_values['i'][i],
                        'z': extinction_values['z'][i],
                        'y': extinction_values['y'][i],
                    }
                    is_valid_flare, modeled_flare = generate_model_flare_file(start_index + number_of_nominal_flares, coordinates[i], galactic_coordinates[i], distances[i], kic_id[i], start_time[i], end_time[i], star_temp[i], flare_temp[i], extinction, output_file)
                    if is_valid_flare:
                        number_of_nominal_flares += 1
                        if to_plot:
                            nominal_flare_indices.append(i)
                            nominal_flare_instance.append(modeled_flare)
    output_file.close()
    print("{:2f} %% of the simulations passed the LSST observation thresholds".format((number_of_nominal_flares / number_of_simulated_flares) * 100))
    if to_plot:
        nominal_coordinates = coord.SkyCoord(coordinates[nominal_flare_indices])
        print("7. Generating plots ...")
        save_simulation_plots(nominal_coordinates, nominal_flare_instance, rng)


def generate_model_flare_file(index, coordinates, galactic_coordinates, distance, KIC_ID, start_time, end_time, star_temp, flare_temp, extinction, output_file):
    """
    Generates the model flare based on the parameters and saves to the LCLIB file if it makes the threshold cuts.

    Args:
        index (int): Index of the flare instance.
        coordinates (Sky coord): Coordinates of the flare instance.
        galactic_coordinates (Sky coord): Galactic Coordinates of the flare instance.
        distance (astropy distance unit): Distance of the flare instance.
        KIC_ID (int): Kepler Input Catalogue ID
        start_time (float): Start time of the flare instance
        end_time (float): End time of the flare instance
        star_temp (float): Star temperature for spectrum modelling.
        flare_temp (float): Flare temperature for flare modelling.
        extinction (dictionary of floats): Extinction values in lsst passbands.
        output_file (string): path of the output file.

    Returns:
        tuple: boolean, modelled flare
    """

    # Loading the orignal flare light curve and normalizing it
    lc = load_light_curve(KIC_ID)
    flare_lc = get_flare_lc_from_time(lc, start_time, end_time)
    new_lc = get_normalized_lc(flare_lc)

    # Computing the stellar luminosity based on Kepler and Gaia Data
    luminosity = get_stellar_luminosity(KIC_ID).si  

    # Modelling the spetra and fitting the flare on the nominal stellar luminosity
    flare_luminosities = get_flare_luminosities_in_lsst_passbands(new_lc, KIC_ID, flare_temp, luminosity)
    baseline_luminosities = get_baseline_luminosity_in_lsst_passband(new_lc, KIC_ID, star_temp, luminosity)
    model_luminosities = fit_flare_on_base(flare_luminosities, baseline_luminosities)

    # Converting luminosity to distances based on distances and applying extinction
    model_mags = get_mags_in_lsst_passbands(model_luminosities, distance)
    model_mags_with_extinction = apply_extinction_to_lsst_mags(model_mags, extinction)

    if is_nominal_flare(model_mags_with_extinction):
        # Writing modelled data to LCLIB file if the flare is nominal
        dump_modeled_data_to_LCLIB(index, galactic_coordinates.l, galactic_coordinates.b, KIC_ID, start_time, end_time, star_temp, flare_temp, distance, model_mags_with_extinction, output_file)
        return True, model_mags_with_extinction
    else:
        return False, model_mags_with_extinction

def is_nominal_flare(flare):
    """
    Checking if the generated flare is under the threshold for max magnitude and has an
    amplitude greater than the threshold

    Args:
        flare (dicitonary): A dictionary of lightcurves in the LSST ugrizy and kep passbands

    Returns:
        [boolean]: True if the flare is nominal, False otherwise.
    """

    passbands = ['u','g','r','i','z','y']
    if all(np.all(flare[passband].flux > PEAK_MAGNITUDE_THRESHOLD[passband]) for passband in passbands):
        # Checking if all bands have magnitude greater than the correspnding PEAK_MAGNITUDE_THRESHOLD values. 
        # Sim is nominal if even one passband has 
        return False
    if np.ptp(flare['u'].flux) <= U_BAND_AMPLITUDE_THRESHOLD:
        # Checking if u band has amplitude less than U_BAND_AMPLITUDE
        return False

    return True

def get_number_of_expected_flares():
    """
    Computes the number of expected flares that can be observed by LSST per day.

    Returns:
        int : An estimate of the number of flares per day that fall within the LSST magnitude limits.
    """
    TOTAL_KEPLER_M_DWARF_COUNT = 4664
    TOTAL_KEPLER_DURATION_IN_DAYS = 1460
    TOTAL_KEPLER_FLARE_COUNT = 103187

    # Estimated number of m dwarfs in solar neighbourhood
    N_M_DWARF = 24.8 * (10 ** 9)

    # Fraction of simulated flares that fall within the LSST thresholds
    LSST_VISIBILITY_FRACTION = .1170

    # Unit: Number of flares per m dwarf per day
    flares_per_star_per_day = TOTAL_KEPLER_FLARE_COUNT / (TOTAL_KEPLER_M_DWARF_COUNT * TOTAL_KEPLER_DURATION_IN_DAYS)

    # Unit: Total number of flares per day
    total_flares_per_day = flares_per_star_per_day * N_M_DWARF

    # Unit: Total number of flares per day potentially visible to LSST
    flare_count = total_flares_per_day * LSST_VISIBILITY_FRACTION

    return int(flare_count)

def get_realistically_distributed_spherical_coordinates(count, rng):
    """
    Uses the the milky way juric 2008 ch_vars module by Malanchev Kostya
    to get realistic distribution of coordinates in the Milky way.

    Args:
        count (int): Number of coordinates that need to be generated
        rng ([type]): [description]

    Returns:
        [type]: [description]
    """
    mw = MWDensity()
    coordinates = mw.sample_eq(count, rng)
    return coordinates


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
    return coordinates

@lru_cache()
def get_flare_data():
    return pd.read_csv(FLARE_DATA_PATH)

def get_random_flare_events(count, rng, threshold = 0):
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

    df = get_flare_data()
    # Filtering flares with fluxes below the threshold
    df_filtered = df[df['flux_amp'] >= threshold]
    df_filtered.reset_index(drop=True, inplace = True)
    indices = rng.integers(low = 0, high = len(df_filtered), size = count)
    print(len(df_filtered), ' out of ', len(df), ' flares are chosen for random selection based on threshold value of:', threshold)
    
    for index in indices:
        KIC.append(df_filtered['KIC'][index])
        St_time.append(df_filtered['St-BKJD'][index])
        End_time.append(df_filtered['End-BKJD'][index])
    
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

    return rng.normal(9000, 1000, count)

if __name__ == "__main__":
    # Getting Arguments
    argparser = argparse.ArgumentParser(
    description='Generates a LCLIB file with simulated flare instances')
    argparser.add_argument('flare_count', type = int,
                            help = 'Number of flares to be generated')
    argparser.add_argument('--file_name', type = str, required = False, default = 'LCLIB_Mdwarf-flare-LSST.TEXT',
                            help = 'Name of the output LCLIB file. Should have a .TEXT extension (Default: LCLIB_Mdwarf-flare-LSST.TEXT)')
    argparser.add_argument('--start_index', type = int, required = False, default = 0,
                            help = 'Use this if you want to start your file with an event number other than 0. LCLIB header is not added for start indices other than 0 (Default: 0)')
    argparser.add_argument('--remove_header', required = False, action = 'store_true',
                            help = 'Use this if you want to remove the LCLIB header. (Default: False)')
    argparser.add_argument('--generate_plots', required = False, action = 'store_true',
                            help = 'Use this if you want to save plots based on the simulations. Please note that this might have memory implications. Plotting is disabled by default (Default: False)')
    argparser.add_argument('--header_only', required = False, action = 'store_true',
                            help = 'Use this if you want only want to generate a LCLIB header. This does not generate any flares and thus cannot be used with --generate_plots to save plots. (Default: False)')

    args = argparser.parse_args()

    # Checking file name
    if not args.file_name.endswith(".TEXT"):
        print('Output file must be a .TEXT file')
        sys.exit(1)
    
    if args.header_only:
        with open(args.file_name, 'w') as output_file:
            add_LCLIB_header(args.flare_count, output_file)
        print('Created a header file. Exiting without flare simulation')
        sys.exit(1)
    else:
        # Starting flare modelling process
        start_time = time.time()
        run_generator(args.flare_count, args.file_name, args.start_index, args.remove_header, args.generate_plots)
        print("--- Simulations completed in %s seconds. File(s) saved. ---" % (int(time.time() - start_time)))
    
