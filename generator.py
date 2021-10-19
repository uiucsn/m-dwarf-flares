from functools import lru_cache, partial
from astropy import units as u
from astropy.modeling.models import BlackBody

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

# Minimum Relative Flux Amplitude of the flares. Flares below this relative flux amplitude will be filtered out. This is a performance optimization.
MIN_RELATIVE_FLUX_AMPLITUDE = 0.01

# Maximum gap between two flare data points in days.
MAX_TIME_GAP_IN_FLARE = 0.1

# Maximum magnitude for a flare in the LSST u passband. Flares above this mag value will be filtered out.
PEAK_MAGNITUDE_TOLERANCE = 0.5
PEAK_MAGNITUDE_THRESHOLD = {
    'u': 23.66 + PEAK_MAGNITUDE_TOLERANCE,
    'g': 24.69 + PEAK_MAGNITUDE_TOLERANCE,
    'r': 24.06 + PEAK_MAGNITUDE_TOLERANCE,
    'i': 23.45 + PEAK_MAGNITUDE_TOLERANCE,
    'z': 22.54 + PEAK_MAGNITUDE_TOLERANCE,
    'y': 21.62 + PEAK_MAGNITUDE_TOLERANCE,
}
# Minimum magnitude amplitude of the simulated flare in the u passband. Flares below this mag amplitude will be filtered out.
BAND_AMPLITUDE_THRESHOLD = 0.2

# Generating new parameters is an expernsive process. One way to speed it up is 
PARAMETER_COUNT_MULTIPLIER = 75

def run_generator(flare_count, file_path, start_index, remove_header, to_plot, use_dpf, spectrum_type):
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
                
                print("4. Building star spectrum functions ...")
                star_spectrum_functions = get_star_spectrum_function(parameter_count, rng)
                
                print("5. Building flare spectrum functions ...")
                flare_spectrum_functions = get_flare_spectrum_function(spectrum_type, parameter_count, rng)

                print("6. Commencing flare modelling ...")
                for i in range(parameter_count):
                    bar.update(number_of_nominal_flares)
                    number_of_simulated_flares += 1 
                    # Breaking out of the loop if the correct number of flares are generated
                    if (number_of_nominal_flares == flare_count):
                        break
                    extinction = {
                        'u': extinction_values['u'][i],
                        'g': extinction_values['g'][i],
                        'r': extinction_values['r'][i],
                        'i': extinction_values['i'][i],
                        'z': extinction_values['z'][i],
                        'y': extinction_values['y'][i],
                    }
                    is_valid_flare, modeled_flare = generate_model_flare_file(start_index + number_of_nominal_flares, coordinates[i], galactic_coordinates[i], distances[i], kic_id[i], start_time[i], end_time[i], star_spectrum_functions[i], flare_spectrum_functions[i], extinction, output_file, use_dpf)
                    if is_valid_flare:
                        number_of_nominal_flares += 1
                        if to_plot:
                            nominal_flare_indices.append(i)
                            nominal_flare_instance.append(modeled_flare)
    output_file.close()
    print(int((flare_count * 100) / number_of_simulated_flares),'%','of the simulated flares passed the threshold cuts')
    if to_plot:
        nominal_coordinates = coord.SkyCoord(coordinates[nominal_flare_indices])
        print("7. Generating plots ...")
        save_simulation_plots(nominal_coordinates, nominal_flare_instance, rng)


def generate_model_flare_file(index, coordinates, galactic_coordinates, distance, KIC_ID, start_time, end_time, star_spectrun_function, flare_spectrum_function, extinction, output_file, use_dpf):
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
    flare_luminosities = get_flare_luminosities_in_lsst_passbands(new_lc, KIC_ID, flare_spectrum_function, luminosity)
    baseline_luminosities = get_baseline_luminosity_in_lsst_passband(new_lc, KIC_ID, star_spectrun_function, luminosity)
    model_luminosities = fit_flare_on_base(flare_luminosities, baseline_luminosities)

    # Converting luminosity to distances based on distances and applying extinction
    model_mags = get_mags_in_lsst_passbands(model_luminosities, distance)
    model_mags_with_extinction = apply_extinction_to_lsst_mags(model_mags, extinction)

    if is_nominal_flare(model_mags_with_extinction, use_dpf) and flare_has_nominal_cadence(model_mags_with_extinction):
        # Writing modelled data to LCLIB file if the flare is nominal

        star_temp = star_spectrun_function.keywords['temp']

        if 'temp' in flare_spectrum_function.keywords.keys():
            # Assigning the same high and low temps for a simple blacbody mmodel
            flare_temp_low = flare_spectrum_function.keywords['temp']
            flare_temp_high = flare_spectrum_function.keywords['temp']
        else:
            # Assigning the different high and low temps for a blacbody mmodel with balmer jump
            flare_temp_low = flare_spectrum_function.keywords['temp_low']
            flare_temp_high = flare_spectrum_function.keywords['temp_high']

        dump_modeled_data_to_LCLIB(index, galactic_coordinates.l, galactic_coordinates.b, KIC_ID, start_time, end_time, star_temp, flare_temp_low, flare_temp_high, distance, model_mags_with_extinction, output_file)
        return True, model_mags_with_extinction
    else:
        return False, model_mags_with_extinction

def flare_has_nominal_cadence(flare):
    """
    Checking if the flare has nominal cadence. Computes the time gap between every pair of observations
    in a flare and returns true if the maximum time gap between two data points is less than 
    MAX_TIME_GAP_IN_FLARE, false otherwise.

    Args:
        flare ([dictionary of lc's]): A dictionary of lc's in LSST passbands.

    Returns:
        [boolean]: A boolean if the flare cadence is normal or not.
    """
    
    time = flare['u'].time
    diff = time[1:] - time[:-1]

    return np.max(diff) < MAX_TIME_GAP_IN_FLARE

def is_nominal_flare(flare, use_dpf):
    """
    Checking if the generated flare is under the threshold for max magnitude and has an
    amplitude greater than the threshold

    Args:
        flare (dicitonary): A dictionary of lightcurves in the LSST ugrizy and kep passbands

    Returns:
        [boolean]: True if the flare is nominal, False otherwise.
    """

    passbands = ['u','g','r','i','z','y']

    if use_dpf:
        min_mag = {}
        max_mag = {}
        for passband in passbands:
            min_mag[passband] = np.amin(flare[passband].flux)
            max_mag[passband] = np.amax(flare[passband].flux)

        dict = {passband: (max_mag[passband] - 
                    2.5 * math.log10(
                        10 ** (
                            0.4 * (max_mag[passband] - min_mag[passband])
                        ) 
                        - 1
                    )
                ) for passband in passbands}

        return any(dict[passband] <= PEAK_MAGNITUDE_THRESHOLD[passband] for passband in passbands)
    
    else:
        # Checking if all bands have atleast one magnitude value greater than the correspnding PEAK_MAGNITUDE_THRESHOLD values. 
        peak_mag_is_bright = (np.any(flare[passband].flux <= PEAK_MAGNITUDE_THRESHOLD[passband]) for passband in passbands)

        # Checking if all passbands have amplitude greater than the BAND_AMPLITUDE_THRESHOLD
        ampl_is_high = (np.ptp(flare[passband].flux) > BAND_AMPLITUDE_THRESHOLD for passband in passbands)

        # Returns true if there is atleast one passband that statisfies both of the above conditions
        return any(peak and ampl for peak, ampl in zip(peak_mag_is_bright, ampl_is_high))

def get_star_spectrum_function(parameter_count, rng):

    star_temp = get_normally_distributed_star_temp(parameter_count, rng)
    function_wrappers = []
    for i in range(parameter_count):
        function_wrappers.append(partial(build_spectrum_bb_simple, temp=star_temp[i]))
    return function_wrappers

def get_flare_spectrum_function(spectrum_type, parameter_count, rng):

    if spectrum_type == 'bb_simple':
        flare_temp = get_normally_distributed_flare_temp_low(parameter_count, rng)
        function_wrappers = []
        for i in range(parameter_count):
            function_wrappers.append(partial(build_spectrum_bb_simple, temp=flare_temp[i]))
        return function_wrappers

    elif spectrum_type == 'bb_balmer_jump':
        flare_temp_low = get_normally_distributed_flare_temp_low(parameter_count, rng)
        flare_temp_high = get_normally_distributed_flare_temp_high(parameter_count, rng)
        function_wrappers = []
        for i in range(parameter_count):
            function_wrappers.append(partial(build_spectrum_bb_balmer_jump, temp_low=flare_temp_low[i], temp_high=flare_temp_high[i]))
        return function_wrappers

    raise ValueError(f'Spectrum type {spectrum_type} is not supported')

def build_spectrum_bb_simple(lmbd, temp):

    bb = BlackBody(temperature=temp*u.K)
    return bb(lmbd)

def build_spectrum_bb_balmer_jump(lmbd, temp_low, temp_high):

    cutoff_wavelenght = 3645 * u.AA
    bb_low = BlackBody(temperature=temp_low*u.K)
    bb_high = BlackBody(temperature=temp_high*u.K)
    flux = np.where(lmbd > cutoff_wavelenght, bb_low(lmbd), bb_high(lmbd))
    return flux

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
    N_M_DWARF = 24.8 * (10e9)

    # Fraction of simulated flares that fall within the LSST thresholds
    LSST_VISIBILITY_FRACTION = .08

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


def get_normally_distributed_flare_temp_low(count, rng):
    """
    Returns a numpy array of flare temperatures modelled after a normal distribution.

    Args:
        count (int): Length of numpy array to be returned

    Returns:
        numpy array: numpy array containing the flare temperatures with length = count
    """

    return rng.normal(9000, 1000, count)

def get_normally_distributed_flare_temp_high(count, rng):
    """
    Returns a numpy array of flare temperatures modelled after a normal distribution.

    Args:
        count (int): Length of numpy array to be returned

    Returns:
        numpy array: numpy array containing the flare temperatures with length = count
    """

    return rng.normal(35000, 3000, count)

if __name__ == "__main__":
    # Getting Arguments
    argparser = argparse.ArgumentParser(
    description='Generates a LCLIB file with simulated flare instances')
    argparser.add_argument('flare_count', type = int,
                            help = 'Number of flares to be generated')
    argparser.add_argument('flare_spectrum_type', type = str,
                            help = 'Type of model used for flare spectral modeling. bb_simple or bb_balmer_jump are currently supported.')
    argparser.add_argument('--file_name', type = str, required = False, default = 'LCLIB_Mdwarf-flare-LSST.TEXT',
                            help = 'Name of the output LCLIB file. Should have a .TEXT extension (Default: LCLIB_Mdwarf-flare-LSST.TEXT)')
    argparser.add_argument('--use_dpf', required = False, action = 'store_true',
                            help = 'Use this if you want to use differential photometry for filtering. Standard filtering is used by default. (Default: False)')
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
        print('Output file must be a .TEXT file. Aborting simulation process.')
        sys.exit(1)
    
    # Checking if the spectral model is supported.
    if args.flare_spectrum_type not in ['bb_simple', 'bb_balmer_jump']:
        print('The spectrum model passed is not supported. Aborting simulation process.')    
        sys.exit(1)
    
    if args.header_only:
        with open(args.file_name, 'w') as output_file:
            add_LCLIB_header(args.flare_count, output_file)
        print('Created a header file. Exiting without flare simulation')
        sys.exit(1)
    else:
        # Starting flare modelling process
        start_time = time.time()
        run_generator(args.flare_count, args.file_name, args.start_index, args.remove_header, args.generate_plots, args.use_dpf,  args.flare_spectrum_type)
        print("--- Simulations completed in %s seconds. File(s) saved. ---" % (int(time.time() - start_time)))
    
