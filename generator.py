import numpy as np
import astropy
from astropy import units as u
import astropy.coordinates as coord
import astropy.io.ascii 
import math
import pandas as pd
from lc_tools import *
from spectra_tools import *
import matplotlib.pyplot as plt
from distance import *
from astropy.table import QTable

FLARE_DATA_PATH = 'filtered_flares.csv'

TOTAL_KEPLER_M_DWARF_COUNT = 4664
TOTAL_KEPLER_DURATION_IN_DAYS = 1460
TOTAL_FLARE_INSTANCES_COUNT = 103187

LOCAL_RADIUS = 20
NUMBER_OF_LOCAL_M_DWARFS = 1082

KEPLER_MEAN_EFFECTIVE_TEMP_FOR_M_DWARFS = 3743.4117647058824
KEPLER_STD_EFFECTIVE_TEMP_FOR_M_DWARFS = 161.37182827551771

def run_generator():
    RADIUS = 5
    DAYS = 365

    print("Computing number of flare instances")
    flare_count = get_number_of_expected_flares(RADIUS, DAYS)
    print("Sampling coordinates of stars")
    coordinates, distances = get_uniformly_distributed_spherical_coordinates(RADIUS, flare_count)
    print("Obtaining reference flare")
    kic_id, start_time, end_time = get_random_flare_events(flare_count)
    print("Sampling star temperature")
    star_temp = get_normally_distributed_star_temp(flare_count)
    print("Sampling flare temperature")
    flare_temp = get_normally_distributed_flare_temp(flare_count)
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

    dump_modeled_data_to_LCLIB(index, coordinates.ra, coordinates.dec, KIC_ID, flare_temp, star_temp, distance, start_time, end_time, model_mags)
   
def get_number_of_expected_flares(radius, duration):

    flare_density = (NUMBER_OF_LOCAL_M_DWARFS / (4/3 * math.pi * LOCAL_RADIUS**2)) * (TOTAL_FLARE_INSTANCES_COUNT / (TOTAL_KEPLER_M_DWARF_COUNT * TOTAL_KEPLER_DURATION_IN_DAYS))
    volume = 4/3 * math.pi * (radius ** 3)
    flare_count = math.floor(flare_density * volume * duration)

    print('Flare Density:', flare_density, 'flares/pc^3/day')
    print('Expected Number of flares in {radius} pc sphere during the {duration} day period:'.format(radius = radius, duration = duration), flare_count)

    return flare_count

def get_uniformly_distributed_spherical_coordinates(radius, count):

    coordinates = []
    dist = []

    while len(coordinates) != count:
        print((len(coordinates) / count) * 100,'% Complete')
        x, y, z = (np.random.uniform(-1 * radius,radius, 3) * u.pc)
        if (x ** 2 + y ** 2 + z ** 2) ** (0.5) < radius * u.pc:
            r, lat, lon = coord.cartesian_to_spherical(x,y,z)
            c = coord.SkyCoord(ra = lon, dec = lat)
            coordinates.append(c)
            dist.append(r.value)
    
    return coordinates, dist

def get_random_flare_events(count):
    KIC = []
    St_time = []
    End_time = []

    df = pd.read_csv(FLARE_DATA_PATH)
    indices = np.random.randint(-1, len(df), count)
    
    for index in indices:
        KIC.append(df['KIC'][index])
        St_time.append(df['St-BKJD'][index])
        End_time.append(df['End-BKJD'][index])
    
    return KIC, St_time, End_time


def get_normally_distributed_star_temp(count):
    return np.random.normal(KEPLER_MEAN_EFFECTIVE_TEMP_FOR_M_DWARFS, KEPLER_STD_EFFECTIVE_TEMP_FOR_M_DWARFS, count)

def get_normally_distributed_flare_temp(count):
    return np.random.normal(9000, 500, count)

run_generator()