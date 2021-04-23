from lc_tools import *
from spectra_tools import *
import matplotlib.pyplot as plt
import pandas as pd
from distance import *
from astropy.table import QTable

# KIC 892376
# KIC 1572802
# KIC 4355503 
# KIC 4470937 
# KIC 4726192 
# KIC 5016904 
# KIC 5597604

def generate_flare_file(KIC_ID, flare_temp, star_temp, distance, start_time, end_time):
    # Data extraction
    lc = load_light_curve(KIC_ID)
    flare_lc = get_flare_lc_from_time(lc, start_time, end_time)

    new_lc = get_normalized_lc(flare_lc)
    luminosity = get_luminosity_with_magnitude(KIC_ID).si  
    # Data modelling
    flare_luminosities = get_flare_luminosities_in_lsst_passbands(new_lc, KIC_ID, flare_temp, luminosity)
    
    baseline_luminosities = get_baseline_luminosity_in_lsst_passband(new_lc, KIC_ID, star_temp, luminosity)
    model_luminosities = fit_flare_on_base(flare_luminosities, baseline_luminosities)

    #model_flare_fluxes = get_spectra_data(new_lc, KIC_ID, flare_temp)
    model_mags = get_mags_in_lsst_passbands(model_luminosities, distance)
    add_LCLIB_header()
    dump_modeled_data_to_LCLIB(0, 0, 0, KIC_ID, flare_temp, star_temp, distance, start_time, end_time, model_mags)
    dump_modeled_data_to_LCLIB(0, 0, 0, KIC_ID, flare_temp, star_temp, distance, start_time, end_time, model_mags)

generate_flare_file(892376, 10000, 3000, 10, 909.013, 909.115)