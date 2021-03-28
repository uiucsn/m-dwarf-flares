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

LC_DATA_PATH = 'lc_data/KIC-{}.csv'
FLARE_DATA_PATH = 'flare_data/apjaa8ea2t3_mrt.txt'
FLARE_INSTANCES_PATH = 'flare_instances/KIC-{}/'
FLARE_INSTANCE = 'flare_instances/KIC-{kic}/{start}-{end}.csv'


def load_light_curve(KIC_ID):
    """
    Function for loading light curve data and constructing a kepler light curve.
    This function uses local data to construct the LightCurve, if possible. If 
    not, the data is loaded using lightkurve and a local copy is saved in csv form.
    Args:
        KIC_ID (int): The ID for the object from the Kepler survey. 
    Returns:
        Kepler Light Curve: A light curve for the object from Kepler.
    """

    KIC = "KIC {}".format(KIC_ID)
    if os.path.isfile(LC_DATA_PATH.format(KIC_ID)):
        print("Loading local light curve data ... " + str(KIC))
        df = pd.read_csv(LC_DATA_PATH.format(KIC_ID))
        lc = lk.LightCurve(time = df['time'], flux = df['flux'], flux_err = df['flux_err'])
    else:
        print('Downloading light curve data ... ' + str(KIC))
        lcf = search_lightcurvefile(KIC, mission = 'Kepler').download_all()
        new_lcf =  LightCurveFileCollection([x for x in lcf if x.targetid == int(KIC_ID)])
        lc = new_lcf.SAP_FLUX.stitch()
        lc.to_csv(LC_DATA_PATH.format(KIC_ID))
    return lc

def save_flare_instances(lc, KIC_ID):
    """
    This function saves every instance of a flare for a given Kepler Input Catalog 
    object in csv form.  It looks at flare instances from the apjaa8ea2t3_mrt.txt file.
    Args:
        lc (Kepler Light Curve Object): A light curve for the object from Kepler.
        KIC_ID (int): The ID for the object from the Kepler survey. 
    """

    # Reading flare data file.
    flare_data = astropy.io.ascii.read(FLARE_DATA_PATH, quotechar="\s")

    #Plotting points identified as flares.
    for flare in flare_data:
        if flare['KIC'] != int(KIC_ID):
            continue
            
        # Obtaining start and end indices for flare instances
        flare_lc = get_flare_lc(lc,flare)

        if not os.path.isdir(FLARE_INSTANCES_PATH.format(KIC_ID)):
            os.mkdir(FLARE_INSTANCES_PATH.format(KIC_ID))

        # Saving indiviual flare instances as csv files
        flare_lc.to_csv(FLARE_INSTANCE.format(kic = KIC_ID, start = str(flare['St-BKJD']), end = str(flare['St-BKJD'])))

def all_flare_stats(lc, KIC_ID):
    """
    Returns the flare statistics for all the flares belonging to a particular KIC.

    Args:
        lc (Kepler Light Curve Object): A light curve for the object from Kepler.
        KIC_ID (int): The ID for the object from the Kepler survey. 

    Returns:
        tuple : tuple contatining the list of amplitude, area, duration, data for all the flares in a particular KIC object
    """
    flare_data = astropy.io.ascii.read(FLARE_DATA_PATH, quotechar="\s")
    duration = []
    area = []
    amplitude = []
    data = []
    for flare in flare_data:
        if flare['KIC'] != int(KIC_ID):
            continue

        flare_lc = get_flare_lc(lc,flare)
        flux_arr = [flux - np.amin(flare_lc.flux) for flux in flare_lc.flux]

        amp = np.amax(flux_arr)

        duration.append(flare['End-BKJD'] - flare['St-BKJD'])
        area.append(flare['Area'])
        amplitude.append(amp)
        data.append((str(KIC_ID),flare['St-BKJD'],flare['End-BKJD']))
    return amplitude, area, duration, data

def display_flare_plots(lc, KIC_ID):
    """
    This function plots the flares for a given object from the Kepler Input Catalog.
    It looks at flare instances from the apjaa8ea2t3_mrt.txt file.  
    Args:
        lc (Kepler Light Curve Object): A light curve for the object from Kepler.
        KIC_ID (int): The ID for the object from the Kepler survey. 
    """

    # Reading flare data file.
    flare_data = astropy.io.ascii.read(FLARE_DATA_PATH, quotechar="\s")

    #Plotting points identified as flares.
    for flare in flare_data:

        if flare['KIC'] != int(KIC_ID):
            continue
        fig, (ax1, ax2, ax3) = plt.subplots(3)
        fig.suptitle('Flare from {start} to {end}'.format(start = str(flare['St-BKJD']), end = str(flare['End-BKJD'])))

        flare_lc = get_flare_lc(lc,flare)

        # Full light curve plot
        ax1.plot(lc.time, lc.flux)
        ax1.set_ylabel('Relative flux')
        ax1.set_xlim([flare['St-BKJD']-50,flare['End-BKJD']+50])
        ax1.axvspan(flare['St-BKJD'], flare['End-BKJD'], color='red', alpha=0.2)
        ax1.set_title('LC for KIC {}'.format(str(KIC_ID)))

        # Light curve of the flare
        ax2.plot(lc.time, lc.flux, linestyle='-', marker='.')
        ax2.axvspan(flare['St-BKJD'], flare['End-BKJD'], color='red', alpha=0.2)
        ax2.set_ylabel('Relative flux')
        ax2.set_xlim([flare['St-BKJD']-2,flare['End-BKJD']+2])

        ax3.plot(flare_lc.time, flare_lc.flux, linestyle='-', marker='.')
        ax3.axvspan(flare['St-BKJD'], flare['End-BKJD'], color='red', alpha=0.2)
        ax3.set_xlabel('Time BKJD')
        ax3.set_ylabel('Relative flux')

        plt.show()

def save_flare_stats(lc, KIC_ID):
    """
    This function plots histograms for flare duration, area and amplitude for
    a given Kelper Input Catalog Object. It also saves the plots. It looks at 
    flare instances from the apjaa8ea2t3_mrt.txt file.
    Args:
        lc (Kepler Light Curve Object): A light curve for the object from Kepler.
        KIC_ID (int): The ID for the object from the Kepler survey. 
    """

    flare_data = astropy.io.ascii.read(FLARE_DATA_PATH, quotechar="\s")
    
    duration = []
    flareArea = []
    amplitude = []

    for flare in flare_data:
        if flare['KIC'] != int(KIC_ID):
            continue

        flare_lc = get_flare_lc(lc,flare)
        amp = np.amax(flare_lc.flux)

        duration.append(flare['End-BKJD'] - flare['St-BKJD'])
        flareArea.append(flare['Area'])
        amplitude.append(amp)

    num_bins = 10

    # Plotting flare amplitude stats.
    fig1, (ax1, ax4) = plt.subplots(2,1,figsize=(15, 12))
    fig1.suptitle('Flare Amplitude Statistics for KIC {}'.format(KIC_ID))

    N_amplitude, bins_amplitude, patches_amplitude = ax1.hist(amplitude, bins=num_bins, alpha = 0.5) 
    ax1.set_title('Normalized Amplitude distribution')
    ax1.set_xlabel('Amplitude')
    ax1.set_ylabel('Number of flares')

    amplitude = np.sort(amplitude)
    cumulativeAmp = np.arange(amplitude.size, 0, -1)
    ax4.plot(amplitude, cumulativeAmp)
    ax4.set_title('Normalized Amplitude distribution (Cumulative)')
    ax4.set_xlabel('Amplitude')
    ax4.set_ylabel('Number of flares with amplitude less than')

    # Plotting flare duration stats.
    fig2, (ax2, ax5) = plt.subplots(2,1,figsize=(15, 12))
    fig2.suptitle('Flare Duration Statistics for KIC {}'.format(KIC_ID))

    N_duration, bins_duration, patches_duration = ax2.hist(duration, bins=num_bins, alpha = 0.5) 
    ax2.set_title('Flare Duration distribution')
    ax2.set_xlabel('Duration in days')
    ax2.set_ylabel('Number of flares')
    
    duration = np.sort(duration)
    cumulativeDuration = np.arange(duration.size, 0, -1)
    ax5.plot(duration, cumulativeDuration) 
    ax5.set_title('Flare Duration distribution (Cumulative)')
    ax5.set_xlabel('Duration in days')
    ax5.set_ylabel('Number of flares with duration less than')


    # Plotting flare area stats.
    fig3, (ax3, ax6) = plt.subplots(2,1,figsize=(15, 12))
    fig3.suptitle('Flare Area Statistics for KIC {}'.format(KIC_ID))

    N_area, bins_area, patches_area = ax3.hist(flareArea, bins=num_bins, alpha=0.5) 
    ax3.set_title('Flare Area Distribution')
    ax3.set_xlabel('Flare Area')
    ax3.set_ylabel('Number of flares')

    flareArea = np.sort(flareArea)
    cumulativeArea = np.arange(flareArea.size, 0, -1)
    ax6.plot(amplitude, cumulativeArea)
    ax6.set_title('Flare Area Distribution (Cumulative)')
    ax6.set_xlabel('Flare Area')
    ax6.set_ylabel('Number of flares with area less than')

    if not os.path.isdir('obj_stats/KIC-{}/'.format(KIC_ID)):
        os.mkdir('obj_stats/KIC-{}/'.format(KIC_ID))

    fig1.savefig('obj_stats/KIC-{}/amplitude'.format(KIC_ID))
    fig2.savefig('obj_stats/KIC-{}/duration'.format(KIC_ID))
    fig3.savefig('obj_stats/KIC-{}/area'.format(KIC_ID))

def plot_all_flare_stats(amplitude, duration, area, num_bins):
    """
    Plots and saves the amplitude, duration and area statistics for all the amplitudes, durations and flare areas from all the flares 
    mentioned in the apjaa8ea2t3_mrt.txt file.

    Args:
        amplitude (list): A list contatining all the amplitudes
        duration (list): A list containing all the durations
        area (list): A list containing all the flare Areas
        num_bins (int): Number of bins for the histograms.
    """
    # Plotting flare amplitude stats.
    fig1, ax1 = plt.subplots(1,figsize=(15, 12))

    amplitude = np.array(amplitude)
    amplitude = np.sort(amplitude)

    N_amplitude, bins_amplitude, patches_amplitude = ax1.hist(amplitude, bins=num_bins, alpha = 0.5) 
    ax1.set_title('Normalized Amplitude distribution')
    ax1.set_yscale('log')
    ax1.set_xlabel('Amplitude')
    ax1.set_ylabel('Number of flares')


    # Plotting flare duration stats.
    fig2, (ax2, ax5) = plt.subplots(2,1,figsize=(15, 12))

    duration = np.array(duration)
    duration = np.sort(duration)

    N_duration, bins_duration, patches_duration = ax2.hist(duration, bins=num_bins, alpha = 0.5) 
    ax2.set_title('Flare Duration distribution')
    ax2.set_yscale('log')
    ax2.set_xlabel('Duration in days')
    ax2.set_ylabel('Number of flares')


    cumulativeDuration = np.arange(duration.size, 0, -1)
    ax5.plot(duration, cumulativeDuration) 
    ax5.set_yscale('log')
    ax5.set_title('Flare Duration distribution (Cumulative)')
    ax5.set_xlabel('Duration in days')
    ax5.set_ylabel('Number of flares with duration less than')


    # Plotting flare area stats.
    fig3, (ax3, ax6) = plt.subplots(2,1,figsize=(15, 12))

    area = np.array(area)
    area = np.sort(area)

    N_area, bins_area, patches_area = ax3.hist(area, bins=num_bins, alpha=0.5) 
    ax3.set_title('Flare Area Distribution')
    ax3.set_yscale('log')
    ax3.set_xlabel('Flare Area')
    ax3.set_ylabel('Number of flares')

    cumulativeArea = np.arange(area.size, 0, -1)
    #ax6.plot(amplitude, cumulativeArea)
    ax6.set_yscale('log')
    ax6.set_title('Flare Area Distribution (Cumulative)')
    ax6.set_xlabel('Flare Area')
    ax6.set_ylabel('Number of flares with area less than')

    fig1.savefig('amplitude')
    fig2.savefig('duration')
    fig3.savefig('area')

    plt.show()

def get_flare_lc(lc, flare):
    """
    This function takes the Light Curve for a KIC object along with a flare instance
    and returns a normalized light curve containing just the flare with one point prior
    to and after it.
    Args:
        lc (Kepler light Curve object): A light curve for the object from Kepler.
        flare (flare instance): A flare instances from the apjaa8ea2t3_mrt.txt file
    Returns:
        Kepler light curve object: A normalized light curve of the flare followed by a 
        trailing and leading observation
    """
    start_index = find_nearest_index(lc.time, flare['St-BKJD']) - 1
    end_index = find_nearest_index(lc.time, flare['End-BKJD']) + 2
    
    if(start_index < 0):
        start_index = 0
    

    flare_time = lc.time[start_index:end_index]
    flare_flux = lc.flux[start_index:end_index]
    flare_err = lc.flux_err[start_index:end_index]

    flare_lc = lk.LightCurve(time = flare_time, flux = flare_flux, flux_err = flare_err) 
    return flare_lc

def get_flare_lc_from_time(lc, start_time, end_time):
    """
    This function takes the Light Curve for a KIC object along with a flare instance
    and returns a normalized light curve containing just the flare with one point prior
    to and after it.
    Args:
        lc (Kepler light Curve object): A light curve for the object from Kepler.
        flare (flare instance): A flare instances from the apjaa8ea2t3_mrt.txt file
    Returns:
        Kepler light curve object: A normalized light curve of the flare followed by a 
        trailing and leading observation
    """
    start_index = find_nearest_index(lc.time, start_time) - 1
    end_index = find_nearest_index(lc.time, end_time) + 2
    
    if(start_index < 0):
        start_index = 0
    

    flare_time = lc.time[start_index:end_index]
    flare_flux = lc.flux[start_index:end_index]
    flare_err = lc.flux_err[start_index:end_index]

    flare_lc = lk.LightCurve(time = flare_time, flux = flare_flux, flux_err = flare_err) 
    return flare_lc

def get_normalized_lc(lc):
    """
    This function takes the Light Curve for a KIC object along with a flare instance
    and returns a normalized light curve containing just the flare with one point prior
    to and after it.
    Args:
        lc (Kepler light Curve object): A light curve for the object from Kepler.
        flare (flare instance): A flare instances from the apjaa8ea2t3_mrt.txt file
    Returns:
        Kepler light curve object: A normalized light curve of the flare followed by a 
        trailing and leading observation
    """

    flare_time = lc.time
    flare_flux = lc.flux
    flare_err = lc.flux_err
    
    min_flux = np.amin(flare_flux)
    flare_flux = [flux - min_flux for flux in flare_flux]

    flare_lc = lk.LightCurve(time = flare_time, flux = flare_flux, flux_err = flare_err) 
    return flare_lc

def find_nearest_index(lc_time, value):
    """
    This function takes a sorted array and a value and returns the index where
    the value is found. If it can't find the value in the array, it returns the index of 
    a neigbour value which is closest to it's magnitude (in absolute terms).
    Args:
        lc_time (numpy array): The array where the index for the value is to be found
        value (float): The value to be found in the array
    Returns:
        int : Index of the value in the array or the index of it's closes neigbour.
    """
    index = np.searchsorted(lc_time, value, side="left")

    if index > 0 and index == len(lc_time):
        return index - 1
    if index > 0 and value - lc_time[index-1] < lc_time[index] - value:
        return index - 1
    else:
        return index

