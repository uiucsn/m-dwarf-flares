import astropy
from astropy.io import ascii
import lightkurve as lk
from lightkurve import search_lightcurvefile
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import os.path
import pandas as pd
import csv
import numpy as np

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
        print("Loading local light curve data ...")
        df = pd.read_csv(LC_DATA_PATH.format(KIC_ID))
        lc = lk.LightCurve(time = df['time'], flux = df['flux'], flux_err = df['flux_err'])
    else:
        print('Downloading light curve data ...')
        lcf = search_lightcurvefile(KIC).download_all()
        lc = lcf.SAP_FLUX.stitch()
        lc.to_csv(LC_DATA_PATH.format(KIC_ID))
    return lc

def get_flare_plots(lc, KIC_ID):
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

        start_index = np.searchsorted(lc.time, flare['St-BKJD']) - 2
        end_index = np.searchsorted(lc.time, flare['End-BKJD']) + 1

        flare_time = lc.time[start_index:end_index]
        flare_flux = lc.flux[start_index:end_index]
        flare_err = lc.flux_err[start_index:end_index]

        min_flux = np.amin(flare_flux)
        # Use numpy instead 
        flare_flux = [flux - min_flux for flux in flare_flux]

        flare_lc = lk.LightCurve(time = flare_time, flux = flare_flux, flux_err = flare_err)

        # Full light curve plot
        ax1.plot(lc.time, lc.flux)
        ax1.set_ylabel('Relative flux')
        ax1.set_xlim([flare['St-BKJD']-200,flare['End-BKJD']+200])
        ax1.axvspan(flare['St-BKJD'], flare['End-BKJD'], color='red', alpha=0.2)
        ax1.set_title('LC for KIC {}'.format(str(KIC_ID)))

        # Light curve of the flare
        ax2.plot(lc.time, lc.flux, linestyle='-', marker='.')
        ax2.axvspan(flare['St-BKJD'], flare['End-BKJD'], color='red', alpha=0.2)
        ax2.set_xlabel('Time BKJD')
        ax2.set_ylabel('Relative flux')
        ax2.set_xlim([flare['St-BKJD']-2,flare['End-BKJD']+2])
        ax2.set_title('Flare Area: {}'.format(str(flare['Area'])))

        ax3.plot(flare_lc.time, flare_lc.flux, linestyle='-', marker='.')
        ax3.axvspan(flare['St-BKJD'], flare['End-BKJD'], color='red', alpha=0.2)

        plt.show()

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

        # Make this part a function that takes the KIC and the start and end date and run the function in parallel to save all the flare instances faster.
            
        # Obtaining start and end indices for flare instances
        start_index = np.searchsorted(lc.time, flare['St-BKJD']) - 2
        end_index = np.searchsorted(lc.time, flare['End-BKJD']) + 1

        flare_time = lc.time[start_index:end_index]
        flare_flux = lc.flux[start_index:end_index]
        flare_err = lc.flux_err[start_index:end_index]

        min_flux = np.amin(flare_flux)
        flare_flux = [flux - min_flux for flux in flare_flux]

        flare_lc = lk.LightCurve(time = flare_time, flux = flare_flux, flux_err = flare_err) 
        flare_lc.plot()

        if not os.path.isdir(FLARE_INSTANCES_PATH.format(KIC_ID)):
            os.mkdir(FLARE_INSTANCES_PATH.format(KIC_ID))

        # Saving indiviual flare instances as csv files
        flare_lc.to_csv(FLARE_INSTANCE.format(kic = KIC_ID, start = str(flare['St-BKJD']), end = str(flare['St-BKJD'])))

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

        start_index = np.searchsorted(lc.time, flare['St-BKJD']) - 2
        end_index = np.searchsorted(lc.time, flare['End-BKJD']) + 1

        flare_flux = lc.flux[start_index:end_index]
        flare_time = lc.time[start_index:end_index]

        min_flux = np.amin(flare_flux)
        flare_flux = [flux - min_flux for flux in flare_flux]
        amp = np.amax(flare_flux)

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