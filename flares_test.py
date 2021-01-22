import astropy
from astropy.io import ascii
import lightkurve as lk
from lightkurve import search_lightcurvefile
import matplotlib.pyplot as plt
import os.path
import pandas as pd
import csv
import numpy as np
# KIC 892376

def loadLightCurve(KIC_ID):
    KIC = "KIC " + KIC_ID
    if os.path.isfile("lc_data/KIC-" + KIC_ID + '.csv'):
        print("Loading local light curve data ...")
        df = pd.read_csv("lc_data/KIC-" + KIC_ID + '.csv')
        lc = lk.LightCurve(time = df['time'], flux = df['flux'], flux_err = df['flux_err'])
        return lc
    else:
        print('Downloading light curve data ...')
        lcf = search_lightcurvefile(KIC).download_all()
        lc = lcf.SAP_FLUX.stitch()
        lc.to_csv("lc_data/KIC-" + KIC_ID + '.csv')
        return lc

def plotFlares(lc):
    # Reading flare data file.
    flare_data = astropy.io.ascii.read('flare_data/apjaa8ea2t3_mrt.txt', quotechar="\s")

    #Plotting points identified as flares.
    for flare in flare_data:
        if flare['KIC'] == int(KIC_ID):
            fig, (ax1, ax2, ax3) = plt.subplots(3)
            fig.suptitle('Flare from ' + str(flare['St-BKJD']) + ' to ' + str(flare['End-BKJD']))

            start_index = np.searchsorted(lc.time, flare['St-BKJD']) - 2
            end_index = np.searchsorted(lc.time, flare['End-BKJD']) + 1

            flare_time = lc.time[start_index:end_index]
            flare_flux = lc.flux[start_index:end_index]
            flare_err = lc.flux_err[start_index:end_index]

            min_flux = np.amin(flare_flux)
            flare_flux = [flux - min_flux for flux in flare_flux]

            flare_lc = lk.LightCurve(time = flare_time, flux = flare_flux, flux_err = flare_err)
  
            # Full light curve plot
            ax1.plot(lc.time, lc.flux)
            ax1.set_ylabel('Relative flux')
            ax1.set_xlim([flare['St-BKJD']-200,flare['End-BKJD']+200])
            ax1.axvspan(flare['St-BKJD'], flare['End-BKJD'], color='red', alpha=0.2)
            ax1.set_title('LC for KIC ' + str(KIC_ID))

            # Light curve of the flare
            ax2.plot(lc.time, lc.flux, linestyle='-', marker='.')
            ax2.axvspan(flare['St-BKJD'], flare['End-BKJD'], color='red', alpha=0.2)
            ax2.set_xlabel('Time BKJD')
            ax2.set_ylabel('Relative flux')
            ax2.set_xlim([flare['St-BKJD']-2,flare['End-BKJD']+2])
            ax2.set_title('Flare Area: '+str(flare['Area']))

            ax3.plot(flare_lc.time, flare_lc.flux, linestyle='-', marker='.')
            ax3.axvspan(flare['St-BKJD'], flare['End-BKJD'], color='red', alpha=0.2)

            plt.show()

def saveFlareData(lc):
    # Reading flare data file.
    flare_data = astropy.io.ascii.read('flare_data/apjaa8ea2t3_mrt.txt', quotechar="\s")

    #Plotting points identified as flares.
    index = 0
    for flare in flare_data:
        if flare['KIC'] == int(KIC_ID):
            
            start_index = np.searchsorted(lc.time, flare['St-BKJD']) - 2
            end_index = np.searchsorted(lc.time, flare['End-BKJD']) + 1

            flare_time = lc.time[start_index:end_index]
            flare_flux = lc.flux[start_index:end_index]
            flare_err = lc.flux_err[start_index:end_index]

            min_flux = np.amin(flare_flux)
            flare_flux = [flux - min_flux for flux in flare_flux]

            flare_lc = lk.LightCurve(time = flare_time, flux = flare_flux, flux_err = flare_err) 

            if not os.path.isdir("flare_instances/KIC-" + KIC_ID + '/'):
                os.mkdir("flare_instances/KIC-" + KIC_ID + '/')

            flare_lc.to_csv("flare_instances/KIC-" + KIC_ID + '/' + str(flare['St-BKJD']) + '-' + str(flare['St-BKJD']) + '.csv')
            
KIC_ID = input("Enter the KIC ID: ")
lc = loadLightCurve(KIC_ID)
#plotFlares(lc)
saveFlareData(lc)
