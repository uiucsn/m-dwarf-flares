import astropy
from astropy.io import ascii
import lightkurve as lk
from lightkurve import search_lightcurvefile
import matplotlib.pyplot as plt
import os.path
import pandas as pd
import csv
# KIC 892376

KIC_ID = input("Enter the KIC ID: ")
KIC = "KIC " + KIC_ID

# Obtaining LightCurve.
if os.path.isfile("lc_data/KIC-" + KIC_ID + '.csv'):
    print("Loading local light curve data ...")
    df = pd.read_csv("lc_data/KIC-" + KIC_ID + '.csv')
    lc = lk.LightCurve(time = df['time'], flux = df['flux'], flux_err = df['flux_err'])
else:
    print('Downloading light curve data ...')
    lcf = search_lightcurvefile(KIC).download_all()
    lc = lcf.SAP_FLUX.stitch()
    lc.to_csv("lc_data/KIC-" + KIC_ID + '.csv')


# Reading flare data file.
flare_data = astropy.io.ascii.read('flare_data/apjaa8ea2t3_mrt.txt', quotechar="\s")

#Plotting points identified as flares.
for flare in flare_data:
    if flare['KIC'] == int(KIC_ID):
        fig, (ax1, ax2) = plt.subplots(2)
        fig.suptitle('Flare from ' + str(flare['St-BKJD']) + ' to ' + str(flare['End-BKJD']))

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

        plt.show()
