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
        plt.plot(lc.time, lc.flux, linestyle='-', marker='o')
        plt.axvspan(flare['St-BKJD'], flare['End-BKJD'], color='red', alpha=0.2)
        plt.xlim([flare['St-BKJD']-0.5,flare['End-BKJD']+0.5])
        plt.show()
