import astropy
from astropy.io import ascii
import lightkurve as lk
from lightkurve import search_targetpixelfile
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
    pixelfile = search_targetpixelfile(KIC, quarter=16).download();
    lc = pixelfile.to_lightcurve(aperture_mask='all');
    lc.to_csv("lc_data/KIC-" + KIC_ID + '.csv')

corrected_lc = lc.remove_nans().flatten()

# Reading flare data file.
flare_data = astropy.io.ascii.read('flare_data/apjaa8ea2t3_mrt.txt', quotechar="\s")

#Plotting points identified as flares.
for flare in flare_data:
    if flare['KIC'] == int(KIC_ID):
        if flare['St-BKJD'] > 1520:
            plt.plot(corrected_lc.time, corrected_lc.flux, linestyle='-', marker='o')
            plt.xlim([flare['St-BKJD']-0.01,flare['End-BKJD']+0.01])
            plt.show()
