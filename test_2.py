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
flare_data = astropy.io.ascii.read(FLARE_DATA_PATH, quotechar="\s")
kic = np.array(flare_data['KIC'])
list = np.unique(kic)

#df = pd.read_csv("kepler_kic_v10.csv.gz", compression='gzip')

table = QTable(names = ('KIC ID' , 'kic_kepmag'))

""" print(df.columns)
for i in range(len(df)):
    for kic in list: 
        if(df['kic_kepler_id'][i] == kic):
            print(df['kic_kepler_id'][i], df['kic_kepmag'][i])
            table.add_row((df['kic_kepler_id'][i], df['kic_kepmag'][i]))
            table.write('mag.csv', format = 'ascii.csv') """

def generate_flare_file(KIC_ID, temperature, start_time, end_time):
    lc = load_light_curve(KIC_ID)
    flare_lc = get_flare_lc_from_time(lc, start_time, end_time)
    new_lc = get_normalized_lc(flare_lc)
    get_spectra_data(new_lc, KIC_ID, temperature)
    luminosity = get_luminosity_with_magnitude(KIC_ID).si
    flare_luminosities = get_flare_luminosities_in_lsst_passbands(new_lc, KIC_ID, temperature, luminosity)
    baseline_luminosities = get_baseline_luminosity_in_lsst_passband(new_lc, KIC_ID, 2000, luminosity)
    fit_flare_on_base(flare_luminosities, baseline_luminosities)
KIC_ID = input("Enter the KIC ID: ")
generate_flare_file(KIC_ID, 10000, 909.013, 909.115)