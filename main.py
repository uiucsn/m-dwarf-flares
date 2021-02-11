from tools import *
import time
import concurrent.futures
import matplotlib.pyplot as plt


flare_data = astropy.io.ascii.read(FLARE_DATA_PATH, quotechar="\s")
kic = np.array(flare_data['KIC'])
list = np.unique(kic)

results = []

amplitude = []
area = []
duration = []

start = time.perf_counter()

def compute(KIC_ID):
    lc = load_light_curve(KIC_ID)
    # display_flare_plots(lc, KIC_ID)
    save_flare_instances(lc, KIC_ID)
    # save_flare_stats(lc, KIC_ID)
    amp_arr, area_array, duration_arr = all_flare_stats(lc, KIC_ID)
    return amp_arr, area_array, duration_arr

with concurrent.futures.ProcessPoolExecutor() as executor:
    results = executor.map(compute, list)
    for result in results:
        amplitude = amplitude + result[0]
        area = area + result[1]
        duration = duration + result[2]

end = time.perf_counter()

plot_all_flare_stats(amplitude, duration, area, 100)
print('Took '+str(end-start)+'s')