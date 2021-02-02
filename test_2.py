from flares import *
import time
import concurrent.futures

list = [892376, 1572802, 4355503, 4470937, 4726192, 5016904]

start = time.perf_counter()

def compute(KIC_ID):
    lc = load_light_curve(KIC_ID)
    #get_flare_plots(lc, KIC_ID)
    save_flare_instances(lc, KIC_ID)
    save_flare_stats(lc, KIC_ID)

with concurrent.futures.ProcessPoolExecutor() as executor:
    results = executor.map(compute, list)

end = time.perf_counter()
print('Took '+str(end-start)+'s')