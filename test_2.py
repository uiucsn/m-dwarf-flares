from flares import *
import time

list = [892376, 1572802, 4355503, 4470937, 4726192, 5016904]

start = time.perf_counter()
for KIC_ID in list:
    lc = load_light_curve(KIC_ID)
    #get_flare_plots(lc, KIC_ID)
    save_flare_instances(lc, KIC_ID)
    save_flare_stats(lc, KIC_ID)
end = time.perf_counter()
print('Took '+str(end-start)+'s')