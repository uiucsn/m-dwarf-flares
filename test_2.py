from flares import *

list = [892376, 1572802, 4355503, 4470937, 4726192, 5016904, 5597604]

for KIC_ID in list:
    lc = loadLightCurve(KIC_ID)
    #plotFlares(lc, KIC_ID)
    getFlareStats(lc, KIC_ID)
    saveFlareData(lc, KIC_ID)