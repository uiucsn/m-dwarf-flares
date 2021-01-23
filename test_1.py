from flares import *
# KIC 892376


KIC_ID = input("Enter the KIC ID: ")
lc = loadLightCurve(KIC_ID)
#plotFlares(lc, KIC_ID)
saveFlareData(lc, KIC_ID)
getFlareStats(lc, KIC_ID)