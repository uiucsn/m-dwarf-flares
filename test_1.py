from flares import *
# KIC 892376
# KIC 1572802
# KIC 4355503 
# KIC 4470937 
# KIC 4726192 
# KIC 5016904 
# KIC 5597604


KIC_ID = input("Enter the KIC ID: ")
lc = loadLightCurve(KIC_ID)
#plotFlares(lc, KIC_ID)
saveFlareData(lc, KIC_ID)
getFlareStats(lc, KIC_ID)