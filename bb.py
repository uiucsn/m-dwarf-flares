from tools import *
import matplotlib.pyplot as plt
# KIC 892376
# KIC 1572802
# KIC 4355503 
# KIC 4470937 
# KIC 4726192 
# KIC 5016904 
# KIC 5597604


KIC_ID = input("Enter the KIC ID: ")
lc = load_light_curve(KIC_ID)
get_spectra_data(lc, KIC_ID)