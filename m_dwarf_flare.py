import pickle
import matplotlib.pyplot as plt
from pygit2 import Repository
import os

PICKLE_PROTOCOL = 4
GIT_COMMIT_NUMBER = Repository(os.path.dirname(os.path.abspath(__file__))).head.peel().hex

class MDwarfFlare:

    def __init__(self, index, lightcurves, coordinates, galactic_coordinates, distance, kic_id, start_time, end_time, star_spectrum_function, flare_spectrum_function, extinction):
        self.index = index
        self.lightcurves = lightcurves
        self.coordinates = coordinates
        self.galactic_coordinate = galactic_coordinates
        self.distance = distance
        self.kic_id = kic_id
        self.start_time = start_time
        self.end_time = end_time

        self.extinction = extinction


        self.star_temp = star_spectrum_function.keywords['temp']

        if 'temp' in flare_spectrum_function.keywords:
            # Assigning the same high and low temps for a simple blacbody mmodel
            self.flare_temp_low = flare_spectrum_function.keywords['temp']
            self.flare_temp_high = flare_spectrum_function.keywords['temp']
        else:
            # Assigning the different high and low temps for a blacbody mmodel with balmer jump
            self.flare_temp_low = flare_spectrum_function.keywords['temp_low']
            self.flare_temp_high = flare_spectrum_function.keywords['temp_high']
        
        self.git_commit = GIT_COMMIT_NUMBER

    def pickle_flare_instance(self, dir_path):

        # Making a directory if it does not already exist
        os.makedirs(dir_path, exist_ok=True)
        filename = os.path.join(f'flare_sim_{self.index}.pkl')
        with open(filename, 'wb') as outp:  # Overwrites any existing file.
            pickle.dump(self, outp, PICKLE_PROTOCOL)


    def plot_flare_instance(self):

        for passband, lc in self.lightcurves.items():
            plt.plot(lc.time, lc.flux, label=passband)
        plt.legend('Passband')
        plt.show()
