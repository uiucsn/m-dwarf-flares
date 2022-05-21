import pickle
import matplotlib.pyplot as plt
import os
import numpy as np
import sqlite3
import json

# 4 is the highest pickle protocol for Python 3.7, which is the minimum Python version we support. 
# This is to ensure that future python versions can make use of the pickled objects.
PICKLE_PROTOCOL = 4

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

    def pickle_flare_instance(self, dir_path):

        # Making a directory if it does not already exist
        pickle_file_dir = os.path.join(dir_path, 'flare_objects')
        os.makedirs(pickle_file_dir, exist_ok=True)
        
        filename = os.path.join(pickle_file_dir, f'flare_sim_{self.index}.pkl')
        with open(filename, 'wb') as outp:  # Overwrites any existing file.
            pickle.dump(self, outp, PICKLE_PROTOCOL)


    def plot_flare_instance(self):

        for passband, lc in self.lightcurves.items():
            plt.plot(lc.time, lc.flux, label=passband)
        plt.legend('Passband')
        plt.show()

    def dump_flare_to_LCLIB(self, output_file):

        event_marker = "#------------------------------\n"
        start = "START_EVENT: {}\n".format(self.index)
        end = "END_EVENT: {}\n".format(self.index)
        nrow = "NROW: {nrow} l: {l:.5f} b: {b:.5f}.\n".format(nrow = len(self.lightcurves['kep'].time), 
                                                                l = self.galactic_coordinate.l.value, 
                                                                b = self.galactic_coordinate.b.value)
        parameters = "PARVAL: {KIC_ID} {start} {end} {f_temp_low:.2f} {f_temp_high:.2f} {s_temp:.2f} {dist:.7f}\n".format(KIC_ID = self.kic_id, 
                                                                                        f_temp_low = self.flare_temp_low,
                                                                                        f_temp_high = self.flare_temp_high,  
                                                                                        s_temp = self.star_temp, 
                                                                                        dist = self.distance.value, 
                                                                                        start = self.start_time, 
                                                                                        end = self.end_time)
        angle_match = "ANGLEMATCH_b: {angle_match:.3f}\n".format(angle_match = np.max([5, 0.5 * np.abs(self.galactic_coordinate.b.value)]))                                                                                    
        readings = ""

        # For loop to add readings of the simulations to the text file
        for i in range(len(self.lightcurves['kep'].flux)):
            if i == 0:
                readings += "T: "
            else:
                readings += "S: "
            readings += "{time:>10.5f} {u:>10.3f} {g:>10.3f} {r:>10.3f} {i:>10.3f} {z:>10.3f} {y:>10.3f}\n".format(time = self.lightcurves['kep'].time[i], 
                                                                                                                    kep = self.lightcurves['kep'].flux[i], 
                                                                                                                    u = self.lightcurves['u'].flux[i], 
                                                                                                                    g = self.lightcurves['g'].flux[i], 
                                                                                                                    r = self.lightcurves['r'].flux[i], 
                                                                                                                    i = self.lightcurves['i'].flux[i], 
                                                                                                                    z = self.lightcurves['z'].flux[i], 
                                                                                                                    y = self.lightcurves['y'].flux[i]) 

        simulation = event_marker + start + nrow + parameters + angle_match + readings + end
        output_file.write(simulation)
    
    def relocate_flare_near_GW_trigger(self, GW_trigger_time, offset):

        for passband in self.lightcurves:

            lc = self.lightcurves[passband]

            # Translating the flare to start 
            lc.time -= lc.time[0]

            # Flare starts exactly at the GW trigger time
            lc.time += GW_trigger_time

            # Adding an offset, such that the flare starts at
            # 'offset' days after the GW trigger time.
            lc.time += offset

            # Replacing the lc
            self.lightcurves[passband] = lc

    def add_survey_start_and_end_obsv(self, survey_start_time, survey_end_time):

        for passband in self.lightcurves:

            lc = self.lightcurves[passband]

            # Adding observation before the flare at survey_start_time
            lc.time = np.insert(lc.time, 0, survey_start_time, axis=0)
            lc.flux = np.insert(lc.flux, 0, lc.flux[0], axis=0)

            # Adding observation after the flare at survey_end_time
            lc.time = np.append(lc.time, survey_end_time)
            lc.flux = np.append(lc.flux, lc.flux[-1])

            # Replacing the lc
            self.lightcurves[passband] = lc
            

    def save_flare_to_sqlite_db(self, db):

        bytes = pickle.dumps(self, protocol=PICKLE_PROTOCOL)
        obj = sqlite3.Binary(bytes)

        db.execute("INSERT INTO flares VALUES (?,?,?,?)", 
        (self.index, self.coordinates.ra.deg, self.coordinates.dec.deg, obj))

        db.commit()        