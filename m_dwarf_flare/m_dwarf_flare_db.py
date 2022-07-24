from curses.ascii import NUL
from distutils import core
from itertools import count
import sqlite3
import pickle
import argparse
import sys
from tkinter import Frame
import numpy as np
import healpy as hp
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS
import astropy_healpix
from astropy import units as u
import matplotlib.pyplot as plt
from m_dwarf_flare.plotting_tools import plotGenricSkyMap
from m_dwarf_flare._version import version

class MDwarfFlareDB:

    def __init__(self, db_path):

        self.db = sqlite3.connect(db_path, timeout=10)
        self.db.execute('''CREATE TABLE IF NOT EXISTS flares
        (flare_index INTEGER, 
        ra REAL,
        dec REAL,
        flare_object BLOB)''')
        self.cur = self.db.cursor()

    def write_all_flares_to_LCLIB(self, lclib_path):

        with open(lclib_path, 'w') as output_file:
            for row in self.cur.execute('SELECT flare_object FROM flares ORDER BY flare_index'):
                flare = pickle.loads(row[0])
                flare.dump_flare_to_LCLIB(output_file)


    def get_confidence_interval_mask(self, skymap, confidence_interval):

        prob = skymap["PROB"]
        sorted_prob_index = np.argsort(prob)

        # Finding the cumulative probability distribution for sorted prob values
        cum_sorted_prob = np.cumsum(prob[sorted_prob_index])

        # Searching for the min probability pixel such that cumulatibe proba is still CI
        threshold_index = np.searchsorted(cum_sorted_prob, 1 - confidence_interval)

        # Setting all pixels to 0
        mask = np.zeros(len(skymap))
        
        # If the pixels lie in high probability region, we set them to 1
        mask[sorted_prob_index[threshold_index:]] = 1

        return mask, sorted_prob_index[threshold_index:]

    def add_NON_PERIODIC_LCLIB_header(self, count, output_file):
        """
        Function to write the header of the lclib file.
        """

        header = ('DOCUMENTATION:\n'
                '  PURPOSE: m Dwarf Flare model, Based on Kepler light curves and estimated distances from Gaia\n'
                '  REF:\n'
                '  - AUTHOR: Ved Shah\n'
                '  USAGE_KEY: GENMODEL\n'
                '  NOTES:\n'
                '  - M dwarf flare simulation based on Kepler data extrapolated for LSST\n'
                '  - Flare instances were taken from Yang et al. (2017)\n'
                '  - Distance data was taken from A Bailer Jones et al. (2021)\n'
                '  - Model Version number: {version_no}\n'
                '  PARAMS:\n'  
                '  - MWEBV - Milkyway E(B-V) from 3D dust model\n'
                '  - KIC_ID - Kepler Input Catalogue ID\n'
                '  - flare_temp_low - Temperature of the flare for spectral modelling with wavelength > balmer (in K)\n'
                '  - flare_temp_high - Temperature of the flare for spectral modelling with wavelength < balmer (in K)\n'
                '  - star_temp - Temperature of the star for spectral modelling (in K)\n'
                '  - distance - Distance to the star (in kpc)\n'
                '  - start_time - Start time of the reference flare (in BKJD)\n'
                '  - end_time - End time of the reference flare (in BKJD)\n'
                'DOCUMENTATION_END:\n\n'
                'SURVEY: LSST\n'
                'FILTERS: ugrizY\n'
                'MODEL: m-Dwarf-Flare-Model\n'
                'RECUR_TYPE: RECUR-NONPERIODIC\n'
                'MODEL_PARNAMES: MWEBV,KIC_ID,start_time,end_time,flare_temp_low,flare_temp_high,star_temp,distance.\n'
                'NEVENT: {count}\n\n').format(count = count, version_no = version)
        output_file.write(header)

    def get_flare_healpix_indices(self, nside):
        
        print('Getting coordinates from the DB....')
        ra, dec = zip(*self.cur.execute('SELECT ra, dec FROM flares'))
        print('Making coordinates objects....')
        coordinates = SkyCoord(ra=list(ra), dec=list(dec), frame=ICRS, unit='deg')

        print('Converting coordinates to helpix')
        map = astropy_healpix.HEALPix(nside, frame=ICRS, order="nested")
        hp_index = map.skycoord_to_healpix(coordinates, return_offsets=False)

        return hp_index

    def get_flares_in_skymap_ci(self, skymap_path, confidence_interval, lclib_path, to_plot):

        skymap = Table.read(skymap_path)
        nside = hp.npix2nside(len(skymap))

        mask, high_prob_flare_indices = self.get_confidence_interval_mask(skymap, confidence_interval)
        healpix_indices = self.get_flare_healpix_indices(nside)

        ra = []
        dec = []
        with open(lclib_path, 'w') as output_file:
            for i in range(len(healpix_indices)):
                if (healpix_indices[i] in high_prob_flare_indices):
                    for row in self.cur.execute('SELECT flare_object FROM flares WHERE flare_index = {} ORDER BY flare_index'.format(i)):
                        flare = pickle.loads(row[0])
                        flare.dump_flare_to_LCLIB(output_file)
                        if to_plot:
                            ra.append(flare.coordinates.ra)
                            dec.append(flare.coordinates.dec)

        if to_plot:

            plotGenricSkyMap(SkyCoord(ra=ra*u.deg, dec=dec*u.deg))

            hp.mollview(np.log10(skymap['PROB']), nest=True, min = -3, max = 0, title="Log Prob for GW event")
            hp.graticule(coord="E")

            hp.mollview(mask , nest=True, cmap='Greys',title="{}% Confidence interval healpix mask".format(confidence_interval * 100))
            hp.graticule(coord="E")

            plt.show()

    def get_flares_for_KN(self, skymap_path, confidence_interval, lclib_path, gw_trigger_time, survey_start_time, survey_end_time, window, to_plot):
        
        rng = np.random.default_rng(42)
        offsets = rng.uniform(low=0, high=window, size=1000000)

        skymap = Table.read(skymap_path, format='fits')
        nside = hp.npix2nside(len(skymap))

        print('Finding high CI healpix indices')
        mask, high_prob_flare_indices = self.get_confidence_interval_mask(skymap, confidence_interval)
        #high_prob_flare_indices = frozenset(high_prob_flare_indices)

        print('Getting the flare helpix indices')
        healpix_indices = self.get_flare_healpix_indices(nside)

        ra = []
        dec = []
        local_flare_count = 0

        with open(lclib_path, 'w') as output_file:
            
            print('Isolating Flares...')

            for i in range(len(healpix_indices)):
                if (healpix_indices[i] in high_prob_flare_indices):
                    for row in self.cur.execute('SELECT flare_object FROM flares WHERE flare_index = {} ORDER BY flare_index'.format(i)):
                        flare = pickle.loads(row[0])

                        # Relocating the flare to fall withing a certain period after the GW trigger
                        flare.relocate_flare_near_GW_trigger(gw_trigger_time, offsets[i])

                        # the -0.1 is a hack to make sure the flare starts before teh survey start
                        flare.add_survey_start_and_end_obsv(survey_start_time - 0.1, survey_end_time)
                        
                        flare.dump_flare_to_LCLIB(output_file)

                        local_flare_count += 1

                        if to_plot:
                            ra.append(flare.coordinates.ra)
                            dec.append(flare.coordinates.dec)

        header_name = lclib_path.split('.')[0] + '_HEADER.TEXT'
        print('Saving the header as' + header_name)
        # Write the LCLIB header
        with open(header_name, 'w') as header_file:
            self.add_NON_PERIODIC_LCLIB_header(local_flare_count, header_file)


        if to_plot:

            hp.mollview(skymap['PROB'], nest=True, title="Prob for GW event")
            hp.graticule(coord="E")

            hp.mollview(mask , nest=True, cmap='Greys',title="{}% Confidence interval healpix mask".format(confidence_interval * 100))
            hp.graticule(coord="E")

            plotGenricSkyMap(SkyCoord(ra=ra*u.deg, dec=dec*u.deg))

            plt.show()


def main():
    def parse_args_main():
    
        # Getting Arguments
        argparser = argparse.ArgumentParser(
            description='Load a SQLite3 database file with simulated flare instances')

        argparser.add_argument('db_path', type=str,
                            help='Path to the DB file.')
        argparser.add_argument('output_file_path', type=str,
                            help='Name of the output LCLIB file. Should have a .TEXT extension')

        args = argparser.parse_args()

        if not args.output_file_path.endswith(".TEXT"):
            print('Output file must be a .TEXT file. Aborting process.')
            sys.exit(1)

        return args

    # Arguments
    args = parse_args_main()

    # Loading the database
    flare_db = MDwarfFlareDB(args.db_path)

    # Writing flares to the LCLIB
    flare_db.write_all_flares_to_LCLIB(args.output_file_path)

def gw_event_localized_flares():
    def parse_args():
    
        # Getting Arguments
        argparser = argparse.ArgumentParser(
            description='Write files from db that lie within the CI of the skymap to a LCLIB file')

        argparser.add_argument('--db_path', type=str, required=True,
                            help='Path to the DB file.')
        argparser.add_argument('--output_file_path', type=str, required=True,
                            help='Path of the output LCLIB file. Should have a .TEXT extension')
        argparser.add_argument('--fits_file_path', type=str, required=True,
                            help='Path to the fits file of the GW event. Should be a single order fits file. Uses the same nside value as fit file for picking the flares using healpix')
        argparser.add_argument('--con_int', type=float, required=True,
                            help='Confidence interval of the area from which the flares are picked. CI should be between 0 and 1 inclusive')
        argparser.add_argument('--make_plots',  required=False, action='store_true',
                            help='Construct plots for GW events, the CI area mask and the flare objects skymap')
        args = argparser.parse_args()

        if not args.output_file_path.endswith(".TEXT"):
            print('Output file must be a .TEXT file. Aborting process.')
            sys.exit(1)
        if args.con_int < 0 or args.con_int > 1:
            print('CI should be between 0 and 1 inclusive. Aborting process.')
            sys.exit(1)
        return args

    # Arguments
    args = parse_args()

    flare_db = MDwarfFlareDB(args.db_path)
    flare_db.get_flares_in_skymap_ci(args.fits_file_path, args.con_int, args.output_file_path, args.make_plots)

def localize_flares_for_KN():
    def parse_args():
    
        # Getting Arguments
        argparser = argparse.ArgumentParser(
            description='Write files from db that lie within the CI of the skymap to a LCLIB file')

        argparser.add_argument('--db_path', type=str, required=True,
                            help='Path to the DB file.')
        argparser.add_argument('--output_file_path', type=str, required=True,
                            help='Path of the output LCLIB file. Should have a .TEXT extension')
        argparser.add_argument('--fits_file_path', type=str, required=True,
                            help='Path to the fits file of the GW (KN) event. Should be a single order fits file. Uses the same nside value as fit file for picking the flares using healpix')
        argparser.add_argument('--con_int', type=float, required=True,
                            help='Confidence interval of the area from which the flares are picked. CI should be between 0 and 1 inclusive')
        argparser.add_argument('--gw_trigger', type=float, required=True,
                            help='The gravitational trigger time in mjd.')
        argparser.add_argument('--survey_start', type=float, required=True,
                            help='The survey start time in mjd.')
        argparser.add_argument('--survey_end', type=float, required=True,
                            help='The survey end time in mjd.')
        argparser.add_argument('--window', type=float, required=False, default=1.0,
                            help='The flares will start within window number of days of the GW trigger [Default is 1 day].')
        argparser.add_argument('--make_plots',  required=False, action='store_true',
                            help='Construct plots for GW events, the CI area mask and the flare objects skymap')
        args = argparser.parse_args()

        if not args.output_file_path.endswith(".TEXT"):
            print('Output file must be a .TEXT file. Aborting process.')
            sys.exit(1)
        if args.con_int < 0 or args.con_int > 1:
            print('CI should be between 0 and 1 inclusive. Aborting process.')
            sys.exit(1)
        return args

    # Arguments
    args = parse_args()

    print("Connecting to the DB.")
    flare_db = MDwarfFlareDB(args.db_path)
    flare_db.get_flares_for_KN(args.fits_file_path, args.con_int, args.output_file_path, args.gw_trigger, args.survey_start, args.survey_end, args.window, args.make_plots)

def db_to_density_map():
    def parse_args():
    
        # Getting Arguments
        argparser = argparse.ArgumentParser(
            description='Use the DB to build a density distribution map stored as a fits file.')

        argparser.add_argument('db_path', type=str,
                            help='Path to the DB file.')
        argparser.add_argument('output_path', type=str,
                            help='Path of the output FITS file. Should have a .fits extension')
        argparser.add_argument('nside', type=int, 
                            help='Healpix NSIDE resolution that should be used to render the density map')
        args = argparser.parse_args()

        return args

    # Arguments
    args = parse_args()
    nside = args.nside

    flare_db = MDwarfFlareDB(args.db_path)
    healpix_indices = flare_db.get_flare_healpix_indices(nside)

    counts = np.zeros(hp.nside2npix(nside))
    for pixel in healpix_indices:
        counts[pixel] += 1
    prob = counts / np.sum(counts)
    print("Writing density map to FITS file...")
    hp.fitsfunc.write_map(args.output_path, prob, column_names=['PROB'], nest=True, coord='C', overwrite=True)
    print("Saved!")