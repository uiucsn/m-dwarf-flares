import argparse
from ast import arg
import os
import time
import sys

from m_dwarf_flare.generator import run_generator
from m_dwarf_flare.lc_tools import add_LCLIB_header


def parse_args():
    # Getting Arguments
    argparser = argparse.ArgumentParser(
        description='Generates a LCLIB file with simulated flare instances')

    argparser.add_argument('--flare_count', type=int, required=True,
                           help='Number of flares to be generated.')
    argparser.add_argument('--spectrum_class', type=str, required=True, choices=['bb_simple', 'bb_balmer_jump'],
                           help='Type of model used for flare spectral modeling. bb_simple or bb_balmer_jump are currently supported.')
    argparser.add_argument('--lc_data', type=str, default='lc_data',
                           help='Path to the directory with Kepler light curve, which can be downloaded by "mdwarf-download-data" script')
    argparser.add_argument('--dir_name', type=str, required=True, default='sample',
                           help='Path to the directory to store all the simulation data. Directory will be created if it does not exist (Default: sample)')
    argparser.add_argument('--file_name', type=str, required=True, default='sample.TEXT',
                           help='Name of the output LCLIB file. Should have a .TEXT extension (Default: LCLIB_Mdwarf-flare-LSST.TEXT)')
    argparser.add_argument('--use_dpf', required=False, action='store_true',
                           help='Use this if you want to use differential photometry for filtering. Standard filtering is used by default. (Default: False)')
    argparser.add_argument('--pickle_sims', required=False, action='store_true',
                           help='Use this if you want to store all the nominal flare simulation objects. (Default: Write to LCLIB only)')
    argparser.add_argument('--generate_plots', required=False, action='store_true',
                           help='Use this if you want to save plots based on the simulations. Please note that this might have memory implications. Plotting is disabled by default (Default: False)')
    argparser.add_argument('--start_index', type=int, required=False, default=0,
                           help='Use this if you want to start your file with an event number other than 0. LCLIB header is not added for start indices other than 0 (Default: 0)')
    argparser.add_argument('--remove_header', required=False, action='store_true',
                           help='Use this if you want to remove the LCLIB header. (Default: False)')
    argparser.add_argument('--header_only', required=False, action='store_true',
                           help='Use this if you want only want to generate a LCLIB header. This does not generate any flares and thus cannot be used with --generate_plots to save plots. (Default: False)')
    argparser.add_argument('--save_db', required=False, action='store_true',
                           help='Use this if you want want to store your flare objects in a sqlite3 database (Default: False)')

    args = argparser.parse_args()
    return args


def main():
    args = parse_args()

    # Checking file name
    if not args.file_name.endswith(".TEXT"):
        print('Output file must be a .TEXT file. Aborting simulation process.')
        sys.exit(1)

    # Creating a directory to store everything
    os.makedirs(args.dir_name, exist_ok=True)

    if args.header_only:
        path = os.path.join(args.dir_name, args.file_name)
        with open(path, 'w') as output_file:
            add_LCLIB_header(args.flare_count, output_file)
        print('Created a header file. Exiting without flare simulation')
        sys.exit(1)
    else:
        # Starting flare modelling process
        start_time = time.time()
        run_generator(flare_count=args.flare_count,
                      spectrum_type=args.spectrum_class,
                      lc_data_path=args.lc_data,
                      dir_path=args.dir_name, file_path=args.file_name,
                      use_dpf=args.use_dpf,
                      pickle_sims=args.pickle_sims,
                      generate_plots=args.generate_plots,
                      start_index=args.start_index,
                      remove_header=args.remove_header,
                      save_db=args.save_db)
        print("--- Simulations completed in %s seconds. File(s) saved. ---" % (int(time.time() - start_time)))


if __name__ == '__main__':
    main()