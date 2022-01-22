import sqlite3
import pickle
import argparse
import sys

class MDwarfFlareDB:

    def __init__(self, db_path):

        self.db = sqlite3.connect(db_path)
        self.cur = self.db.cursor()

    def write_flares_to_LCLIB(self, lclib_path):

        with open(lclib_path, 'w') as output_file:
            for row in self.cur.execute('SELECT flare_object FROM flares'):
                flare = pickle.loads(row[0])
                flare.dump_flare_to_LCLIB(output_file)



def parse_args():
    
    # Getting Arguments
    argparser = argparse.ArgumentParser(
        description='Load a SQLite3 database file with simulated flare instances')

    argparser.add_argument('--db_path', type=str, required=True,
                        help='Path to the DB file.')
    argparser.add_argument('--output_file_path', type=str, required=True,
                        help='Name of the output LCLIB file. Should have a .TEXT extension')

    args = argparser.parse_args()
    return args

def main():

    # Arguments
    args = parse_args()

    if not args.output_file_path.endswith(".TEXT"):
        print('Output file must be a .TEXT file. Aborting process.')
        sys.exit(1)

    # Loading the database
    flare_db = MDwarfFlareDB(args.db_path)

    # Writing flares to the LCLIB
    flare_db.write_flares_to_LCLIB(args.output_file_path)

