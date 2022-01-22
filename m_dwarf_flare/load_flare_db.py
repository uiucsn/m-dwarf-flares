import sqlite3
import pickle
import argparse
# import m_dwarf_flare.m_dwarf_flare
import sys

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


def write_flares_to_LCLIB():

    # Arguments
    args = parse_args()

    if not args.output_file_path.endswith(".TEXT"):
        print('Output file must be a .TEXT file. Aborting process.')
        sys.exit(1)

    # Loading the database
    db = sqlite3.connect(args.db_path)
    cur = db.cursor()
    with open(args.output_file_path, 'w') as output_file:
        for row in cur.execute('SELECT flare_object FROM flares'):
            flare = pickle.loads(row[0])
            flare.dump_flare_to_LCLIB(output_file)

if __name__ == '__main__':
    write_flares_to_LCLIB()