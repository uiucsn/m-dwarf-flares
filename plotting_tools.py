import astropy.coordinates as coord
from astropy.coordinates import ICRS, FK5, SkyCoord
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
from astropy.table import Table
from astropy.coordinates import SkyCoord, CylindricalRepresentation
from astropy import units as u
import astropy
import numpy as np
from pyvo.dal import TAPService
from astroquery.vizier import Vizier
import math
import pandas as pd
from astropy.table import QTable

def plotSkyMapFromSDSSData():

    data = fits.getdata('dr7_lowmass_av.fits', 1)

    t = Table(data)
    c = CylindricalRepresentation(rho = t['RGAL'], phi = t['THETA'] * u.rad, z = t['ZGAL'])
    coords = SkyCoord(c)

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection="mollweide")
    scatter = ax.scatter(-coords.ra.wrap_at(180 * u.deg).radian, coords.dec.wrap_at(180 * u.deg).radian, s=3, vmin=0)
    ax.grid(True)
    ax.set_xticklabels(['10h', '8h', '6h', '4h', '2h', '0h', '22h', '20h', '18h', '16h', '14h'])
    plt.show()

def plotSkyMapFromSUPERBLINK():

    df = astropy.io.ascii.read('aj403664t1_mrt.txt', quotechar="\s")
    coords = SkyCoord(ra = df['RAdeg'], dec = df['DEdeg'])
    distance = 1.0 / df['plx']

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection="mollweide")
    scatter = ax.scatter(-coords.ra.wrap_at(180 * u.deg).radian, coords.dec.wrap_at(180 * u.deg).radian,c = distance, s=3, vmin=0)
    ax.grid(True)
    ax.set_xticklabels(['10h', '8h', '6h', '4h', '2h', '0h', '22h', '20h', '18h', '16h', '14h'])
    plt.colorbar(scatter, label = "Distance in pc")
    plt.show()

def plotLBDistributionSUPERBLINK():

    df = astropy.io.ascii.read('aj403664t1_mrt.txt', quotechar="\s")
    coord = SkyCoord(df['RAdeg'], df['DEdeg'], frame='icrs', unit='deg')
    galactic = coord.transform_to('galactic')

    fig1, (ax1, ax2, ax3) = plt.subplots(1,3)

    ax1.hist(galactic.b.value, bins=100, alpha = 0.5) 
    ax1.set_title('Distribution of m dwarfs as a function of b')
    ax1.set_xlabel('b in degrees')
    ax1.set_ylabel('Number of m dwarfs')

    ax2.hist(galactic.l.value, bins=100, alpha = 0.5) 
    ax2.set_title('Distribution of m dwarfs as a function of b')
    ax2.set_xlabel('l in degrees')
    ax2.set_ylabel('Number of m dwarfs')

    ax3.hist2d(galactic.l.value, galactic.b.value, bins = 20)
    ax3.set_xlabel("Galactic Longitude in degrees")
    ax3.set_ylabel("Galactic Latitude in degrees")

    plt.show()

def plotRaDecDistributionSUPERBLINK():

    df = astropy.io.ascii.read('aj403664t1_mrt.txt', quotechar="\s")
    coord = SkyCoord(df['RAdeg'], df['DEdeg'], frame='icrs', unit='deg')

    fig1, (ax1, ax2, ax3) = plt.subplots(1,3)

    ax1.hist(coord.ra.value, bins=100, alpha = 0.5) 
    ax1.set_title('Distribution of m dwarfs as a function of b')
    ax1.set_xlabel('RA in degrees')
    ax1.set_ylabel('Number of m dwarfs')

    ax2.hist(coord.dec.value, bins=100, alpha = 0.5) 
    ax2.set_title('Distribution of m dwarfs as a function of b')
    ax2.set_xlabel('Dec in degrees')
    ax2.set_ylabel('Number of m dwarfs')

    ax3.hist2d(coord.ra.value, coord.dec.value, bins = 20)
    ax3.set_xlabel("RA in degrees")
    ax3.set_ylabel("Dec in degrees")

    plt.show()

def plotDistanceDistributionSUPERBLINK():

    table = astropy.io.ascii.read('aj403664t1_mrt.txt')
    distance = np.sort(1.0 / table[(table['plx'] < 9) & (~table['plx'].mask)]['plx'])

    fig1, ax1 = plt.subplots(1,1)
    N_amplitude, bins_amplitude, patches_amplitude = ax1.hist(distance, bins=100, alpha = 0.5) 
    ax1.set_title('Distribution of m dwarfs as a function of distance')
    ax1.set_xlabel('distance in pc')
    ax1.set_ylabel('Number of m dwarfs')
    plt.show()

def plot_effective_kepler_temps():

    df = pd.read_csv('eff_temp.csv')
    print(np.mean(df['teff']))
    print(np.std(df['teff']))
    count, bins, ignored = plt.hist(df['teff'], bins=10, alpha = 0.5) 
    plt.xlabel('Effective Temperature in Kelvin')
    plt.ylabel('Number of m dwarfs')
    plt.show()





