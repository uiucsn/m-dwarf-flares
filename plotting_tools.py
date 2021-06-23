from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord, CylindricalRepresentation
from astropy import units as u
from ch_vars.spatial_distr import MilkyWayDensityJuric2008 as MWDensity

import astropy
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def save_simulation_plots(coordinates, flares, rng):
    """
    Generates and saves the plots for the simulated flares that are saved to the LCLIB file.

    Args:
        coordinates (Astropy sky coords): Array of coords that are saved to the LCLIB
        flares (dictionary of lightcurves): Array of flares containing magnitude time series in lsst passbands.
        rng (numpy rng): Random number generator
    """

    save_peak_mag_distribution(flares)
    save_mag_amp_distribution(flares)
    save_simulated_distance_distribution(coordinates, rng)
    save_simulated_l_distribution(coordinates, rng)
    save_simulated_b_distribution(coordinates, rng)
    save_sky_map_with_distances(coordinates)

def save_peak_mag_distribution(flares):
    """
    Generatres and saves a plot for distribution of the Peak magnitudes in all 6 LSST passbands 

    Args:
        flares (dictionary of lightcurves): Array of flares containing magnitude time series in lsst passbands.
    """
    
    u_peaks = []
    g_peaks = []
    r_peaks = []
    i_peaks = []
    z_peaks = []
    y_peaks = []

    for flare in flares:
        u_peaks.append(np.amax(flare['u'].flux))
        g_peaks.append(np.amax(flare['g'].flux))
        r_peaks.append(np.amax(flare['r'].flux))
        i_peaks.append(np.amax(flare['i'].flux))
        z_peaks.append(np.amax(flare['z'].flux))
        y_peaks.append(np.amax(flare['y'].flux))
    
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot()
    ax.set_title('Distribution peak magnitudes in LSST passbands')
    ax.hist(u_peaks, bins=20, histtype='step', facecolor='m', label = 'u band') 
    ax.hist(g_peaks, bins=20, histtype='step', facecolor='g', label = 'g band') 
    ax.hist(r_peaks, bins=20, histtype='step', facecolor='r', label = 'r band') 
    ax.hist(i_peaks, bins=20, histtype='step', facecolor='c', label = 'i band') 
    ax.hist(z_peaks, bins=20, histtype='step', facecolor='b', label = 'z band') 
    ax.hist(y_peaks, bins=20, histtype='step', facecolor='y', label = 'y band') 
    ax.set_xlabel('Peak Magnitude')
    ax.set_ylabel('Number of m dwarfs')
    plt.legend(loc='upper right')
    plt.savefig("simulation_stats/peak_mag_distribution.pdf")

def save_mag_amp_distribution(flares):
    """
    Generatres and saves a plot for distribution of the magnitude amplitudes in all 6 LSST passbands 

    Args:
        flares (dictionary of lightcurves): Array of flares containing magnitude time series in lsst passbands.
    """

    u_amps = []
    g_amps = []
    r_amps = []
    i_amps = []
    z_amps = []
    y_amps = []

    for flare in flares:
        u_amps.append(np.amax(flare['u'].flux) - np.amin(flare['u'].flux))
        g_amps.append(np.amax(flare['g'].flux) - np.amin(flare['g'].flux))
        r_amps.append(np.amax(flare['r'].flux) - np.amin(flare['r'].flux))
        i_amps.append(np.amax(flare['i'].flux) - np.amin(flare['i'].flux))
        z_amps.append(np.amax(flare['z'].flux) - np.amin(flare['z'].flux))
        y_amps.append(np.amax(flare['y'].flux) - np.amin(flare['y'].flux))

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot()
    ax.set_title('Distribution magnitude amplitudes in LSST passbands')
    ax.hist(u_amps, bins=20, histtype='step', facecolor='m', label = 'u band') 
    ax.hist(g_amps, bins=20, histtype='step', facecolor='g', label = 'g band') 
    ax.hist(r_amps, bins=20, histtype='step', facecolor='r', label = 'r band') 
    ax.hist(i_amps, bins=20, histtype='step', facecolor='c', label = 'i band') 
    ax.hist(z_amps, bins=20, histtype='step', facecolor='b', label = 'z band') 
    ax.hist(y_amps, bins=20, histtype='step', facecolor='y', label = 'y band') 
    ax.set_xlabel('Magnitude amplitude')
    ax.set_ylabel('Number of m dwarfs')
    plt.legend(loc='upper right')
    plt.savefig("simulation_stats/mag_amp_distribution.pdf")


def save_simulated_distance_distribution(coordinates, rng):
    """
    Generatres and saves a plot for distribution of the simulated distances vs the MW modelled distances.

    Args:
        coordinates (Astropy Sky coords): Array of coords saved to the LCLIB file.
        rng (Numpy rng): Random number generator.
    """

    mw = MWDensity()
    ideal_coordinates = SkyCoord(mw.sample_eq(len(coordinates), rng))
    ideal_coordinates = ideal_coordinates.galactic

    gal = coordinates.galactic
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot()
    ax.set_title('Distribution of distances for simulated flares vs MW model')
    ax.hist(gal.distance.value, bins=20, histtype='step', label = 'Simulation generated') 
    ax.hist(ideal_coordinates.distance.value, bins=20, histtype='step', label = 'MW Model generated')
    ax.set_xlabel('Distance in kpc')
    ax.set_ylabel('Number of m dwarfs')
    plt.legend(loc='upper right')
    plt.savefig("simulation_stats/dist_distribution.pdf")

def save_simulated_l_distribution(coordinates, rng):
    """
    Generatres and saves a plot for distribution of the simulated galactic longitude vs the MW modelled galactic longitude.

    Args:
        coordinates (Astropy Sky coords): Array of coords saved to the LCLIB file.
        rng (Numpy rng): Random number generator.
    """

    mw = MWDensity()
    ideal_coordinates = SkyCoord(mw.sample_eq(len(coordinates), rng))
    ideal_coordinates = ideal_coordinates.galactic

    gal = coordinates.galactic
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot()
    ax.set_title('Distribution of galactic longitude for simulated flares vs MW model')
    ax.hist(gal.l.value, bins=20, histtype='step', label = 'Simulation generated') 
    ax.hist(ideal_coordinates.l.value, bins=20, histtype='step', label = 'MW Model generated')
    ax.set_xlabel('Galactic longitude in degrees')
    ax.set_ylabel('Number of m dwarfs')
    plt.legend(loc='upper right')
    plt.savefig("simulation_stats/l_distribution.pdf")

def save_simulated_b_distribution(coordinates, rng):
    """
    Generatres and saves a plot for distribution of the simulated galactic latitude vs the MW modelled galactic latitude.

    Args:
        coordinates (Astropy Sky coords): Array of coords saved to the LCLIB file.
        rng (Numpy rng): Random number generator.
    """

    mw = MWDensity()
    ideal_coordinates = SkyCoord(mw.sample_eq(len(coordinates), rng))
    ideal_coordinates = ideal_coordinates.galactic

    gal = coordinates.galactic
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot()
    ax.set_title('Distribution of galactic latitude for simulated flares vs MW model')
    ax.hist(gal.b.value, bins=20, histtype='step', label = 'Simulation generated') 
    ax.hist(ideal_coordinates.b.value, bins=20, histtype='step', label = 'MW Model generated')
    ax.set_xlabel('Galactic latitude in degrees')
    ax.set_ylabel('Number of m dwarfs')
    plt.legend(loc='upper right')
    plt.savefig("simulation_stats/b_distribution.pdf")

def save_sky_map_with_distances(coordinates):
    """
    Plots the sky map fo simulated m dwarfs with distances which are written to the LCLIB

    Args:
        coordinates (Astropy Sky coords): Coordintaes of the m dwarf flares.
    """

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection="mollweide")
    ax.set_title('Skymap for simulated m dwarf flare instances')
    scatter = ax.scatter(-coordinates.ra.wrap_at(180 * u.deg).radian, coordinates.dec.wrap_at(180 * u.deg).radian,c = coordinates.distance.value, s=3, vmin=0)
    ax.grid(True)
    ax.set_xticklabels(['10h', '8h', '6h', '4h', '2h', '0h', '22h', '20h', '18h', '16h', '14h'])
    plt.colorbar(scatter, label = "Distance in kpc")
    plt.savefig("simulation_stats/sky_map.pdf")

def plotSkyMapFromSUPERBLINK():
    """
    Plots the Sky Map for m dwarf population based on SUPERBLINK Survey
    """

    df = astropy.io.ascii.read('data_files/aj403664t1_mrt.txt', quotechar="\s")
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
    """
    Plots a histogram for distribution of stars as a function of galactic latitude, longitude
    and a 2D histogram heat map based on data from the SUPERBLINK survey.
    """

    df = astropy.io.ascii.read('data_files/aj403664t1_mrt.txt', quotechar="\s")
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

    ax3.hist2d(galactic.l.value, galactic.b.value, cmap = 'plasma', bins = 20)
    ax3.set_xlabel("Galactic Longitude in degrees")
    ax3.set_ylabel("Galactic Latitude in degrees")

    plt.show()

def plotRaDecDistributionSUPERBLINK():
    """
    Plots a histogram for distribution of stars as a function of RA, dec
    and a 2D histogram heat map based on data from the SUPERBLINK survey.
    """

    df = astropy.io.ascii.read('data_files/aj403664t1_mrt.txt', quotechar="\s")
    coord = SkyCoord(df['RAdeg'], df['DEdeg'], frame='icrs', unit='deg')

    fig1, ax1, ax2, ax3 = plt.subplots(1,3)

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
    """
    Plots ditance distribution based on data form the SUPERBLINK survey.
    """
    
    table = astropy.io.ascii.read('data_files/aj403664t1_mrt.txt')
    distance = np.sort(1.0 / table[(table['plx'] < 9) & (~table['plx'].mask)]['plx'])

    fig1, ax1 = plt.subplots(1,1)
    N_amplitude, bins_amplitude, patches_amplitude = ax1.hist(distance, bins=100, alpha = 0.5) 
    ax1.set_title('Distribution of m dwarfs as a function of distance')
    ax1.set_xlabel('distance in pc')
    ax1.set_ylabel('Number of m dwarfs')
    plt.show()

def plot_effective_kepler_temps():
    """
    Plots the histogram for Teff of stars based on the Kepler Input Catalouge data for the
    541 flaring m dwarfs form the Yang 2017 paper.
    """

    df = pd.read_csv('data_files/eff_temp.csv')
    print(np.mean(df['teff']))
    print(np.std(df['teff']))
    count, bins, ignored = plt.hist(df['teff'], bins=10, alpha = 0.5) 
    plt.xlabel('Effective Temperature in Kelvin')
    plt.ylabel('Number of m dwarfs')
    plt.show()

def plotGenericHistogram(data):
    """
    A generic function to plot a histogram.

    Args:
        data ([Numpy array]): The 1D numpy array for which the histogram is to be plotted
    """

    fig1, ax = plt.subplots(1,1)
    N_amplitude, bins_amplitude, patches_amplitude = ax.hist(data, bins=100, alpha = 0.5) 
    ax.set_ylabel('Number of m dwarfs')
    plt.show()

def plotGeneric2DHistogram(d1, d2):
    """
    A generic function to plot a 2d histogram

    Args:
        d1 (numpy array): 1d numpy array along the x axis
        d2 (numpy array): 1d numpy array along the y axis
    """

    fig1, ax = plt.subplots(1,1)
    ax.hist2d(d1, d2, bins = 20)
    plt.show()

def plotGenricSkyMap(coords):
    """
    A generic function to plot a skymap for the given Sky coord array.

    Args:
        coords (numpy array): A numpy array of skycoord objects
    """
    
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection="mollweide")
    scatter = ax.scatter(-coords.ra.wrap_at(180 * u.deg).radian, coords.dec.wrap_at(180 * u.deg).radian, s=3, vmin=0)
    ax.grid(True)
    ax.set_xticklabels(['10h', '8h', '6h', '4h', '2h', '0h', '22h', '20h', '18h', '16h', '14h'])
    plt.show()

def plotGenricSkyMapWithDistances(coords):
    """
    A genric function to plot a skymap for the given Sky coord array, with distances.

    Args:
        coords (numpy array): A numpy array of skycoord objects
    """

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection="mollweide")
    scatter = ax.scatter(-coords.ra.wrap_at(180 * u.deg).radian, coords.dec.wrap_at(180 * u.deg).radian,c = coords.distance.value, s=3, vmin=0)
    ax.grid(True)
    ax.set_xticklabels(['10h', '8h', '6h', '4h', '2h', '0h', '22h', '20h', '18h', '16h', '14h'])
    plt.colorbar(scatter, label = "Distance in pc")
    plt.show()

def plotFlareFluxAmplitude():
    """
    Plots the distribution for the reletive flux amplitude of all the flares in filtered_flares.csv
    """

    df = pd.read_csv('data_files/filtered_flares.csv')
    amp = df['flux_amp']
    fig1, ax = plt.subplots(1,1)
    N_amplitude, bins_amplitude, patches_amplitude = ax.hist(amp, bins=10, alpha = 0.5) 
    ax.set_xlabel('Relative Flux Amplitude')
    ax.set_ylabel('Number of m dwarfs')
    ax.set_yscale('log')
    plt.show()
    