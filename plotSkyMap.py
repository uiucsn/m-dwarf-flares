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
    ra_list = []
    dec_list = []

    for i in range(len(t['RGAL'])):
        c = CylindricalRepresentation(rho = t['RGAL'][i], phi = t['THETA'][i] * u.rad, z = t['ZGAL'][i])
        s = SkyCoord(c)
        ra_list.append(s.ra.to(u.degree))
        dec_list.append(s.dec.to(u.degree))


    print(ra_list)
    print(dec_list)
    ra = coord.Angle(ra_list)
    ra = ra.wrap_at(180*u.degree)
    dec = coord.Angle(dec_list)

    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111, projection="mollweide")
    ax.scatter(ra.radian, dec.radian)
    ax.grid(True)
    plt.show()

def plotSkyMap():

    ra_list = []
    dec_list = []

    df = astropy.io.ascii.read('aj403664t1_mrt.txt', quotechar="\s")
    ra_list = df['RAdeg']
    dec_list = df['DEdeg']
    print(df)

    ra = coord.Angle(ra_list)
    ra = ra.wrap_at(180*u.degree)
    dec = coord.Angle(dec_list)

    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111, projection="mollweide")
    sm = ax.scatter(ra.radian, dec.radian)
    ax.grid(True)
    plt.colorbar(sm)
    plt.show()

def plotBDistribution():

    ra_list = []
    dec_list = []

    df = astropy.io.ascii.read('aj403664t1_mrt.txt', quotechar="\s")
    ra_list = df['RAdeg']
    dec_list = df['DEdeg']

    print(df)



    ra = coord.Angle(ra_list)
    ra = ra.wrap_at(180*u.degree)
    dec = coord.Angle(dec_list)
    l = []
    b = []

    for i in range(len(ra)):
        print(i)
        c = SkyCoord(ra[i], dec[i], frame='icrs', unit='deg')
        g = c.transform_to('galactic')
        l.append(g.l.value)
        b.append(g.b.value)


    fig1, ax1 = plt.subplots(1,1)

    N_amplitude, bins_amplitude, patches_amplitude = ax1.hist(b, bins=100, alpha = 0.5) 
    ax1.set_title('Distribution of m dwarfs as a function of b')
    ax1.set_xlabel('b in degrees')
    ax1.set_ylabel('Number of m dwarfs')
    plt.show()

def plotLDistribution():

    ra_list = []
    dec_list = []

    df = astropy.io.ascii.read('aj403664t1_mrt.txt', quotechar="\s")
    ra_list = df['RAdeg']
    dec_list = df['DEdeg']

    print(df)



    ra = coord.Angle(ra_list)
    ra = ra.wrap_at(180*u.degree)
    dec = coord.Angle(dec_list)
    l = []
    b = []

    for i in range(len(ra)):
        print(i)
        c = SkyCoord(ra[i], dec[i], frame='icrs', unit='deg')
        g = c.transform_to('galactic')
        l.append(g.l.value)
        b.append(g.b.value)


    fig1, ax1 = plt.subplots(1,1)

    N_amplitude, bins_amplitude, patches_amplitude = ax1.hist(l, bins=100, alpha = 0.5) 
    ax1.set_title('Distribution of m dwarfs as a function of l')
    ax1.set_xlabel('l in degrees')
    ax1.set_ylabel('Number of m dwarfs')
    plt.show()

def plotLBDistribution():

    ra_list = []
    dec_list = []

    df = astropy.io.ascii.read('aj403664t1_mrt.txt', quotechar="\s")
    ra_list = df['RAdeg']
    dec_list = df['DEdeg']

    print(df)

    ra = coord.Angle(ra_list)
    ra = ra.wrap_at(180*u.degree)
    dec = coord.Angle(dec_list)
    l = []
    b = []

    for i in range(len(ra)):
        print(i)
        c = SkyCoord(ra[i], dec[i], frame='icrs', unit='deg')
        g = c.transform_to('galactic')
        l.append(g.l.value)
        b.append(g.b.value)

    l = np.array(l)
    b = np.array(b)

    plt.hist2d(l, b, bins = 20)
    plt.title("Distribution of m - dwarfs in the galaxy")
    plt.xlabel("Galactic Longitude in degrees")
    plt.ylabel("Galactic Latitude in degrees")
    plt.colorbar(label = "Number of m dwarfs")
    plt.show()

def plotSkyMapWithDistances():

    ra_list = []
    dec_list = []

    df = astropy.io.ascii.read('aj403664t1_mrt.txt', quotechar="\s")
    ra_list = df['RAdeg']
    dec_list = df['DEdeg']
    dist = distance = 1.0 / df['plx']

    ra = coord.Angle(ra_list)
    ra = ra.wrap_at(180*u.degree)
    dec = coord.Angle(dec_list)
    
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111, projection="mollweide")
    sm = ax.scatter(ra.radian, dec.radian, c = dist)
    ax.grid(True)
    plt.colorbar(sm, label = "Distance in pc")
    plt.show()

def plotLDistCorrelation():

    ra_list = []
    dec_list = []

    df = astropy.io.ascii.read('aj403664t1_mrt.txt', quotechar="\s")
    ra_list = df['RAdeg']
    dec_list = df['DEdeg']
    dist = distance = 1.0 / df['plx']

    ra = coord.Angle(ra_list)
    ra = ra.wrap_at(180*u.degree)
    dec = coord.Angle(dec_list)

    l = []
    b = []

    for i in range(len(ra)):
        c = SkyCoord(ra[i], dec[i], frame='icrs', unit='deg')
        g = c.transform_to('galactic')
        l.append(g.l.value)
        b.append(g.b.value)
    l = np.array(l)
    b = np.array(b)
    
    plt.scatter(b,dist)
    plt.show()

def plotDensityDistribution():

    table = astropy.io.ascii.read('aj403664t1_mrt.txt')
    distance = np.sort(1.0 / table[(table['plx'] < 9) & (~table['plx'].mask)]['plx'])
    plt.plot(distance, np.arange(1, distance.size+1), label='data')
    plt.plot(np.linspace(2, 20, 100), np.linspace(2, 20, 100)**3 / 20**3 * 1200,label=r'$1200 \times (r / 20 pc)^3$')
    plt.xlabel('r, pc')
    plt.ylabel('N')
    plt.legend()
    plt.xscale('log')
    plt.yscale('log')
    plt.show()

def plotDistanceDistribution():

    table = astropy.io.ascii.read('aj403664t1_mrt.txt')
    distance = np.sort(1.0 / table[(table['plx'] < 9) & (~table['plx'].mask)]['plx'])

    fig1, ax1 = plt.subplots(1,1)
    N_amplitude, bins_amplitude, patches_amplitude = ax1.hist(distance, bins=100, alpha = 0.5) 
    ax1.set_title('Distribution of m dwarfs as a function of distance')
    ax1.set_xlabel('distance in pc')
    ax1.set_ylabel('Number of m dwarfs')
    plt.show()

def plotDistanceDistributionWithBellCurve():

    table = astropy.io.ascii.read('aj403664t1_mrt.txt')
    distance = np.sort(1.0 / table[(table['plx'] < 9) & (~table['plx'].mask)]['plx'])

    mu = np.mean(distance)
    sig = np.std(distance)
    d = np.random.normal(mu, sig, 1000)

    count, bins, ignored = plt.hist(distance, bins=100, alpha = 0.5, density=True) 
    plt.xlabel('distance in pc')
    plt.ylabel('Number of m dwarfs')
    plt.plot(bins, 1/(sig * np.sqrt(2 * np.pi)) *
               np.exp( - (bins - mu)**2 / (2 * sig**2)),
         linewidth=2, color='r')
    plt.show()

def plotBDistributionWithBellCurve():

    ra_list = []
    dec_list = []

    df = astropy.io.ascii.read('aj403664t1_mrt.txt', quotechar="\s")
    ra_list = df['RAdeg']
    dec_list = df['DEdeg']

    print(df)



    ra = coord.Angle(ra_list)
    ra = ra.wrap_at(180*u.degree)
    dec = coord.Angle(dec_list)
    l = []
    b = []

    for i in range(len(ra)):
        c = SkyCoord(ra[i], dec[i], frame='icrs', unit='deg')
        g = c.transform_to('galactic')
        l.append(g.l.value)
        b.append(g.b.value)

    mu = np.mean(b)
    sig = np.std(b)
    k = np.random.normal(mu, sig, 1000)

    count, bins, ignored = plt.hist(b, bins=100, alpha = 0.5, density=True) 
    plt.xlabel('distance in pc')
    plt.ylabel('Number of m dwarfs')
    plt.plot(bins, 1/(sig * np.sqrt(2 * np.pi)) *
               np.exp( - (bins - mu)**2 / (2 * sig**2)),
         linewidth=2, color='r')
    plt.show()

def test_spherical_distribution():
    count = 1000
    R = 1000

    ra = []
    dec = []
    dist = []


    while len(ra) != 1000:
        x = (np.random.uniform(-1 * R,R) * u.pc)
        y = (np.random.uniform(-1 * R,R) * u.pc)
        z = (np.random.uniform(-1 * R,R) * u.pc)

        if (x ** 2 + y ** 2 + z ** 2) ** (0.5) < R * u.pc:
            r, lat, lon = coord.cartesian_to_spherical(x,y,z)
            ra.append(lon)
            dec.append(lat)
            dist.append(r.value)
    
    ra = coord.Angle(ra)
    ra = ra.wrap_at(180*u.degree)
    dec = coord.Angle(dec)

    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111, projection="mollweide")
    sm = ax.scatter(ra.radian, dec.radian, c = dist)
    ax.grid(True)
    plt.colorbar(sm)
    plt.show()

def test_dist_distribution():
    count = 1000
    R = 1000

    ra = []
    dec = []
    dist = []


    while len(ra) != 1000:
        x = (np.random.uniform(-1 * R,R) * u.pc)
        y = (np.random.uniform(-1 * R,R) * u.pc)
        z = (np.random.uniform(-1 * R,R) * u.pc)

        if (x ** 2 + y ** 2 + z ** 2) ** (0.5) < R * u.pc:
            r, lat, lon = coord.cartesian_to_spherical(x,y,z)
            ra.append(lon)
            dec.append(lat)
            dist.append(r.value)
    
    count, bins, ignored = plt.hist(dist, bins=10, alpha = 0.5) 
    plt.xlabel('distance in pc')
    plt.ylabel('Number of m dwarfs')
    plt.show()

def test_dec_distribution():
    count = 1000
    R = 1000

    ra = []
    dec = []
    dist = []


    while len(ra) != 1000:
        x = (np.random.uniform(-1 * R,R) * u.pc)
        y = (np.random.uniform(-1 * R,R) * u.pc)
        z = (np.random.uniform(-1 * R,R) * u.pc)

        if (x ** 2 + y ** 2 + z ** 2) ** (0.5) < R * u.pc:
            r, lat, lon = coord.cartesian_to_spherical(x,y,z)
            ra.append(lon.degree)
            dec.append(lat.degree)
            dist.append(r.value)
    
    count, bins, ignored = plt.hist(dec, bins=100, alpha = 0.5) 
    plt.xlabel('Declination')
    plt.ylabel('Number of m dwarfs')
    plt.show()

def test_ra_distribution():
    count = 1000
    R = 1000

    ra = []
    dec = []
    dist = []


    while len(ra) != 1000:
        x = (np.random.uniform(-1 * R,R) * u.pc)
        y = (np.random.uniform(-1 * R,R) * u.pc)
        z = (np.random.uniform(-1 * R,R) * u.pc)

        if (x ** 2 + y ** 2 + z ** 2) ** (0.5) < R * u.pc:
            r, lat, lon = coord.cartesian_to_spherical(x,y,z)
            ra.append(lon.degree)
            dec.append(lat.degree)
            dist.append(r.value)
    
    count, bins, ignored = plt.hist(ra, bins=10, alpha = 0.5) 
    plt.xlabel('RA')
    plt.ylabel('Number of m dwarfs')
    plt.show()

def save_effective_kepler_temps():

    mag = pd.read_csv('mag.csv')
    kepler_cataloug = pd.read_csv('kepler_kic_v10.csv.gz')
    
    df = pd.DataFrame()
    k = kepler_cataloug.loc[kepler_cataloug['kic_kepler_id'].isin(mag['KIC ID'].astype(int))]
    
    df['KIC ID'] = k['kic_kepler_id']
    df['teff'] = k['kic_teff']

    df.to_csv('eff_temp.csv')

def plot_effective_kepler_temps():

    df = pd.read_csv('eff_temp.csv')
    print(np.mean(df['teff']))
    print(np.std(df['teff']))
    count, bins, ignored = plt.hist(df['teff'], bins=10, alpha = 0.5) 
    plt.xlabel('Effective Temperature in Kelvin')
    plt.ylabel('Number of m dwarfs')
    plt.show()

def remove_incomplete_entries_from_flare_data():
    dist = pd.read_csv('dist_new.csv')
    flare_data = astropy.io.ascii.read('flare_data/apjaa8ea2t3_mrt.txt', quotechar="\s")
    df = flare_data.to_pandas()

    filtered = df.loc[df['KIC'].isin(dist['KIC ID'].astype(int))]
    df.to_csv('filtered_flares.csv')

plot_effective_kepler_temps()
remove_incomplete_entries_from_flare_data()

    

    






