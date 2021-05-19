import astropy.coordinates as coord
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
from astropy.table import Table
from astropy.coordinates import SkyCoord, CylindricalRepresentation
from astropy import units as u
import astropy
import numpy as np
from pyvo.dal import TAPService

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



    # distance = []

    # for i in range(10):
    #     search_radius_arcsec = 10
    #     string_query = f'''SELECT distance(ra, dec, {ra_list[i]}, {dec_list[i]}) as d, r_med_geo, r_lo_geo, r_hi_geo, r_med_photogeo,r_lo_photogeo, r_hi_photogeo, phot_g_mean_mag FROM gedr3dist.main JOIN gaia.edr3lite USING (source_id) WHERE distance(ra, dec, {ra_list[i]}, {dec_list[i]}) < {search_radius_arcsec / 3600.0}'''
        
    #     tap = TAPService('https://dc.zah.uni-heidelberg.de/tap')
    #     response = tap.search(string_query)
    #     r = response.to_table()

    #     if len(r) > 1:
    #         index = np.where(r['r_med_geo'] == np.amin(r['r_med_geo']))
    #         print(r[index])
    #         distance.append(r[index]['r_med_geo'])
    #         continue
    #     elif len(r) == 0:
    #         print('no data')
    #         continue
    #     else:
    #         print(r[0])
    #         distance.append(r[0]['r_med_geo'])

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

plotLBDistribution()
#plotSkyMap()