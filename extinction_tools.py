from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from dustmaps.bayestar import BayestarQuery
from dustmaps.sfd import SFDQuery, fetch
from ch_vars.spatial_distr import MilkyWayDensityJuric2008 as MWDensity
from plotting_tools import plotGenricSkyMapWithDistances

def plot_extinction_map_with_flipped_galactic_longitude_and_latitude():
    bayestar = BayestarQuery(version='bayestar2019')

    ra, dec = np.meshgrid(np.linspace(-180, 180, 361, endpoint=False), np.linspace(-90, 90, 181, endpoint=False))
    coordinates =  SkyCoord(ra = ra * u.degree, dec = dec * u.degree, distance = 2000 * u.pc)
    gal = coordinates.galactic
    inv_coordinates = SkyCoord(l = - 1 * gal.l, b = -1 * gal.b, distance = 2000 * u.pc, frame = 'galactic').icrs

    reddening = bayestar(coordinates, mode='median')
    reddening_flipped = bayestar(inv_coordinates, mode='median')

    reddening_interpolated = np.where(np.isnan(reddening), reddening_flipped, reddening)

    fig1 = plt.figure(figsize=(8, 6))
    ax1 = fig1.add_subplot(111, projection="mollweide")
    scatter1 = ax1.scatter(-coordinates.ra.wrap_at(180 * u.deg).radian, coordinates.dec.wrap_at(180 * u.deg).radian,c = reddening, cmap = 'plasma', s=3, vmin=0)
    ax1.grid(True)
    ax1.set_xticklabels(['10h', '8h', '6h', '4h', '2h', '0h', '22h', '20h', '18h', '16h', '14h'])
    plt.colorbar(scatter1, label = "Extinction")
    plt.show()

    fig2 = plt.figure(figsize=(8, 6))
    ax2 = fig2.add_subplot(111, projection="mollweide")
    scatter2 = ax2.scatter(-coordinates.ra.wrap_at(180 * u.deg).radian, coordinates.dec.wrap_at(180 * u.deg).radian,c = reddening_flipped, cmap = 'plasma', s=3, vmin=0)
    ax2.grid(True)
    ax2.set_xticklabels(['10h', '8h', '6h', '4h', '2h', '0h', '22h', '20h', '18h', '16h', '14h'])
    plt.colorbar(scatter2, label = "Extinction")
    plt.show()

    fig3 = plt.figure(figsize=(8, 6))
    ax3 = fig3.add_subplot(111, projection="mollweide")
    scatter3 = ax3.scatter(-coordinates.ra.wrap_at(180 * u.deg).radian, coordinates.dec.wrap_at(180 * u.deg).radian,c = reddening_interpolated, cmap = 'gist_stern', s=3, vmin=0)
    ax3.grid(True)
    ax3.set_xticklabels(['10h', '8h', '6h', '4h', '2h', '0h', '22h', '20h', '18h', '16h', '14h'])
    plt.colorbar(scatter3, label = "Extinction")
    plt.show()

def two_dim_interpolation_test():
    bayestar = BayestarQuery(version='bayestar2019')
    sfd = SFDQuery()
    bayestar = BayestarQuery(version='bayestar2019')

    ra, dec = np.meshgrid(np.linspace(-180, 180, 361, endpoint=False), np.linspace(-90, 90, 181, endpoint=False))
    coordinates =  SkyCoord(ra = ra * u.degree, dec = dec * u.degree, distance = 2000 * u.pc)
    gal = coordinates.galactic
    inv_coordinates = SkyCoord(l = - 1 * gal.l, b = -1 * gal.b, distance = coordinates.distance, frame = 'galactic').icrs

    reddening = bayestar(coordinates, mode='median')
    reddening_flipped = bayestar(inv_coordinates, mode='median')
    reddening_sfd_norm = sfd(coordinates)
    reddening_sfd_flip = sfd(inv_coordinates)

    interpolated_reddening_flip = (reddening_sfd_norm / reddening_sfd_flip) * reddening_flipped
    reddening_interpolated = np.where(np.isnan(reddening), interpolated_reddening_flip, reddening)

    fig2 = plt.figure(figsize=(8, 6))
    ax2 = fig2.add_subplot(111, projection="mollweide")
    scatter2 = ax2.scatter(-coordinates.ra.wrap_at(180 * u.deg).radian, coordinates.dec.wrap_at(180 * u.deg).radian,c = reddening_sfd_norm, cmap = 'gist_stern', s=3, vmin=0)
    ax2.grid(True)
    ax2.set_xticklabels(['10h', '8h', '6h', '4h', '2h', '0h', '22h', '20h', '18h', '16h', '14h'])
    plt.colorbar(scatter2, label = "Extinction")
    plt.show()

    fig3 = plt.figure(figsize=(8, 6))
    ax3 = fig3.add_subplot(111, projection="mollweide")
    scatter3 = ax3.scatter(-coordinates.ra.wrap_at(180 * u.deg).radian, coordinates.dec.wrap_at(180 * u.deg).radian,c = reddening_interpolated, cmap = 'gist_stern', s=3, vmin=0)
    ax3.grid(True)
    ax3.set_xticklabels(['10h', '8h', '6h', '4h', '2h', '0h', '22h', '20h', '18h', '16h', '14h'])
    plt.colorbar(scatter3, label = "Extinction")
    plt.show()

two_dim_interpolation_test()