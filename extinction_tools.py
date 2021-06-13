from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from dustmaps.bayestar import BayestarQuery
from dustmaps.sfd import SFDQuery, fetch
from ch_vars.spatial_distr import MilkyWayDensityJuric2008 as MWDensity

def get_extinction_after_symmetric_interpolation(coordinates):

    bayestar = BayestarQuery(version='bayestar2019')

    gal = coordinates.galactic
    mirrored_coordinates = SkyCoord(l = - 1 * gal.l, b = -1 * gal.b, distance = gal.distance, frame = 'galactic').icrs

    reddening = bayestar(coordinates, mode='median')
    reddening_mirrored = bayestar(mirrored_coordinates, mode='median')

    reddening_interpolated = np.where(np.isnan(reddening), reddening_mirrored, reddening)

    return reddening_interpolated

def get_extinction_after_symmetric_interpolation_with_sfd_factor(coordinates):

    bayestar = BayestarQuery(version='bayestar2019')
    sfd = SFDQuery()

    gal = coordinates.galactic
    mirrored_coordinates = SkyCoord(l = - 1 * gal.l, b = -1 * gal.b, distance = gal.distance, frame = 'galactic').icrs

    reddening_bayes = bayestar(coordinates, mode='median')
    reddening_bayes_mirrored = bayestar(mirrored_coordinates, mode='median')
    reddening_sfd_norm = sfd(coordinates)
    reddening_sfd_flip = sfd(mirrored_coordinates)

    interpolated_reddening_3D = (reddening_sfd_norm / reddening_sfd_flip) * reddening_bayes_mirrored
    reddening_interpolated = np.where(np.isnan(reddening_bayes), interpolated_reddening_3D, reddening_bayes)

    return reddening_interpolated

def plot_extinction_after_symmetric_interpolation():

    ra, dec = np.meshgrid(np.linspace(-180, 180, 361, endpoint=False), np.linspace(-90, 90, 181, endpoint=False))
    coordinates =  SkyCoord(ra = ra * u.degree, dec = dec * u.degree, distance = 2000 * u.pc)
    reddening = get_extinction_after_symmetric_interpolation(coordinates)

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection="mollweide")
    scatter = ax.scatter(-coordinates.ra.wrap_at(180 * u.deg).radian, coordinates.dec.wrap_at(180 * u.deg).radian,c = reddening, cmap = 'gist_stern', s=3, vmin=0)
    ax.grid(True)
    ax.set_xticklabels(['10h', '8h', '6h', '4h', '2h', '0h', '22h', '20h', '18h', '16h', '14h'])
    plt.colorbar(scatter, label = "Extinction")
    plt.show()

def plot_extinction_after_symmetric_interpolation_with_sfd_factor():

    ra, dec = np.meshgrid(np.linspace(-180, 180, 361, endpoint=False), np.linspace(-90, 90, 181, endpoint=False))
    coordinates =  SkyCoord(ra = ra * u.degree, dec = dec * u.degree, distance = 2000 * u.pc)
    reddening = get_extinction_after_symmetric_interpolation_with_sfd_factor(coordinates)

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection="mollweide")
    scatter = ax.scatter(-coordinates.ra.wrap_at(180 * u.deg).radian, coordinates.dec.wrap_at(180 * u.deg).radian,c = reddening, cmap = 'gist_stern', s=3, vmin=0)
    ax.grid(True)
    ax.set_xticklabels(['10h', '8h', '6h', '4h', '2h', '0h', '22h', '20h', '18h', '16h', '14h'])
    plt.colorbar(scatter, label = "Extinction")
    plt.show()