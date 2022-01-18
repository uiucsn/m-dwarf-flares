from functools import lru_cache
from importlib.resources import open_binary

import numpy as np

import m_dwarf_flare.filters


@lru_cache()
def get_transmission(band, n_grid=200):
    """
    Returns the transmission for a given lsst passband.

    Args:
        band (char): Letter corresponding to an lsst passband.
        n_grid (int): Number of grid points

    Returns:
        [type]: [description]
    """
    with open_binary(m_dwarf_flare.filters, f'LSST_LSST.{band}.dat') as fh:
        lmbd, t = np.genfromtxt(fh, unpack=True)
    step = lmbd.size // n_grid or 1
    return lmbd[::step], t[::step]


@lru_cache()
def get_kepler_transmission(n_grid=200):
    """
    Returns the transmission for the kepler passband.

    Args:
        n_grid (int): Number of grid points

    Returns:
        [type]: [description]
    """
    with open_binary(m_dwarf_flare.filters, 'Kepler_Kepler.K.dat') as fh:
        lmbd, t = np.genfromtxt(fh, unpack=True)
    step = lmbd.size // n_grid or 1
    return lmbd[::step], t[::step]