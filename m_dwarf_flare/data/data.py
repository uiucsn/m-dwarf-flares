from importlib.resources import open_binary
from functools import lru_cache

import astropy.io.ascii
import pandas as pd

import m_dwarf_flare.data


@lru_cache(maxsize=1)
def dist():
    with open_binary(m_dwarf_flare.data, 'dist.csv') as fh:
        return astropy.io.ascii.read(fh, format='csv')


@lru_cache(maxsize=1)
def mag():
    with open_binary(m_dwarf_flare.data, 'mag.csv') as fh:
        return astropy.io.ascii.read(fh, format='csv')


@lru_cache(maxsize=1)
def get_LSST_extinction():
    with open_binary(m_dwarf_flare.data, 'LSST_rv_3_1.csv') as fh:
        df = pd.read_csv(fh)
    rv = pd.Series(df['value'].values,index=df['lsst_passband']).to_dict()
    return rv


@lru_cache(maxsize=1)
def yang_table3():
    with open_binary(m_dwarf_flare.data, 'apjaa8ea2t3_mrt.txt') as fh:
        return astropy.io.ascii.read(fh, quotechar="\s")


@lru_cache(maxsize=1)
def lepine_gaidos_table1():
    with open_binary(m_dwarf_flare.data, 'aj403664t1_mrt.txt') as fh:
        return astropy.io.ascii.read(fh, quotechar="\s")


@lru_cache(maxsize=1)
def filtered_flares():
    with open_binary(m_dwarf_flare.data, 'filtered_flares.csv') as fh:
        return pd.read_csv(fh)