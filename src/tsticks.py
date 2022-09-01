# !/usr/bin/env python
# -*- coding utf-8 -*-
#
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Author: Markus Ritschel
# eMail:  kontakt@markusritschel.de
# Date:   22/10/2020
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
import logging

import pandas as pd
import xarray as xr
from deprecation import deprecated
import tqdm

logger = logging.getLogger(__name__)


def read_tsticks(file):
    """Reads a log file of a T-Stick designed by Leif Riemenschneider. The log file is in tabular form either of the format
        stick sec time t0 t1 t2 t3 t4 t5 t6 t7 (old version from Leif's logger)
    or
        module, hex_id, timestamp, t0, t1, t2, t3, t4, t5, t6, t7 (new version from Arduino logger)

    Data are returned as an xarray.Dataset. Users can select the sticks via the `module` number.
    The new data format comprises a hex id that is unique for each T-stick.
    The hex id is added as a non-dimensional coordinate, and the mapping can be also found in the attributes of the dataset.
    If that hex id is desired as the selector, please apply a `swap_dims` beforehand.
    """
    with open(file) as myfile:
        head = '\n'.join([next(myfile) for x in range(10)])
        if not ',' in head:
            logger.debug("Inferred string format matches the old data format from Leif's logger")
            read_routine = read_tsticks_old
        else:
            logger.debug("Inferred string format matches the new data format")
            read_routine = read_tsticks_new
        return read_routine(file)


@deprecated("Please use the new routine with the new file format.")
def read_tsticks_old(file):
    """Reads a log file of a T-Stick designed by Leif Riemenschneider. The log file is in tabular form of the format
        stick sec time t0 t1 t2 t3 t4 t5 t6 t7

    Data get converted into and returned as an xarray.Dataset.
    """
    col_names = ['stick', 'sec', 'time', 't0', 't1', 't2', 't3', 't4', 't5', 't6', 't7']
    df = pd.read_csv(file, names=col_names, sep=r"\s+", comment='#',
                     engine='python', error_bad_lines=False)

    # cleanup data frame
    df['stick'] = df['stick'].str.strip(':').astype(int)
    df.drop(['sec'], axis=1, inplace=True)

    # convert timme and set index
    df['time'] = pd.to_datetime(df['time'], format='%y%m%d_%H%M_%S', errors='coerce')
    df['time'] = df.time.dt.tz_localize(None)
    df = df.dropna(subset=['time'])
    df.set_index(['time', 'stick'], inplace=True)

    # add column index level
    df.columns = pd.MultiIndex.from_product([['temperature'], df.columns])

    # ensure numeric values
    df = df.apply(pd.to_numeric, errors='coerce')

    df = df.unstack(level=1).reorder_levels([0, 2, 1], axis=1)

    lastcol_idx = df.dropna(subset=[df.columns[-1]]).index
    firstcol_idx = df.dropna(subset=[df.columns[0]]).index

    # ... and fill NANs in all columns backwards such that in each row with `lastcol_idx` are
    # all values of one block
    df.fillna(method='bfill', inplace=True)

    # now limit data frame to those respective rows
    # df = df.loc[firstcol_idx][1:]  # first row has only one single value in column 1
    df = df.loc[lastcol_idx][:-1]  # first row has only one single value in column 1

    # convert back to multi index
    df = df.stack(level=[1, 2])
    df.index = df.index.rename('sensor', level=2)
    sensors = df.index.get_level_values('sensor').unique()
    df.index.set_levels(range(len(sensors)), level=2, inplace=True)  # replace 't0', 't1' etc with 0, 1, etc

    ds = df.to_xarray()

    return ds


def read_tsticks_new(file):
    """Reads a log file of a T-Stick designed by Leif Riemenschneider. The log file is in tabular form of the format
        module, hex_id, timestamp, t0, t1, t2, t3, t4, t5, t6, t7

    Data get converted into and returned as an xarray.Dataset.
    Users can select the sticks via the `module` number. If the hex id is desired, apply a `swap_dims` beforehand.
    The hex id is added as a non-dimensional coordinate, and the mapping can be also found in the attributes of the dataset.
    """
    col_names = ['module', 'hex_id', 'time', 't0', 't1', 't2', 't3', 't4', 't5', 't6', 't7']
    parser = lambda x : pd.to_datetime(x, errors='coerce')
    df = pd.read_csv(file, names=col_names, sep=r",\s+", comment='#',
                     parse_dates=['time'], date_parser=parser,
                     engine='python', error_bad_lines=False)
                     
    df.dropna(subset=['time'], inplace=True)
    
    # convert timme and set index
    df['time'] = pd.to_datetime(df['time']).dt.tz_localize(None)
    df = df.dropna(subset=['time'])
    # df.set_index(['time', 'module'], inplace=True)

    # get unique module <-> hex_id mappings
    module_dict = dict(df.groupby(['module', 'hex_id']).groups.keys())
    df.drop(['hex_id'], axis=1, inplace=True)

    # df['identifier'] = df['module'].astype(str).str.cat(df['hex_id'].astype(str), sep=':')
    df.set_index(['time', 'module'], inplace=True)

    # add column index level
    df.columns = pd.MultiIndex.from_product([['temperature'], df.columns])

    # ensure numeric values
    df = df.apply(pd.to_numeric, errors='coerce')

    # df = df.unstack(level=1).reorder_levels([0, 3, 2, 1], axis=1)
    
    df = df[~df.index.duplicated(keep='first')]
    df = df.unstack(level=1).reorder_levels([0, 2, 1], axis=1)

    lastcol_idx = df.dropna(subset=[df.columns[-1]]).index
    firstcol_idx = df.dropna(subset=[df.columns[0]]).index

    # ... and fill NANs in all columns backwards such that in each row with `lastcol_idx` are
    # all values of one block
    df.fillna(method='bfill', inplace=True)

    # now limit data frame to those respective rows
    # df = df.loc[firstcol_idx][1:]  # first row has only one single value in column 1
    df = df.loc[lastcol_idx][:-1]  # first row has only one single value in column 1

    # convert back to multi index
    df = df.stack(level=[1, 2])
    df.index = df.index.rename('sensor', level=2)
    sensors = df.index.get_level_values('sensor').unique()
    df.index.set_levels(range(len(sensors)), level=2, inplace=True)  # replace 't0', 't1' etc with 0, 1, etc

    ds = df.to_xarray()

    # Add hex_id as non-dimensional coordinate so that users can switch between the two selectors by using `swap_dims` before
    ds.coords['stick_id'] = ('module', list(module_dict.values()))
    ds.attrs['module_dict'] = module_dict

    return ds
