# -*- coding: utf-8 -*-
'''
Routine for reading saltharp data and obtaining retrievals
Based on Notz et al. (2005) and Master thesis Fuchs (2022)
Compiled by Niels Fuchs 2022-08-15
niels.fuchs@uni-hamburg.de

ToDo: check script to make it work with several modules
'''
import datetime as dt
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
import xarray as xr
from saltharp_parameters import *
from saltharp_functions import *
from matplotlib.dates import date2num

def read_teraterm_logfile(file_path, max_timedelta_reading=30):
    '''
    Read teraterm logfile format, skip bailout lines, SD errors and lines not starting with a timestamp "["
    '''
    
    exclude = [i for i, line in enumerate(open(file_path)) if (not line.startswith('[') or line.split(']')[1].startswith(' F') or line.split(']')[1].startswith(' ba'))] # exclude all lines without timestamp, SD errors and bailouts, not using bash commands here to make the script also working under Windows
    
    date_parser = lambda x, y: dt.datetime.strptime(f"{x} {y}",  '[%Y-%m-%d %H:%M:%S.%f]') # parser for strange teraterm time log format
    
    df=pd.read_csv(file_path,delimiter=' ', error_bad_lines=False, usecols=[0,1,2,4,6,8,9], names=['dates','times', 'device', 'R_2', 'R_16', 'T_deg', 'T_box'], parse_dates={'date': ['dates', 'times']}, date_parser=date_parser, skiprows=exclude, 
    dtype={'R_2':'f4', 'R_16':'f4', 'T_deg':'f4', 'T_box':'f4'})
    
    # split device column into module and wire pairs
    
    df['module']=[int(device[0]) for device in df.device] # module line
    df['wire_pair']=[int(device[2]) for device in df.device] # sensors line
    
    df.drop('device', axis='columns', inplace=True)
    
    # Add new timeaxis with one timestamp per profile, check function definition for further information
    
    df['nominal_date'] = get_nominal_date(df, max_timedelta_reading)
    
    # set multiindex
    
    df.set_index(['nominal_date', 'wire_pair', 'module'], inplace=True)
    
    df.sort_index(inplace=True)
    
    # convert to xarray
    
    ds=df.to_xarray()
    
    # add some meta information to the file

    ds.T_deg.attrs=attrs_definition_dict['T_deg']
    ds.T_box.attrs=attrs_definition_dict['T_box']
    ds.R_2.attrs=attrs_definition_dict['R_2']
    ds.R_16.attrs=attrs_definition_dict['R_16']
    ds.wire_pair.attrs=attrs_definition_dict['wire_pair']
    ds.module.attrs=attrs_definition_dict['module']
    ds.date.attrs=attrs_definition_dict['date']
    ds.nominal_date.attrs=attrs_definition_dict['nominal_date']
    
    ds.attrs=attrs_definition_dict['general_attrs']
    
    return ds
    
    
def get_nominal_date(df, max_timedelta_reading):
    '''
    Each wire pair reading has its own timestamp. However, to facilitate evaluation it is easier when all readings of a profile are on the same time. 
    There would be two ways to solve this: interpolation or simple reassignement. Interpolation has the problem of rendering at least 3 values 
    invalid in the case of NaNs/ false outliers, which occur repeatedly in the saltharp data. Hence we chose reassignement ot timestamps, 
    which means all profile measurements are assigned to the time of the reading of the first sensor. We call that the "nominal date". 
    If sensor 0 was skipped by any reason, we use sensor 1 and so on. There is also a new timestamp for every new harp module reading 
    or when more time elapsed than specified by max_timedelta_reading (default=30s).
    '''
    next_profile_bool = (df.wire_pair[1:].values-df.wire_pair[:-1].values) <= 0 # True whenever wire pair number smaller or equal to the one before
    
    next_module_bool = np.not_equal(df.module[1:].values,df.module[:-1].values) # True whenever module changes
    
    timedelta_too_large_bool = (df.date[1:].values-df.date[:-1]) > np.timedelta64(max_timedelta_reading,'s') # True when too much time passes before reading
    
    new_profile_bool = np.logical_or(
                                    np.logical_or(next_profile_bool, 
                                                next_module_bool), 
                                    timedelta_too_large_bool)
                                    
    new_profile_bool = np.concatenate([[True],new_profile_bool])    # add value for the first measurement
    
    new_profile_left_idx = np.where(new_profile_bool)[0] # indexes of the dates to use
    
    new_profile_idx = new_profile_left_idx[np.cumsum(new_profile_bool)-1] # compile array of full length with indexes of the dates to use, subtract 1 as cumsum starts with 1, but first index is 0
    
    nominal_date = df.iloc[new_profile_idx]['date'].values
    
    return nominal_date

def depth_axis(ds, flip_ud=False, vertical_spacing = 0.02, upper_wire_pair_depth = 0.01):
    '''
    create sensor depth axis (without considering surface melt)
    requires:
    - vertical spacing between wire pairs in [m]
    - depth of the uppermost wire pair in [m]
    - boolean that is True, when harp was flipped with wire pair 7 at top. False when normal aligned with 0 at top and 7 at bottom. 
    '''
    if flip_ud:
        ds.coords['depth'] = ('wire_pair',[upper_wire_pair_depth+i*vertical_spacing for i in range(7,-1,-1)])
    else:
        ds.coords['depth'] = ('wire_pair',[upper_wire_pair_depth+i*vertical_spacing for i in range(8)])
        
    return ds
    
def adjust_gamma_coefficients(a, b, c, T_ref, S_ref, R_ref):
    '''
    The gamma term basically calibrates the temperature and salinity dependency of the conductivity measurement of the salinity harps. 
    There are no up-to-date measurements. The last ones where measured in Fuchs (2017) for previous harp models. Hence, we need to 
    assume, that the relationships remain the same. However, we can at least adjust the zeros value, which describes the conductivity 
    at the time, when R_ref (in the euations Z_0) is determined. 
    This function is written to adapt a, b and c from Fuchs (2017), equation (28) in Section 7.2 to the conductivity of the used harp model
    '''
    cond_ref = 1./R_ref
    
    # the equation is: gamma(T,S) = a + b * exp(-c*S) - 0.00167 * (0.5-T)
    # there are two solutions: 
    #           1. gamma(T,S=0) = 0, we assume that the conductivity is zero when S=0.
    #              We furthermore assume that del_gamma(T)/del_T << del_gamma(S)/del_S, and the equation is therefore valid for all T~T_freezing and thus: 
    #              gamma(T_ref, S=0) = 0, independent of the exact T_ref. 
    #           2. gamma(T_ref, S_ref) = cond_ref
    # we keep the temperature and exponential term constant, hence two solutions are enough to recalibrate coefficients a and b
    ######
    # this gives:
    #           1. 0 = a + b - temperature_term -> a = -b + temperature_term(T=T_ref)
    #           2. cond_ref = -b + temperature_term + b * exponential_term(S=S_ref) - temperature_term
    #           -> b = cond_ref/(exponential_term(S=S_ref)-1)
    #####
    
    temperature_term = 0.00167*(0.5-T_ref)
    
    exponential_term = np.exp(-c_gamma * S_ref)
    
    b_new = cond_ref/(exponential_term-1.)
    a_new = -b_new + temperature_term
    
    return a_new, b_new
    
def replace_T_with_T_Sticks(ds, module_dict, stick_file_path):
    '''
    Replace Temperature measurements of all harp modules given in module dict with T_stick data. 
    
    Module_dict:
        keys: harp modules
        vals: [T_stick number, Stick_flipud_bool]
                                with: Stick_flipud_bool: True, when Tstick vertical sensor counting is different to harp counting (example: wire_pair:0 equals Tstick sensor:7)
        stick_file_path: file path to T_Stick log file
    
    Function uses saltharp "date" array for interpolation instead of "nominal_date" to be more precise.
    '''
    import tsticks
    
    ds_sticks = tsticks.read_tsticks_new(stick_file_path)
    
    for key in module_dict:
        for wire_pair in range(8):
            harp_locator = dict(module=key, wire_pair=wire_pair)
            stick_locator = dict(module=module_dict[key][0], sensor=(wire_pair-7*module_dict[key][1])*(1-2*module_dict[key][1])) # flip the stick sensor numbering if necessary
        
            ds['T_deg'].loc[harp_locator] = interp1d(date2num(ds_sticks['time']), 
                                                     ds_sticks['temperature'].loc[stick_locator], 
                                                     kind='linear', bounds_error=False, fill_value=np.nan)(date2num(ds['date'].loc[harp_locator]))
    
    return ds

    
def obtain_retrievals(ds, start, end, salt='SW', gamma_method='gamma_warm', resistance_value='R_2', sea_salinity=None, correct_gamma_coefficients=True):
    '''
    Calculates brine and bulk salinity and liquid volume fraction from saltharp measurement.
    Inputs:
    Start, end: time range from which reference resistance (impedance) R_0 (Z_0) should be chosen.
    salt: NaCl or SW for sea water. Defines Brine salinity function
    gamma_method: no_gamma if gamma term in equation (2) Notz et al 2005 should be supressed, 
                  gamma if gamma term is considered with formulas fitted to saltharp models 2016 used in Fuchs (2017), 
                  gamma_warm if gamma term is only to be considered when the desalinated ice becomes warmer than the ocean freezing temperature, as motivated in Fuchs (2017)
    resistance_value: resistance measurements to use for retrieval, depends on excitation frequency. Theoretically higher values should perfom better. But somehow 2kHz performed best in Fuchs (2017).
    sea_salinity: Salinity of the sea water in [g kg] for gamma_0 retrieval; if not specified, sea_salinity is determined from brine salinity during the start, end time range 
    
    default settings: salt='SW', gamma_method='gamma_warm', resistance_value='R2'
    '''
    
    # calculate brine salinity
    
    if salt == 'NaCl':
        S_brine = S_brine_NaCl
    elif salt == 'SW':
        S_brine = S_brine_SW

    ds['S_brine'] = (ds['T_deg'].dims, S_brine(ds['T_deg']))

    # initialize data arrays
    
    ds['R_ref'] = xr.full_like(ds['T_deg'], fill_value=np.nan)
    ds['phi_l'] = xr.full_like(ds['T_deg'], fill_value=np.nan)
    ds['phi_s'] = xr.full_like(ds['T_deg'], fill_value=np.nan)
    ds['phi_l_mass'] = xr.full_like(ds['T_deg'], fill_value=np.nan)
    ds['S_bulk'] = xr.full_like(ds['T_deg'], fill_value=np.nan)
    
    if gamma_method in ['gamma', 'gamma_warm']:
        ds['gamma'] = xr.full_like(ds['T_deg'], fill_value=np.nan)
        ds['gamma_0'] = xr.full_like(ds['T_deg'], fill_value=np.nan)
        
        # add some metadata
        ds.gamma_0.attrs={'units':'1/Ohm', 'long_name':'Conductivity at freezing'}
        ds.gamma.attrs={'units':'1/Ohm', 'long_name':'Conductivity of the liquid brine'}
    
    # obtain retrievals
    
    for module in ds.module:
        for wire_pair in ds.wire_pair:
            
            device_selector = dict(wire_pair=wire_pair, module=module)
            time_selector = dict(nominal_date=slice(start,end))
            
            ds['R_ref'].loc[device_selector] = np.nanmean(ds[resistance_value].loc[{**device_selector, **time_selector}])*\
                                                (1+rise_percent/100.)
                        
            if gamma_method in ['gamma', 'gamma_warm']:
                
                if not sea_salinity:
                    sea_salinity = ds['S_brine'].loc[{**device_selector, **time_selector}]
                
                if correct_gamma_coefficients:
                    T_ref = np.nanmean(ds['T_deg'].loc[{**device_selector, **time_selector}])
                    S_ref = np.nanmean(sea_salinity)
                    R_ref = np.nanmean(ds['R_ref'].loc[device_selector])
                    
                    a_gamma, b_gamma = adjust_gamma_coefficients(a_gamma, b_gamma, c_gamma, T_ref, S_ref, R_ref) # look into the function for more information
                
                ds['gamma_0'].loc[device_selector] = np.nanmean(a_gamma + b_gamma * np.exp(-c_gamma*\
                    sea_salinity)-0.00167*\
                    (0.5-ds['T_deg'].loc[{**device_selector, **time_selector}]))
                
                ds['gamma'].loc[device_selector] = a_gamma + b_gamma * np.exp(-c_gamma*\
                    ds['S_brine'].loc[device_selector])-0.00167*(0.5-ds['T_deg'].loc[device_selector])
                                    
                if gamma_method == 'gamma':
                    ds['phi_l'].loc[device_selector] = (ds['R_ref'].loc[device_selector]*ds['gamma_0'].loc[device_selector]) / \
                                                        (ds[resistance_value].loc[device_selector]*ds['gamma'].loc[device_selector])
                elif gamma_method == 'gamma_warm':
                    
                    # no gamma correction when T below freezing:
                    
                    cold_times = ds['T_deg'].loc[device_selector] <= np.nanmean(ds['T_deg'].loc[{**device_selector, **time_selector}])
                    
                    ds['phi_l'].loc[device_selector] = xr.where(cold_times, ds['R_ref'].loc[device_selector]/ds[resistance_value].loc[device_selector], ds['phi_l'].loc[device_selector])
                    
                    # gamma correction when desalinated ice becomes warmer than initial freezing:
                    
                    warm_times = ds['T_deg'].loc[device_selector] > np.nanmean(ds['T_deg'].loc[{**device_selector, **time_selector}]) 
                    
                    ds['phi_l'].loc[device_selector] = xr.where(warm_times, 
                                                                (ds['R_ref'].loc[device_selector]*ds['gamma_0'].loc[device_selector]) / \
                                                                        (ds[resistance_value].loc[device_selector]*ds['gamma'].loc[device_selector]), 
                                                                ds['phi_l'].loc[device_selector])
                                                                
                    
            
            elif gamma_method == 'no_gamma':
                ds['phi_l'].loc[device_selector] = ds['R_ref'].loc[device_selector]/ds[resistance_value].loc[device_selector]
        
            ds['phi_s'].loc[device_selector] = 1. - ds['phi_l'].loc[device_selector]  # assuming there is no gas.
            
            ds['phi_l_mass'].loc[device_selector] = liquid_volume_fraction_2_liquid_mass_fraction(ds['phi_l'].loc[device_selector], 
                                                                                   ds['T_deg'].loc[device_selector], 
                                                                                   ds['S_brine'].loc[device_selector],
                                                                                   rho_s)

        
            ds['S_bulk'].loc[device_selector] = ds['phi_l_mass'].loc[device_selector]*ds['S_brine'].loc[device_selector]
    
    
    # update metadata
    
    ds.S_bulk.attrs=attrs_definition_dict['S_bulk']
    ds.S_brine.attrs=attrs_definition_dict['S_brine']
    ds.phi_l.attrs=attrs_definition_dict['phi_l']
    ds.phi_s.attrs=attrs_definition_dict['phi_s']
    ds.phi_l_mass.attrs=attrs_definition_dict['phi_l_mass']
    ds.R_ref.attrs=attrs_definition_dict['R_ref']

    
    # last polishing: change dimension to depth
    
    ds = ds.swap_dims({'wire_pair':'depth'})    
    
    return ds

def liquid_volume_fraction_2_liquid_mass_fraction(liquid_volume_fraction, T, S_br, rho_s):
    # equ 3.20 Dirk PhD, Seite 44 bzw. 60 (pdf) modified
    return 1.-(1./(1.+((1./(1.-liquid_volume_fraction))-1.)*rho_l(T, S_br)/rho_s))

def rho_l(T, S):
    ### Chiareli, could be improved ;)
    return 1028.

def standard_conversion(file_path, start, end, salt='SW', gamma_method='no_gamma', sea_salinity=35., resistance_value='R_2', flip_ud=False, \
                        correct_gamma_coefficients=True, vertical_spacing=0.02, upper_wire_pair_depth=0.0, T_Stick_module_dict=None, stick_file_path=None):
    '''
    Standard conversion pipeline for saltharp data
    
    file_path (str):                                file path of saltharp log file
    start (date string 'YYYY-MM-DD hh:mm:ss.sss'):                     
    start (date string 'YYYY-MM-DD hh:mm:ss.sss'):
    salt
    gamma_method
    sea_salinity
    resistance_value
    flip_ud
    correct_gamma_coefficients
    vertical_spacing
    upper_wire_pair_depth
    T_Stick_module_dict                    
    '''       
    # read data to xarray
    ds = saltharp.read_teraterm_logfile(file_path)
    
    # add ice depth column
    ds = saltharp.depth_axis(ds, flip_ud, vertical_spacing, upper_wire_pair_depth)

    # replace temperature data with Tsticks if wanted
    if T_Stick_module_dict and stick_file_path:
        ds = saltharp.replace_T_with_T_Sticks(ds, T_Stick_module_dict, stick_file_path)

    # calculate liquid fraction and salinity
    ds = saltharp.obtain_retrievals(ds, start, end, salt, resistance_value, gamma_method, sea_salinity, correct_gamma_coefficients)
    
    # write NetCDF file
    ds.to_netcdf(file_path.rsplit('.',1)[0]+'_retrievals.nc')

if __name__ == '__main__':

    file_path='/Users/nielsfuchs/Library/CloudStorage/OneDrive-Personal/UniHH/Labor/Lab_projects/Salzharfe_EllensExperimente/data/SaltharpLog_GroßerTank_Versuch2.log'
    
    ds = read_teraterm_logfile(file_path)
    
    ds = depth_axis(ds, flip_ud, vertical_spacing, upper_wire_pair_depth)

    ds = obtain_retrievals(ds, '2022-07-28 16:07:00', '2022-07-28 16:25:00', salt='NaCl', gamma_method='no_gamma', sea_salinity=35.)
    
    ds.to_netcdf(file_path.rsplit('.',1)[0]+'_retrievals.nc')

