# -*- coding: utf-8 -*-
'''
Some fundamental functions used for saltharp retrievals
Niels Fuchs 2022-08-15
niels.fuchs@uni-hamburg.de
'''
import numpy as np


def S_brine_SW(T):
    '''
    Brine salinity of sea water solution in [ppt] from Temperature in [°C]
    from Eq.6, Notz and Worster, 2009
    '''
    if not isinstance(T, np.ndarray):
        T=np.array(T)
    ref = np.float64(np.zeros(T.shape)) # avoid S<0 for T>0 values
    return np.max((-0.0170*(T**3) - 0.886*(T**2) - 21.4*(T),ref),axis=0)
    
    
def S_brine_NaCl(T):
    '''
    Brine salinity of NaCl solution in [ppt] from Temperature in [°C]
    Dirks PhD thesis, equ 3.2
    '''
    if not isinstance(T, np.ndarray):
        T=np.array(T)
    ref = np.float64(np.zeros(T.shape)) # avoid S<0 for T>0 values
    return np.max((-0.00362*(T**3) - 0.389*(T**2) - 17.6*(T),ref),axis=0)
    
