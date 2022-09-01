# -*- coding: utf-8 -*-
'''
Fixed parameters and variable definitions used for saltharp retrievals
Niels Fuchs 2022-08-15
niels.fuchs@uni-hamburg.de
'''

# Liquid fraction retrieval parameter

rise_percent = 3. # rise the initial R_ref (Z_0 in equations) by this percentage. Mentionned by Dirk before, confirmed by Niels to be the best fit

# gamma term values valid for old harp models with 2cm vertical spacing. Retrieved in Fuchs (2017)

a_gamma = 0.101 
b_gamma = -0.101
c_gamma = 0.102 # note the minus before c in the gamma term in equ. (28) in Fuchs (2017)

# Density of solid ice

rho_s = 920.

# Metadata definitions

attrs_definition_dict={
    'S_bulk':{'units':'ppt', 'long_name':'Bulk salinity'},
    'S_brine':{'units':'ppt', 'long_name':'Brine salinity'},
    'T_deg':{'units':'degree_Celsius', 'long_name':'Temperature [°C]', 'standard_name':'sea_ice_temperature'},
    'phi_l':{'units':'1', 'long_name':'liquid volume fraction'},
    'phi_s':{'units':'1', 'long_name':'solid volume fraction'},
    'phi_l_mass':{'units':'1', 'long_name':'liquid mass fraction'},
    'R_ref':{'units':'Ohm', 'long_name':'Reference resistance at freezing, equals Z_0 in harp equation'},
    'T_box':{'units':'degree_Celsius', 'long_name':'Temperature in the data acquisition unit [°C]'},
    'R_2':{'units':'Ohm', 'long_name':'Measured wire pair resistance at 2kHz excitation AC'},
    'R_16':{'units':'Ohm', 'long_name':'Measured wire pair resistance at 16kHz excitation AC'},
    'wire_pair':{'long_name':'horizontal wire pair'},
    'module':{'long_name':'slot number of the salt harp module on the acquisition unit'},
    'date':{'long_name':'teraterm generated timestamp at the serial receive of the measurement'},
    'nominal_date':{'long_name':'unified timestamp of the profile measurement'},
    'general_attrs':{'description':"Salinity harp retrievals. Based on Notz et al. (2005), Fuchs (2017)",
                     'script_version':"1.0, (N. Fuchs 2022-08-15)"}
}