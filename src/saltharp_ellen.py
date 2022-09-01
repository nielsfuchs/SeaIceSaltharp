# -*- coding: utf-8 -*-
'''
Obtain saltharp retrievals for Ellens Master project
Niels Fuchs 2022-08-15
niels.fuchs@uni-hamburg.de
'''
import saltharp
import matplotlib.pyplot as plt
from saltharp_cmaps import *


#file_path='../data/SaltharpLog_GroßerTank_Versuch1_vollstandig.log'
file_path='../data/SaltharpLog_GroßerTank_Versuch3.log'

    
#stick_file_path= '../data/TsticksLog_GroßerTank_Versuch1.log'
stick_file_path= '../data/TsticksLog_GroßerTank_Versuch3.log'

ds = saltharp.read_teraterm_logfile(file_path)
    
ds = saltharp.depth_axis(ds, flip_ud=True, vertical_spacing = 0.02, upper_wire_pair_depth = 0.01)

ds = saltharp.replace_T_with_T_Sticks(ds, {0:[4, True]}, stick_file_path)

#ds = saltharp.obtain_retrievals(ds, '2022-07-28 16:07:00', '2022-07-28 16:25:00', salt='SW', gamma_method='no_gamma', sea_salinity=35.)
#ds = saltharp.obtain_retrievals(ds, '2022-08-07 13:55:00', '2022-08-07 14:00:00', salt='SW', gamma_method='no_gamma', sea_salinity=35.)
ds = saltharp.obtain_retrievals(ds, '2022-08-17 06:00:00', '2022-08-17 11:00:00', salt='SW', gamma_method='no_gamma', sea_salinity=25.)
    
ds.to_netcdf(file_path.rsplit('.',1)[0]+'_retrievals.nc')

# plot some data

ds = ds.loc[dict(module=0)] # choose module 0
fig1, ax1 = plt.subplots(2, constrained_layout=True)
ds.S_bulk.T.plot.contourf(ax=ax1[0], cmap=S_cmap,vmin=0,vmax=40,levels=100, cbar_kwargs={'label':'Bulk Salinity'})
ax1[0].invert_yaxis()
ds.phi_l.T.plot.contourf(ax=ax1[1], cmap=psi_cmap,vmin=0,vmax=1,levels=100, cbar_kwargs={'label':'liquid fraction'})
ax1[1].invert_yaxis()
plt.show()


### script bin. Functions used before, uncomment if necessary

'''
### read data again to check for SD and bailout errors, if they correlate with other outtakes

import pandas as pd
import matplotlib.pyplot as plt
import datetime as dt

date_parser = lambda x, y: dt.datetime.strptime(f"{x} {y}",  '[%Y-%m-%d %H:%M:%S.%f]')
data_nans = pd.read_csv(file_path, delimiter=' ', error_bad_lines=False, usecols=[0,1,2], names=['dates','times','dummy'], 
            parse_dates={'date': ['dates', 'times']}, date_parser=date_parser, skiprows=[i for i, line in enumerate(open(file_path)) if not line.startswith('[')])

nan_dates=data_nans.iloc[[i for i, line in data_nans.iterrows() if line['dummy'].startswith('b') or line['dummy'].startswith('F')]]['date']

fig, ax = plt.subplots(1)
for i in range(8):
    ax.plot(data.loc[data['sen']==i, 'date'], data.loc[data['sen']==i, 'R2'], zorder=2)
ax.scatter(nan_dates.values, np.ones(len(nan_dates)), color='grey', alpha=0.5, zorder=1)
plt.show()
'''
