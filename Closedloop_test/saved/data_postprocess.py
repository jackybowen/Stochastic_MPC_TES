import pandas as pd
import numpy as np
import matplotlib
import matplotlib.dates as mdates
import datetime
import matplotlib.pyplot as plt
import csv
# from sklearn.preprocessing import PolynomialFeatures
# from sklearn.linear_model import LinearRegression
# from sklearn.model_selection import train_test_split
# from sklearn.metrics import mean_squared_error
# from regex import F
baseline = pd.read_csv('../../Energyplus/data_process/baseline.csv')
mpc_data = pd.read_csv('sol_T60.csv')[1:]
raw_data = pd.read_csv('../../Energyplus/data_process/eplusout.csv')

# VAV102 ~VAV121/CORRIDOR => AHU1 (17 zones)
# VAV104/RESTROOM ~VAV116 => AHU3 (6 zones)
raw_data['ts'] = pd.date_range(start='8/1/2017 00:05:00', end='8/15/2017 00:00:00',freq='5min')
start_ts = '8/1/2017 00:05:00'
end_ts = '8/2/2017 00:00:00'
raw_data = raw_data.loc[raw_data["ts"].between(start_ts, end_ts)]

raw_data = raw_data.groupby(np.arange(len(raw_data))//12).mean()
# print(len(raw_data))
# print(len(pd.date_range(start=start_ts, end=end_ts,freq='1H')))
ts = pd.date_range(start=start_ts, end=end_ts,freq='1H')
raw_data['ts'] = ts

price_1day = pd.read_csv('../profile/daily_prices.csv').iloc[::60,:]
price_14days = pd.concat([price_1day]*14, ignore_index=True)
price_14days['ts'] = pd.date_range(start='8/1/2017 00:05:00', end='8/15/2017 00:00:00',freq='1H')
price_1day = price_14days.loc[price_14days["ts"].between(start_ts, end_ts)]
print(price_1day)

zonelist  = ['VAV-102', 'VAV-118', 'VAV-119','VAV-120','VAV-123A','VAV-123B','VAV-127A','VAV-127B','VAV-129','VAV-131','VAV-133','VAV-136','VAV-142','VAV-143','VAV-150','VAV-100','VAV-121','VAV-104','VAV-105','VAV-107','VAV-108','VAV-112','VAV-116','AHU-002','AHU-004']
zonelist1  = ['VAV-102', 'VAV-118', 'VAV-119','VAV-120','VAV-123A','VAV-123B','VAV-127A','VAV-127B','VAV-129','VAV-131','VAV-133','VAV-136','VAV-142','VAV-143','VAV-150','VAV-CORRIDOR','VAV-RESTROOM','VAV-104','VAV-105','VAV-107','VAV-108','VAV-112','VAV-116','AHU-002','AHU-004']
zonelist1_1 = ['VAV-102', 'VAV-118', 'VAV-119','VAV-120','VAV-123A','VAV-123B','VAV-127A','VAV-127B','VAV-129','VAV-131','VAV-133','VAV-136','VAV-142','VAV-143','VAV-150','VAV-CORRIDOR','VAV-RESTROOM','VAV-104','VAV-105','VAV-107','VAV-108','VAV-112','VAV-116']
zonelist1_2 = ['AHU-002','AHU-004']

pr = ['price']
oat = ['OutdoorTemperature']
tempsetpoints1 = ['sol_tzon_'+str(z+1) for z in range(25)]
massflowrates1 = ['sol_mflow_'+str(z+1) for z in range(25)]
slack_u = ['Tlow_'+str(z+1) for z in range(25)]
slack_d = ['Tupp_'+str(z+1) for z in range(25)]
# Tupp = ['Tupp_'+str(z+1) for z in range(25)]
# Tlow = ['Tlow_'+str(z+1) for z in range(25)]

tempsetpoints = ['ZONE-'+z+':Zone Thermostat Cooling Setpoint Temperature [C](TimeStep)' for z in zonelist1]
massflowrates = ['ZONE-'+z+' VAV BOX OUTLET NODE:System Node Mass Flow Rate [kg/s](TimeStep)' for z in zonelist1_1] + ['ZONE-'+z+' DIRECT AIR INLET NODE NAME:System Node Mass Flow Rate [kg/s](TimeStep)' for z in zonelist1_2]
tempdischarge = ['ZONE-'+z+' VAV BOX OUTLET NODE:System Node Temperature [C](TimeStep)' for z in zonelist1_1] + ['ZONE-'+z+' DIRECT AIR INLET NODE NAME:System Node Temperature [C](TimeStep)' for z in zonelist1_2]


dates =  matplotlib.dates.date2num(raw_data['ts'].tolist())

# Tset = raw_data[tempsetpoints]
# print(Tzone)

zonelist_AHU1 = ['VAV-102', 'VAV-118', 'VAV-119','VAV-120','VAV-123A','VAV-123B','VAV-127A','VAV-127B','VAV-129','VAV-131','VAV-133','VAV-136','VAV-142','VAV-143','VAV-150','VAV-CORRIDOR']
zonelist_AHU2 = ['AHU-002']
zonelist_AHU3 = ['VAV-RESTROOM','VAV-104','VAV-105','VAV-107','VAV-108','VAV-112','VAV-116']
zonelist_AHU4 = ['AHU-004']

zonelist_AHUzone1 = ['ZONE-'+z+' VAV BOX OUTLET NODE:System Node Mass Flow Rate [kg/s](TimeStep)' for z in zonelist_AHU1]
zonelist_AHUzone2 = ['ZONE-'+z+' DIRECT AIR INLET NODE NAME:System Node Mass Flow Rate [kg/s](TimeStep)' for z in zonelist_AHU2]
zonelist_AHUzone3 = ['ZONE-'+z+' VAV BOX OUTLET NODE:System Node Mass Flow Rate [kg/s](TimeStep)' for z in zonelist_AHU3]
zonelist_AHUzone4 = ['ZONE-'+z+' DIRECT AIR INLET NODE NAME:System Node Mass Flow Rate [kg/s](TimeStep)' for z in zonelist_AHU4]

zonelist_AHUzone = [zonelist_AHUzone1,zonelist_AHUzone2,zonelist_AHUzone3,zonelist_AHUzone4]

# print(zonelist_AHUzone[0])
##  for individual AHU
for z in range(25):
    fig, (ax1,ax2) = plt.subplots(2,1, sharex=True)
    ax1.plot_date(dates,price_1day['price'],linestyle='-',label='Price',color='k')
    ax2.plot_date(dates,raw_data[tempsetpoints[z]],linestyle='-',label='MPC-Tset-{}'.format(zonelist1[z]),color='g')
    ax2.plot_date(dates,mpc_data[tempsetpoints1[z]],linestyle='-',label='MPC-Tpred-{}'.format(zonelist1[z]),color='c')
    ax2.plot_date(dates,mpc_data[slack_d[z]],linestyle='--',label='T_lower-{}'.format(zonelist1[z]),color='r')
    ax2.plot_date(dates,mpc_data[slack_u[z]],linestyle='--',label='T_upper-{}'.format(zonelist1[z]),color='r')
    ax2.plot_date(dates,baseline['ZONE-{}:Zone Thermostat Cooling Setpoint Temperature [C](TimeStep)'.format(zonelist1[z])][:len(dates)],linestyle='-',label='Baseline-Tset-{}'.format(zonelist1[z]),color='b')
    # ax1.plot_date(dates,prediction['sol_tzon_{}'.format(i+1)][:24-start_hour-1],linestyle='-.',label='pre',marker='.')
    # ax1.plot_date(dates[start_hour+1:],prediction['sol_tzon_{}_corr'.format(i+1)][:24-start_hour-1],linestyle='-.',label='cor',marker='.')
    # ax1.plot_date(dates,baseline['{}:Schedule Value [](TimeStep)'.format(points[0].upper())]-2,linestyle='-.',label='control',marker='')
    # ax1.plot_date(dates,baseline['{}:Schedule Value [](TimeStep)'.format(points[0].upper())]+2,linestyle='-.',label='control',marker='')
    # ax2.plot_date(dates,baseline['{}:Schedule Value [](TimeStep)'.format(points[0].upper())]-2,linestyle='-.',label='control',marker='')
    # ax2.plot_date(dates,baseline['{}:Schedule Value [](TimeStep)'.format(points[0].upper())]+2,linestyle='-.',label='control',marker='')
    ax1.legend()
    ax2.legend()
    myFmt = mdates.DateFormatter('%H:%M')
    ax1.xaxis.set_major_formatter(myFmt)
    plt.savefig('Zone-{}.png'.format(zonelist[z]))

ahu_mflow = ['AHU-00'+str(z)+' SUPPLY EQUIPMENT OUTLET NODE:System Node Mass Flow Rate [kg/s](TimeStep)' for z in [1,2,3,4]]

fanpower = ['AHU-00'+str(z)+' FAN:Fan Electric Power [W](TimeStep)' for z in [1,2,3,4]]

tmixs = ['AHU-00'+str(z)+' OA COOLCNODE:System Node Temperature [C](TimeStep)' for z in [1,2,3,4]]

tdiss = ['AHU-00'+str(z)+' SUPPLY EQUIPMENT OUTLET NODE:System Node Temperature [C](TimeStep)' for z in [1,2,3,4]]

chiller = 'SEB CHILLER:Chiller Electric Power [W](TimeStep)'




# ################################################################### Energy saving analysis ###########################################################################################
# tt_power_mpc = 0
# tt_power_baseline = 0
# tt_power_df = pd.DataFrame(np.zeros([14,2]),index = pd.date_range(start='8/1/2017', end='8/14/2017',freq='1D'), columns = ['MPC total power Integration', 'Baseline total power Integration'])
# for f in range(4):
#     tt_power_mpc = tt_power_mpc + sum(raw_data['AHU-00'+str(f+1)+' FAN:Fan Electric Power [W](TimeStep)'])
#     tt_power_baseline = tt_power_baseline + sum(baseline['AHU-00'+str(f+1)+' FAN:Fan Electric Power [W](TimeStep)'])
#     fig, (ax1) = plt.subplots(1, sharex=True)
#     ax1.plot_date(dates,raw_data['AHU-00'+str(f+1)+' FAN:Fan Electric Power [W](TimeStep)'],linestyle='-',label='MPC',color='b')
#     ax1.plot_date(dates,baseline['AHU-00'+str(f+1)+' FAN:Fan Electric Power [W](TimeStep)'],linestyle='-',label='Baseline',color='r')
#     ax1.legend()
#     myFmt = mdates.DateFormatter('%H:%M')
#     ax1.xaxis.set_major_formatter(myFmt)
#     plt.savefig('AHU-{}.png'.format(str(f+1)))
# tt_power_mpc = tt_power_mpc + sum(raw_data['SEB CHILLER:Chiller Electric Power [W](TimeStep)'])
# tt_power_baseline = tt_power_baseline + sum(baseline['SEB CHILLER:Chiller Electric Power [W](TimeStep)'])
# print('Total power of MPC: ', tt_power_mpc,'\n')
# print('Total power of Baseline: ', tt_power_baseline,'\n')

# hindex = range(288)
# for h in range(14):
#     for f in range(4):
#         tt_power_df.iloc[h,0] = tt_power_df.iloc[h,0] + raw_data.loc[hindex,'AHU-00'+str(f+1)+' FAN:Fan Electric Power [W](TimeStep)'].sum()
#         tt_power_df.iloc[h,1] = tt_power_df.iloc[h,1] + baseline.loc[hindex,'AHU-00'+str(f+1)+' FAN:Fan Electric Power [W](TimeStep)'].sum()
#     print('AHU ',f,' and Day ',h,' is \n',tt_power_df)
#     tt_power_df.iloc[h,0] = tt_power_df.iloc[h,0] + raw_data.loc[hindex,'SEB CHILLER:Chiller Electric Power [W](TimeStep)'].sum()
#     tt_power_df.iloc[h,1] = tt_power_df.iloc[h,1] + baseline.loc[hindex,'SEB CHILLER:Chiller Electric Power [W](TimeStep)'].sum()
#     hindex = [hind + 288 for hind in hindex]
# tt_power_df['Energy saving Ratio'] = tt_power_df['MPC total power Integration']/tt_power_df['Baseline total power Integration']
# tt_power_df.to_csv('Daily_energy_saving.csv')
