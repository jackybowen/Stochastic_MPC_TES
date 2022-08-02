import pandas as pd
import numpy as np
import matplotlib
import matplotlib.dates as mdates
import datetime
import matplotlib.pyplot as plt
import csv

from regex import F

raw_data = pd.read_csv('baseline.csv')


raw_data['ts'] = pd.date_range(start='8/1/2017 00:00:00', end='8/5/2017 23:59:00',freq='5min')
print(len(raw_data))
zonelist  = ['VAV-102', 'VAV-118', 'VAV-119','VAV-120','VAV-123A','VAV-123B','VAV-127A','VAV-127B','VAV-129','VAV-131','VAV-133','VAV-136','VAV-142','VAV-143','VAV-150','VAV-100','VAV-121','VAV-104','VAV-105','VAV-107','VAV-108','VAV-112','VAV-116','AHU-002','AHU-004']
zonelist1  = ['VAV-102', 'VAV-118', 'VAV-119','VAV-120','VAV-123A','VAV-123B','VAV-127A','VAV-127B','VAV-129','VAV-131','VAV-133','VAV-136','VAV-142','VAV-143','VAV-150','VAV-CORRIDOR','VAV-RESTROOM','VAV-104','VAV-105','VAV-107','VAV-108','VAV-112','VAV-116','AHU-002','AHU-004']
zonelist1_1 = ['VAV-102', 'VAV-118', 'VAV-119','VAV-120','VAV-123A','VAV-123B','VAV-127A','VAV-127B','VAV-129','VAV-131','VAV-133','VAV-136','VAV-142','VAV-143','VAV-150','VAV-CORRIDOR','VAV-RESTROOM','VAV-104','VAV-105','VAV-107','VAV-108','VAV-112','VAV-116']
zonelist1_2 = ['AHU-002','AHU-004']
tempsetpoints = ['ZONE-'+z+':Zone Thermostat Cooling Setpoint Temperature [C](TimeStep)' for z in zonelist1]
massflowrates = ['ZONE-'+z+' VAV BOX OUTLET NODE:System Node Mass Flow Rate [kg/s](TimeStep)' for z in zonelist1_1] + [z+' SUPPLY EQUIPMENT OUTLET NODE:System Node Mass Flow Rate [kg/s](TimeStep)' for z in zonelist1_2]

# print(len(massflowrates))
# print(raw_data[massflowrates[10]])
# Tsetpoint = raw_data[tempsetpoints.append('Date/Time')]
# Tsetpoint.to_csv('baseline_setpoint.csv',index=False)

Mflowrate = raw_data[massflowrates]
Minmax = Mflowrate.agg([min, max])
print(Minmax)
Mflowrates = pd.concat([Mflowrate, Minmax])
Minmax.to_csv('baseline_mflowbounds.csv')
