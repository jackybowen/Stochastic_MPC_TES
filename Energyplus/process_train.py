import pandas as pd
import numpy as np


tab = pd.read_csv('train.csv')

vavs = []

tout = 'Environment:Site Outdoor Air Drybulb Temperature [C](TimeStep)'

zonelist  = ['VAV-102', 'VAV-118', 'VAV-119','VAV-120','VAV-123A','VAV-123B','VAV-127A','VAV-127B','VAV-129','VAV-131','VAV-133','VAV-136','VAV-142','VAV-143','VAV-150','VAV-100','VAV-121','VAV-104','VAV-105','VAV-107','VAV-108','VAV-112','VAV-116','AHU-002','AHU-004']
zonelist1  = ['VAV-102', 'VAV-118', 'VAV-119','VAV-120','VAV-123A','VAV-123B','VAV-127A','VAV-127B','VAV-129','VAV-131','VAV-133','VAV-136','VAV-142','VAV-143','VAV-150','VAV-CORRIDOR','VAV-RESTROOM','VAV-104','VAV-105','VAV-107','VAV-108','VAV-112','VAV-116','AHU-002','AHU-004']
zonelist1_1 = ['VAV-102', 'VAV-118', 'VAV-119','VAV-120','VAV-123A','VAV-123B','VAV-127A','VAV-127B','VAV-129','VAV-131','VAV-133','VAV-136','VAV-142','VAV-143','VAV-150','VAV-CORRIDOR','VAV-RESTROOM','VAV-104','VAV-105','VAV-107','VAV-108','VAV-112','VAV-116']
zonelist1_2 = ['AHU-002','AHU-004']
tempsetpoints = ['ZONE-'+z+':Zone Thermostat Cooling Setpoint Temperature [C](TimeStep)' for z in zonelist1]
massflowrates = ['ZONE-'+z+' VAV BOX OUTLET NODE:System Node Mass Flow Rate [kg/s](TimeStep)' for z in zonelist1_1] + ['ZONE-'+z+' DIRECT AIR INLET NODE NAME:System Node Mass Flow Rate [kg/s](TimeStep)' for z in zonelist1_2]
tempdischarge = ['ZONE-'+z+' VAV BOX OUTLET NODE:System Node Temperature [C](TimeStep)' for z in zonelist1_1] + ['ZONE-'+z+' DIRECT AIR INLET NODE NAME:System Node Temperature [C](TimeStep)' for z in zonelist1_2]

zonetemps =  ['ZONE-'+z+':Zone Mean Air Temperature [C](TimeStep)' for z in zonelist1]

for i in range(len(zonetemps)):
    temp = tab[zonetemps[i]].to_list()
    tempsupply = tab[tempdischarge[i]].to_list()
    mdot = tab[massflowrates[i]].to_list()
    tempset = tab[tempsetpoints[i]].to_list()
    touttemp = tab[tout].to_list()
    tabs = pd.DataFrame()
    tabs['temp'] = temp
    tabs['tempsupply'] = tempsupply    
    tabs['tempset'] = tempset
    tabs['mdot'] = mdot  
    tabs['tout'] = touttemp 
    tabs = tabs.groupby(np.arange(len(tabs))//12).mean()
    tabs2 = pd.DataFrame()
    tabs2['temp_pre'] = tabs['temp'][1:].to_list()
    tabs2['tempsupply'] = tabs['tempsupply'][:-1].to_list()    
    tabs2['tempset'] = tabs['tempset'][:-1].to_list()   
    tabs2['mdot'] = tabs['mdot'][:-1].to_list()    
    tabs2['tout'] = tabs['tout'][:-1].to_list()    
    tabs2.to_csv('{}.csv'.format(zonelist1[i]))
    
ahu_mflow = ['AHU-00'+str(z)+' SUPPLY EQUIPMENT OUTLET NODE:System Node Mass Flow Rate [kg/s](TimeStep)' for z in [1,2,3,4]]

fanpower = ['AHU-00'+str(z)+' FAN:Fan Electric Power [W](TimeStep)' for z in [1,2,3,4]]

tmixs = ['AHU-00'+str(z)+' OA COOLCNODE:System Node Temperature [C](TimeStep)' for z in [1,2,3,4]]

tdiss = ['AHU-00'+str(z)+' SUPPLY EQUIPMENT OUTLET NODE:System Node Temperature [C](TimeStep)' for z in [1,2,3,4]]

for i in range(len(ahu_mflow)):
    flow = tab[ahu_mflow[i]].to_list()
    power = tab[fanpower[i]].to_list()
    tinlet = tab[tmixs[i]].to_list()
    toutlet = tab[tdiss[i]].to_list()
    tabs = pd.DataFrame()
    tabs['power'] = power
    tabs['flow'] = flow    
    tabs['tinlet'] = tinlet
    tabs['toutlet'] = toutlet 
    tabs = tabs.groupby(np.arange(len(tabs))//12).mean()
    tabs.to_csv('ahu_{}.csv'.format(i+1))