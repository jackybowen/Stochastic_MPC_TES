import pandas as pd
import numpy as np
# import matplotlib
# import matplotlib.dates as mdates
# import datetime
import matplotlib.pyplot as plt
import csv
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
from regex import F

raw_data = pd.read_csv('../Closedloop_test/profile/baseline.csv')

# VAV102 ~VAV121/CORRIDOR => AHU1 (17 zones)
# VAV104/RESTROOM ~VAV116 => AHU3 (6 zones)
raw_data['ts'] = pd.date_range(start='8/1/2017 00:00:00', end='8/5/2017 23:59:00',freq='5min')
print(len(raw_data))
zonelist  = ['VAV-102', 'VAV-118', 'VAV-119','VAV-120','VAV-123A','VAV-123B','VAV-127A','VAV-127B','VAV-129','VAV-131','VAV-133','VAV-136','VAV-142','VAV-143','VAV-150','VAV-100','VAV-121','VAV-104','VAV-105','VAV-107','VAV-108','VAV-112','VAV-116','AHU-002','AHU-004']
zonelist1  = ['VAV-102', 'VAV-118', 'VAV-119','VAV-120','VAV-123A','VAV-123B','VAV-127A','VAV-127B','VAV-129','VAV-131','VAV-133','VAV-136','VAV-142','VAV-143','VAV-150','VAV-CORRIDOR','VAV-RESTROOM','VAV-104','VAV-105','VAV-107','VAV-108','VAV-112','VAV-116','AHU-002','AHU-004']
zonelist1_1 = ['VAV-102', 'VAV-118', 'VAV-119','VAV-120','VAV-123A','VAV-123B','VAV-127A','VAV-127B','VAV-129','VAV-131','VAV-133','VAV-136','VAV-142','VAV-143','VAV-150','VAV-CORRIDOR','VAV-RESTROOM','VAV-104','VAV-105','VAV-107','VAV-108','VAV-112','VAV-116']
zonelist1_2 = ['AHU-002','AHU-004']
tempsetpoints = ['ZONE-'+z+':Zone Thermostat Cooling Setpoint Temperature [C](TimeStep)' for z in zonelist1]
massflowrates = ['ZONE-'+z+' VAV BOX OUTLET NODE:System Node Mass Flow Rate [kg/s](TimeStep)' for z in zonelist1_1] + ['ZONE-'+z+' DIRECT AIR INLET NODE NAME:System Node Mass Flow Rate [kg/s](TimeStep)' for z in zonelist1_2]
tempdischarge = ['ZONE-'+z+' VAV BOX OUTLET NODE:System Node Temperature [C](TimeStep)' for z in zonelist1_1] + ['ZONE-'+z+' DIRECT AIR INLET NODE NAME:System Node Temperature [C](TimeStep)' for z in zonelist1_2]

Mflowrate = raw_data[massflowrates]
Minmax = Mflowrate.agg([min, max])

zonelist_AHU1 = ['VAV-102', 'VAV-118', 'VAV-119','VAV-120','VAV-123A','VAV-123B','VAV-127A','VAV-127B','VAV-129','VAV-131','VAV-133','VAV-136','VAV-142','VAV-143','VAV-150','VAV-CORRIDOR']
zonelist_AHU2 = ['AHU-002']
zonelist_AHU3 = ['VAV-RESTROOM','VAV-104','VAV-105','VAV-107','VAV-108','VAV-112','VAV-116']
zonelist_AHU4 = ['AHU-004']

zonelist_AHUzone1 = ['ZONE-'+z+' VAV BOX OUTLET NODE:System Node Mass Flow Rate [kg/s](TimeStep)' for z in zonelist_AHU1]
zonelist_AHUzone2 = ['ZONE-'+z+' DIRECT AIR INLET NODE NAME:System Node Mass Flow Rate [kg/s](TimeStep)' for z in zonelist_AHU2]
zonelist_AHUzone3 = ['ZONE-'+z+' VAV BOX OUTLET NODE:System Node Mass Flow Rate [kg/s](TimeStep)' for z in zonelist_AHU3]
zonelist_AHUzone4 = ['ZONE-'+z+' DIRECT AIR INLET NODE NAME:System Node Mass Flow Rate [kg/s](TimeStep)' for z in zonelist_AHU4]

zonelist_AHUzone = [zonelist_AHUzone1,zonelist_AHUzone2,zonelist_AHUzone3,zonelist_AHUzone4]

print(zonelist_AHUzone[0])
# for individual AHU
for ahu in range(4):
    Mflowrate = raw_data[zonelist_AHUzone[ahu]]
    Mflowsum = Mflowrate.sum(axis=1)

    X, y = Mflowsum.to_numpy(), raw_data['AHU-00{} FAN:Fan Electric Power [W](TimeStep)'.format(ahu+1)]


    poly = PolynomialFeatures(degree=3, include_bias=False)
    poly_features = poly.fit_transform(X.reshape(-1,1))
    X_train, X_test, y_train, y_test = train_test_split(poly_features, y, test_size=0.3, random_state=42)

    poly_reg_model = LinearRegression()
    poly_reg_model.fit(X_train, y_train)

    yfit = poly_reg_model.predict(X_train)
    y_result = poly_reg_model.predict(X_test)

    RMSE_train= np.sqrt(mean_squared_error(y_train, yfit))
    RMSE_test = np.sqrt(mean_squared_error(y_test, y_result))
    
    #prints relevant error values
    print('\nRMSE: ', RMSE_train, RMSE_test)

    #initializes arrays for plotting
    ts = np.linspace(1, 5, len(y)) #array of days for training data
    t1 = ts[:len(y_train)]
    t2 = ts[len(y_test):] #array of days for test data

    #plots linear regression training results
    plt.figure(ahu)
    plt.plot(t1, y_train, 'o', mfc ='none', label='training data')
    plt.plot(t1, yfit, label='Polynomial regression results')
    plt.xlabel('Day')
    plt.ylabel('Fan electricity power(W)')
    plt.title('Linear regression training results')
    plt.legend()
    plt.show()

    plt.figure(ahu+4)
    print(len(t2))
    print(len(y_test))
    plt.plot(t2, y_test, 'o', mfc ='none', label='test data')
    plt.plot(t2, y_result, label='Polynomial regression results')
    plt.xlabel('Day')
    plt.ylabel('Fan electricity power(W)')
    plt.title('Linear regression test results')
    plt.legend()
    plt.show()


