#!/usr/bin/env python
# coding: utf-8

# In[26]:


# Author: Bowen Huang
# Created in May 23
#Import necessary modules
import pandas as pd
import numpy as np
import sys
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score
from sklearn.metrics import mean_squared_error

def F2Cfun(f_list):
    "Return list in Celsius degree given Fahrenheit degree"
    return [(f -32)*5/9 for f in f_list]
def Cfm2Kgs(v_list):
    "Return list in Kg/s given Cubic-feet/min(water mass)"
    return [0.47*v for v in v_list]
def AvgLst(v_list, dT):
    return [sum(v_list[i:i+dT])/dT for i in range(0,len(v_list),dT)]
def RMSEp(y_true, y_pred):
    # Normalized RMSE, normalized by (max - min) or mean
    y_t = np.array(y_true)
    y_p = np.array(y_pred)
#     return np.sqrt(np.mean(np.square((y_t - y_p)))) 
    mask = (y_t >= 0.001)
    return np.sqrt(np.mean(np.square((y_t[mask] - y_p[mask])/y_t[mask]))) #np.max(y_t)-np.min(y_t)

def linReg(Xtrain, ytrain, Xtest, ytest, Xtest_m, Xtest_tdis): #linear regression model
    flag = False # Print output or not
    global modelflag
    def predict_acc(Xt):
        yt = []
        for t in range(len(Xt)):
            Tz = reg.predict([Xt[t]])[0]
            if t < len(Xt)-1:
                Xt[t+1][0] = Tz
                Xt[t+1][2] = Xtest_m[t+1]*(Xtest_tdis[t] - Tz)
            yt.append(Tz)
        return np.array(yt)
    
    #performs linear regression for training data
    reg = LinearRegression().fit(Xtrain, ytrain)
    if flag:
        #outputs r^2 scores
        print('Linear regression scores: ', reg.score(Xtrain, ytrain), reg.score(Xtest,ytest))
    
    #predicts output for training data based on regression
    yfit = reg.predict(Xtrain)
    
    err_train = yfit-ytrain #calculates error for training data
    RMSE_train=RMSEp(ytrain,yfit)#np.sqrt(mean_squared_error(ytrain,yfit)) #calculates RMSE for training
    
    #predicts output for test data based on regression
    ######### Compute 1 day multiple days averaged(check if significant change)
    coeff = np.concatenate((reg.coef_,[reg.intercept_]))
    
    if modelflag == 1:
        yresults = predict_acc(Xtest)
    elif modelflag == 2:
        yresults = reg.predict(Xtest)
    
    if flag:
        for i in range(len(coeff)):
            print("{:.16f}".format(coeff[i]))
    
    err_test = yresults-ytest #error for test data
    if modelflag == 1:
        RMSE_test = np.sqrt(mean_squared_error(ytest,yresults))
    elif modelflag == 2:
        RMSE_test= RMSEp(ytest,yresults)#np.sqrt(mean_squared_error(ytest,yresults)) #RMSE for test data
    
    if flag:
        #printts relevant error values
        print('\nRMSE: ', RMSE_train, RMSE_test)
        print('\nMax % error: ', 100*max(np.abs(err_train/ytrain)), 100*max(np.abs(err_test/ytest)))
        print('\nMean % error: ', 100*np.mean(err_train/ytrain), 100*np.mean(err_test/ytest))
        print('\nMean mag % error: ', 100*np.mean(np.abs(err_train/ytrain)), 100*np.mean(np.abs(err_test/ytest)))
    
        #initializes arrays for plotting
        t1 = np.linspace(1, 21, len(ytrain)) #array of days for training data
        t2 = np.linspace(1, 7, len(ytest)) #array of days for test data

        #plots linear regression training results
        plt.figure(1)
        plt.plot(t1, ytrain, 'o', mfc ='none', label='training data')
        plt.plot(t1, yfit, label='linear regression results')
        plt.xlabel('Day')
        plt.ylabel('Zone temperature ($^\circ$C)')
        plt.title('Linear regression training results')
        plt.legend()


        #plots linear regression training % error
        plt.figure(2)
        plt.plot(t1, 100*err_train/ytrain, label='linear regression results')
        plt.xlabel('Day')
        plt.ylabel('% Error')
        plt.title('Linear regression training results')
        plt.legend()

        #plots linear regression test results
        plt.figure(3)
        plt.plot(t2, ytest, 'o', mfc ='none', label='test data')
        plt.plot(t2, yresults, label='linear regression results')
        plt.title('Linear regression test results')
        plt.xlabel('Day')
        plt.ylabel('Zone temperature ($^\circ$C)')
        plt.legend()
        plt.show()
        #plots linear regression test % error
        plt.figure(4)
        plt.plot(t2, 100*err_test/ytest, label='linear regression results')
        plt.title('Linear regression test results')
        plt.xlabel('Day')
        plt.ylabel('% Error')
        plt.legend()
    
    return yresults, coeff, RMSE_test


# In[2]:


rawdata = pd.read_csv('building1.csv')
RoomName = ['002','004','100','102','104','105','107','108',             '112','116','118','119','120','121','123A','123B',            '127A','127B','129','133','136','142','143','150']
# RoomName = RoomName[1:9]
# 002 004 100 121 missing lighting power
N = len(RoomName)
print("Total number of rooms:",N)


# In[3]:


Ttpower = rawdata['WholeBuildingPower']
Ts  = rawdata['Date'].tolist()
TOA = rawdata['OurdoorAirTemperature'].tolist()
TOA = F2Cfun(TOA)
Tzon = []
Tset = []
Mdot = []
Pfan = []
for z in range(N):
    rName = RoomName[z]
    Tzon.append(rawdata['ZoneTemperature_{0}'.format(rName)].tolist()) # unit in Fahrenheit degree
    Tset.append(rawdata['ZoneCoolingTemperatureSetpoint_{0}'.format(rName)].tolist()) # unit in Fahrenheit degree
    Mdot.append(rawdata['SupplyAirFlow_{0}'.format(rName)].tolist()) # unit in cubic-feet/min
#     Pfan.append(rawdata['WholeBuildingPower'] - rawdata['LightingPower_{0}'.format(rName)])
Tzon = [F2Cfun(flist) for flist in Tzon]
Tset = [F2Cfun(flist) for flist in Tset]
Mdot = [Cfm2Kgs(vlist) for vlist in Mdot]
Tdis1 = rawdata['DischargeAirTemperature_112'].tolist()
Tdis1 = F2Cfun(Tdis1)
Tdis2 = rawdata['DischargeAirTemperature_116'].tolist()
Tdis1 = F2Cfun(Tdis2)


# In[28]:


get_ipython().run_line_magic('matplotlib', 'inline')
#qt

Ttotal = 40320
##########################################################
## Set time resolution to 1 or 60(min)
##########################################################
dT = 1# Time interval
##########################################################
## Set model option 1 or 2
##########################################################
modelflag = 2
print('Total Number of Rooms is',N)
RMSE_Ttt = 0
RMSE_Mtt = 0
debugflag = 0 # No printout in linear regression

T0 = []
Tstart = [1,4,10,24,35,78,86,123,137,143]
for tt in Tstart:
    T0.append(int(24*(60/dT)*21+tt*(60/dT)))
Tlen = 24*int(60/dT)
Ttest = []
for t0 in T0:
    Ttest.append(list(range(t0,t0+Tlen)))
for Zid in range(N):
    Mdot_z= Mdot[Zid]
    Tset_z= Tset[Zid]
    Tzon_z= Tzon[Zid]
    Tdis  = [12.8]*len(Tdis1) # Constant discharge temperature 12.8 Degree Celsius

    Tzonm = [sum(Tzon_z[i:i+dT])/dT for i in range(0,len(Tzon_z),dT)]
    TOAm  = [sum(TOA[i:i+dT])/dT for i in range(0,len(TOA),dT)]
    Mdotm = [sum(Mdot_z[i:i+dT])/dT for i in range(0,len(Mdot_z),dT)]
    Tsetm = [sum(Tset_z[i:i+dT])/dT for i in range(0,len(Tset_z),dT)]
    Tdism = [sum(Tdis[i:i+dT])/dT for i in range(0,len(Tdis),dT)]

    Xtrain = []
    Xtest = []
    Xtest_m = []
    Xtest_tdis = []
    # Model 1: Tz' = a1*Tz + a2*Toa + a3*mdot(Tdis - Tz) + a4
    # Model 2: mdot = a1*Tz + a2*Toa + a3*Tset + a4
    for t in range(1,int(Ttotal/dT)):
        if (t < 24*(60/dT)*21): # Use first 3 weeks as training data
            if modelflag == 1:
                Xtrain.append([Tzonm[t-1], TOAm[t-1], Mdotm[t]*(Tdism[t-1] - Tzonm[t-1])])
            elif modelflag == 2:
                Xtrain.append([Tzonm[t-1], TOAm[t-1], Tsetm[t-1]])
        else:
            if modelflag == 1:
                if t in Ttest[9]: # Manually change the Ttest then average RMSE_T
                    Xtest.append([Tzonm[t-1], TOAm[t-1], Mdotm[t]*(Tdism[t-1] - Tzonm[t-1])])
            elif modelflag == 2:
                Xtest.append([Tzonm[t-1], TOAm[t-1], Tsetm[t-1]])
            Xtest_m.append(Mdotm[t])
            Xtest_tdis.append(Tdism[t-1])
    if modelflag == 1:
        ym = Tzonm
    elif modelflag == 2:
        ym = Mdotm
    ytrain= ym[1:(len(Xtrain)+1)]
    ytest = ym[(len(Xtrain)+1):len(Xtrain)+len(Xtest)+1]

    yrlr, coeff, RMSE_test = linReg(Xtrain, ytrain, Xtest, ytest, Xtest_m, Xtest_tdis)
    if modelflag == 1:
        RMSE_Ttt = RMSE_Ttt + RMSE_test
#         Xt_np = np.array(Xtest)
#         Xt1_np = np.array(Xtest_T)
#         Ypred_np = np.array(yrlr)
#         Mz_pred = (Ypred_np - coeff[0]*Xt_np[:,0] - coeff[1]*Xt_np[:,1]-coeff[3])/coeff[2]/Xt1_np
        RMSE_Mtt = 0 #RMSE_Mtt + RMSEp(Xt_np[:,2]/Xt1_np,Mz_pred)#np.sqrt(mean_squared_error(Xt_np[:,2]/Xt1_np,Mz_pred))
    elif modelflag == 2:
        RMSE_Mtt = RMSE_Mtt + RMSE_test
        Xt_np = np.array(Xtest)
        Ypred_np = np.array(yrlr)
        Tz_pred = (Ypred_np - coeff[1]*Xt_np[:,1] - coeff[2]*Xt_np[:,2]-coeff[3])/coeff[0]
        RMSE_Ttt = RMSE_Ttt + RMSEp(Xt_np[:,0],Tz_pred)#np.sqrt(mean_squared_error(Xt_np[:,0],Tz_pred))
#     RMSE_Ptt = RMSE_Ptt + RMSE_P

RMSE_Ttt = RMSE_Ttt/N # Average over all rooms
RMSE_Mtt = RMSE_Mtt/N
print("Zone Temperature prediciton RMSE is",RMSE_Ttt)
print("Mass flow prediciton RMSE is",RMSE_Mtt)


# |dT= 1|              |
# |---------|---------|
# |Time period 1:|    2.036119848146391|
# |Time period 2:|   3.3266630454612582|
# |Time period 3:|    3.5430330659625793|
# |Time period 4:|    1.7355829560805711|
# |Time period 5:|    3.8228522465531594|
# |Time period 6:|    3.069715525288094|
# |Time period 7:|    3.5215623197062196|
# |Time period 8:|    2.5448713493573774|
# |Time period 9:|    7.247195672668298|
# |Time period 10:|    2.5190519429581557|
# |Average: | |

# |dT= 60   |         |
# |---------|---------|
# |Time period 1:|    1.0775063688227344|
# |Time period 2:|    1.516429913271723|
# |Time period 3:|    2.397786684591065|
# |Time period 4:|    1.005062503422326|
# |Time period 5:|    2.382451326587289|
# |Time period 6:|    2.3608151169382148|
# |Time period 7:|    2.013007142159245|
# |Time period 8:|    1.3671103166714278|
# |Time period 9:|    2.032621355378222|
# |Time period 10:|    1.0870957111098922|
# |Average: | |

# In[ ]:


# # Training the Electric power consumption model
# Mtot1 = [0]*len(Mdot_z)
# Mtot2 = [0]*len(Mdot_z)
# Mtot3 = [0]*len(Mdot_z)
# Tr = 

# for Zid in range(N):
#     Mdot_z= Mdot[Zid]
#     Mtot1 = [Mtot1[t] + Mdot_z[t] for t in range(len(Mtot1))
#     Mtot2 = [Mtot2[t] + Mdot_z[t]**2 for t in range(len(Mtot2))
#     Mtot3 = [Mtot3[t] + Mdot_z[t]**3 for t in range(len(Mtot3))
#     Tset_z= Tset[Zid]
#     Tzon_z= Tzon[Zid]
#     Tdis  = [12.8]*len(Tdis1) # Constant discharge temperature 12.8 Degree Celsius


# In[19]:


# import numpy as np
# from sklearn.linear_model import LinearRegression
# X = np.array([[1, 1], [1, 2], [2, 2], [2, 3]])
# # y = 1 * x_0 + 2 * x_1 + 3
# y = np.dot(X, np.array([1, 2])) + 3
# reg = LinearRegression().fit(X, y)
# reg.score(X, y)

# reg.coef_

# reg.intercept_

# yy = reg.predict([X[-1]])
# aa = False
# print(aa)
a = [1,2,4,5,10,11,-1]
a = np.array(a)
mask = (a<=1)
b = a[mask]
print(b)

