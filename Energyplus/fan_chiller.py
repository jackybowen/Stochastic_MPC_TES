import pandas as pd
import numpy as np
from sklearn import datasets, linear_model
import sys
import numpy as np
import json
regr = linear_model.LinearRegression()


for i in range(1,5):
    tab=pd.read_csv('ahu_{}.csv'.format(i))
    tab['flow2'] = tab['flow']*tab['flow']
    tab['flow3'] = tab['flow2']*tab['flow']    
    x=tab.loc[:,('flow', 'flow2', 'flow3')]
    y = tab['power']
    
    if i == 1:
        load =  tab['flow']*(tab['tinlet']-tab['toutlet'])
    else:
        load =  load + tab['flow']*(tab['tinlet']-tab['toutlet'])    

    regr.fit(x, y)
    coeffs = {}
    coeffs['flow'] = regr.coef_[0]
    coeffs['flow2'] = regr.coef_[1]    
    coeffs['flow3'] = regr.coef_[2] 
    coeffs['intercept'] = regr.intercept_
    print(coeffs)

    with open('coeff_ahu_{}.json'.format(i), 'w', encoding='utf-8') as f:
        json.dump(coeffs, f, indent=4)    
    y_pre1=regr.predict(x)

    tab2=pd.DataFrame()
    tab2['prediction']=y_pre1
    tab2['real']=y
    tab2.to_csv('ahu{}_comparison.csv'.format(i))
    
tabs = pd.DataFrame()
tabs['load'] = load
tabs['power'] = tab['chiller_power']
tabs.to_csv('chiller.csv')
tabs['load2'] = tabs['load']*tabs['load']
x=tabs.loc[:,('load', 'load2')]
y = tabs['power']
regr.fit(x, y)
coeffs = {}
coeffs['load'] = regr.coef_[0]
coeffs['load2'] = regr.coef_[1]     
coeffs['intercept'] = regr.intercept_
print(coeffs)

with open('coeff_chiller.json'.format(i), 'w', encoding='utf-8') as f:
    json.dump(coeffs, f, indent=4)    
y_pre1=regr.predict(x)
tab2=pd.DataFrame()
tab2['prediction']=y_pre1
tab2['real']=y
tab2.to_csv('chiller_comparison.csv')