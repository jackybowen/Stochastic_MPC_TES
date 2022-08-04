import pandas as pd
import numpy as np
from sklearn import datasets, linear_model
import sys
import numpy as np
import json
regr = linear_model.LinearRegression()

zonelist  = ['AHU-002','AHU-004']

for i in range(len(zonelist)):
    tab=pd.read_csv('{}.csv'.format(zonelist[i]))

    x=tab.loc[:,('temp_pre', 'tempset', 'tout')]
    y = tab['tempsupply']

    regr.fit(x, y)
    coeffs = {}
    coeffs['temp_pre'] = regr.coef_[0]
    coeffs['tempset'] = regr.coef_[1]    
    coeffs['tout'] = regr.coef_[2] 
    coeffs['intercept'] = regr.intercept_
    print(coeffs)   

    with open('coeff_{}.json'.format(zonelist[i]), 'w', encoding='utf-8') as f:
        json.dump(coeffs, f, indent=4)    
    y_pre1=regr.predict(x)

    tab2=pd.DataFrame()
    tab2['prediction']=y_pre1
    tab2['real']=y
    tab2.to_csv(zonelist[i]+'_comparison.csv')
    
    