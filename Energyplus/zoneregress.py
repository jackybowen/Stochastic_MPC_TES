import pandas as pd
import numpy as np
from sklearn import datasets, linear_model
import sys
import numpy as np
import json
regr = linear_model.LinearRegression()

zonelist  = ['VAV-102', 'VAV-118', 'VAV-119','VAV-120','VAV-123A','VAV-123B','VAV-127A','VAV-127B','VAV-129','VAV-131','VAV-133','VAV-136','VAV-142','VAV-143','VAV-150','VAV-CORRIDOR','VAV-RESTROOM','VAV-104','VAV-105','VAV-107','VAV-108','VAV-112','VAV-116']

for i in range(len(zonelist)):
    tab=pd.read_csv('{}.csv'.format(zonelist[i]))

    x=tab.loc[:,('temp_pre', 'tempset', 'tout')]
    y = tab['mdot']

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
    
    