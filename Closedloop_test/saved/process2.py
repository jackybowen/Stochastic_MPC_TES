import pandas as pd
import numpy as np
import matplotlib
import matplotlib.dates as mdates
import datetime
import matplotlib.pyplot as plt
import csv


tt_power_df = pd.DataFrame(np.zeros([48,2]), columns = ['upp_cost','low_cost'])

# tt_power_df['tt_power_cost'] = mpc_data[''].sum(axis=1)
for t in range(48):
    mpc_data = pd.read_csv('sol_T{}.csv'.format(str((t+1)*60)))
    tt_power_df['upp_cost'][t] = sum(sum([mpc_data['upp_penalty_cost_{}'.format(str(z+1))] for z in range(25)]))
    tt_power_df['low_cost'][t] = sum(sum([mpc_data['low_penalty_cost_{}'.format(str(z+1))] for z in range(25)]))
tt_power_df.to_csv('result1.csv')