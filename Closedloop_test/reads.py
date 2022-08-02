import sqlite3
import pandas as pd
import datetime
import os.path
import numpy as np
import matplotlib.pyplot as plt
import json

query_template = query_template = "SELECT ts, value_string FROM data " \
                                    "INNER JOIN topics ON data.topic_id = topics.topic_id " \
                                    "WHERE topics.topic_name = \"record/{tnc}/{bldg}/{device}/{point}\""

def gen_query(tnc, bldg, device, point):
    return query_template \
        .replace("{tnc}", tnc) \
        .replace("{bldg}", bldg) \
        .replace("{device}", device) \
        .replace("{point}", point)  
  
s = 0

for i in range(1,13):

   for j in range(1,32):
         
       if i<10:
       
             month = '0{}'.format(i)
             
          
       else:
       
             month = '{}'.format(i) 

       if j <10:

             day =  '0{}'.format(j)

       else:

             day =  '{}'.format(j)
       
   
       fname = 'small_office.historian.sqlite'.format(month,day,month,day)
   
       if os.path.isfile(fname):

            database = sqlite3.connect(fname)
        
            # query = gen_query('PNNL', 'SMALL_OFFICE', 'METERS','WholeBuildingPower')
            # for z in [5]:
            zoneid = 'HP'+str(1)
            # query = gen_query('PNNL', 'SMALL_OFFICE', zoneid,'ZoneTemperature')
            # query = gen_query('PNNL', 'SMALL_OFFICE', zoneid,'ZoneCoolingTemperatureSetPoint')
            query = gen_query('PNNL', 'SMALL_OFFICE','SMALL OFFICE','BuildingPower')

            # query = gen_query('PNNL', 'SMALL_OFFICE', zoneid,'Power')
            # print(query)
            
            if s <1:
                  out_t = pd.read_sql_query(query, database)
            else:
                  temp = pd.read_sql_query(query, database)
                  out_t = out_t.merge(temp, how='outer')
                  # out_t = out_t.rename({'value_string': zoneid+'ZoneTemperature'}, axis = 'columns')
            
            s = s + 1
            print(len(out_t))
            
out_t['ts'] =  pd.to_datetime(out_t['ts'], errors='coerce')

def fun1(str_row):
      out_str = json.loads(str_row)
      out_str = out_str[0]["Target"]
      return out_str

def fun2(str_row):
      out_str = json.loads(str_row)
      out_t = out_str[0]["Timestamp"]
      return out_t

out_t['ts'] = out_t['value_string'].apply(fun2)
out_t['value_string'] = out_t['value_string'].apply(fun1)
out_t = out_t[~out_t['value_string'].str.contains("None",na=False)]
out_t.drop_duplicates()

out_t.groupby(out_t["ts"].dt.hour)["value_string"].mean().plot(kind='bar', rot=0, ax=axs)
plt.xlabel("Hour of the day")
plt.ylabel("$Building power$")
# out_column = out_t['value_string']
# out_str = ""
# for ind, val in out_column.items():
#       m = re.search('"Target": (.*)},',out_t.at[ind, 'value_string'])
#       out_str = m.group(1)
#       print(out_str)
#       out_t.at[ind, 'value_string'] = out_str
#       print(out_t.at[ind, 'value_string'])
#       # print(f"Index: {ind}, Value: {val}")


# out_avg = out_t['value_string'][1:].values.reshape(-1,60)
# # out_t['ts'] =  out_t['ts']-datetime.timedelta(hours=7)
# out_avg = np.squeeze(np.average(out_avg.astype(np.float), axis=1))

# pd.DataFrame(out_avg).to_csv('smalloffice_power_hourly.csv')
out_t.to_csv('smalloffice_power_record.csv')

# out_t.to_csv('smalloffice_zonetemp_'+zoneid+'.csv')
# out_t.to_csv('smalloffice_setpoint_'+zoneid+'.csv')
# out_t.to_csv('smalloffice_power_'+zoneid+'.csv')
# pd.DataFrame(out_avg).to_csv('smalloffice_zonetemp_'+zoneid+'.csv')
# pd.DataFrame(out_avg).to_csv('smalloffice_setpoint_'+zoneid+'.csv')
plt