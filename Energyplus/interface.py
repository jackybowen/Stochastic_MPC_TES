from genericpath import exists
import json
from attr import NOTHING
import requests
import time

class DataSampler:
    def __init__(self):
        self.skipheader = 5 # Number of lines to skip in the header
        self.Windowlen = 12  # length of control window for data sampling
        self.Numzones = 25  # Number of zones
        self.NumAHUs = 4     # Number of AHUs
        key_list = ["Time","OutdoorAirTemperature"]
        for i in range(self.Numzones):
            key_list.append("Tzon"+str(i+1))
        for i in range(self.Numzones):
            key_list.append("Mflow"+str(i+1))
        for i in range(self.Numzones):
            key_list.append("Tzondis"+str(i+1))
        for i in range(self.NumAHUs):
            key_list.append("Tahudis"+str(i+1))
            key_list.append("Sflow"+str(i+1))
            key_list.append("Tret"+str(i+1))
            key_list.append("Tmix"+str(i+1))
        # key_list.append("OATforecast")
        self.keylist = key_list
        ##########################################################
        # 0~4:      header
        # 5:        time
        # 6~30:     Zone mean air temperature
        # 31~56:    Zone Air flow
        # 57~81:    Zone discharge Air Temperature
        # 82~97:    AHU level DischargeAirTemp, SupplyAirFlow, ReturnAirTemp, MixedAirTemp.
        ##########################################################

    def work(self, url):
        ''' Run the model in its directory. '''
        # Delete output file every run if it exists
    #	u={}
        with open('./input.json') as f:
            data = f.read()
        u = json.loads(data)
    #	json_object = json.dumps(u, indent = 4,default=str)
        result = requests.post('{0}/calculate'.format(url), json=u, timeout=300).json()
        # while True:
        #     time.sleep(60) # sixty second delay
        #     try:
        #         isexist(file_name)
        #         break
        #     except Error:
        #         print ("Result not ready, trying again in one minute")
        # with open('./output.json','w') as f:
        # 	json.dump(result, f)
        # print result
        return result

    def sample2dict(self, Outputdict, Inputlist, Windowlen = NOTHING):
        keys = Outputdict.keys()
        key_list = list(keys)

        if Windowlen is NOTHING:
            Windowlen = self.Windowlen
        for i in range(len(key_list)):
            Outputdict[key_list[i]].append(float(Inputlist[i+self.skipheader]))
        n = len(Outputdict[key_list[0]])
            
        return Outputdict