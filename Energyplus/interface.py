import json
import requests

class DataSampler:
    def __init__(self):
        self.skipheader = 5 # Number of lines to skip in the header
        self.Windowlen = 1  # length of window for data sampling
        self.Numzones = 25  # Number of zones
        self.NumAHUs = 4     # Number of AHUs
        key_list = ["time","OutdoorAirTemperature"]
        for i in range(self.Numzones):
            key_list.append("Tzon"+str(i))
        for i in range(self.Numzones):
            key_list.append("Mflow"+str(i))
        for i in range(self.NumAHUs):
            key_list.append("Tahudis"+str(i))
            key_list.append("Sflow"+str(i))
            key_list.append("Tret"+str(i))
            key_list.append("Tmix"+str(i))
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
        result = requests.post('{0}/calculate'.format(url), json=u).json()
        # with open('./output.json','w') as f:
        # 	json.dump(result, f)
        # print result
        return result

    def sample2dict(self, Outputdict, Inputlist, Windowlen):
        keys = Outputdict.keys()
        key_list = list(keys)
        n = len(Outputdict[key_list[0]])
        if len(Windowlen) == 0:
            Windowlen = self.Windowlen
        if n < Windowlen:
            for i in range(len(key_list)):
                Outputdict[key_list[i]].append(float(Inputlist[i+self.skipheader]))
        else:
            with open('./input.json','w') as f:
                json.dump(Outputdict, f)
            Outputdict = {k:[] for k in self.keylist}
            result = self.work('http://127.0.0.1:5000')
        return Outputdict, result