import os
import subprocess
from shutil import copyfile

import time
start_time = time.time()




modelPath = 'SEB.idf'
weatherPath = 'USA_WA_Pasco-Tri.Cities.AP.727845_TMY3.epw'

os.environ['BCVTB_HOME'] = 'bcvtb'

print(os.name)
if os.name == 'nt':
   cmdStr = "C:\EnergyPlusV9-0-1\energyplus -w \"%s\" -r \"%s\"" % (weatherPath, modelPath)
else:
   cmdStr = "energyplus -w \"%s\" -r \"%s\"" % (weatherPath, modelPath)

sock = subprocess.Popen('python master_nn_train.py'+' '+modelPath+' '+'ep_SEB.config', shell=True)
simulation = subprocess.Popen(cmdStr, shell=True)
simulation.wait()
sock.terminate()
copyfile('eplusout.csv', 'train.csv')
print("--- %s seconds ---" % (time.time() - start_time))