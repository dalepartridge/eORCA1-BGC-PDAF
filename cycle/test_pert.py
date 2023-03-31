import yaml
import numpy as np
import os

dirname = '/work/n01/n01/ymchen/eORCA1-BGC-PDAF/RUN/EXP-MEDUSA_2000-2020/'
params = ['xaln', 'xald', 'xnln', 'xnld',
          'xvpn', 'xvpd', 'xmetapn', 'xmetapd', 
          'xmpn', 'xmpd', 'xkphn', 'xkphd', 
          'xkmi', 'xkme', 'xgmi', 
          'xgme', 'xbetac', 'xkc']
for i in range(1, 31):
    stream = open(f'{dirname}/ensemble_{i}/fabm.yaml','r')
    f = yaml.safe_load(stream)
    stream.close()
    params = {p:f['instances']['pelagic']['parameters'][p] for p in params}
    print (params)
