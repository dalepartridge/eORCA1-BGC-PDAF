import numpy as np
import pyPDAF.PDAF as PDAF
import config

dim = 30
dim_ens = 1000
element = 0
oens = np.random.random((dim, dim_ens))
obs = oens[:, 8]
CRPS, reli, resol, uncert, status = PDAF.diag_crps(element, oens, obs)
print (status)