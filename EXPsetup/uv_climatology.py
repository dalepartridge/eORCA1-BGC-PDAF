import os
import numpy as np
import datetime
import netCDF4

_datapath = '/work/n01/n01/dapa/NCEO/eORCA1-BGC-PDAF/OUTPUTS/'
def mask():
    """This function is used to mask the filled values in the climatology files
    The climatology files were generated in the following way:
    1. link all files into one directory
    2. use `cdo ensmean inputs outputfile` command to get the mean value
    The code is not provided here
    """
    for var in ['u', 'v']:
  	    f0 = netCDF4.Dataset(f'{var}_skeb_climatology.nc', 'r+')
        f0[f'{var}o'][:] = np.ma.masked_where(np.abs(f0[f'{var}o'][:])>1e10, f0[f'{var}o'][:])
	    f0.close()

if __name__ == '__main__':
	mask()
