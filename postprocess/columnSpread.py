import numpy as np
import netCDF4
import matplotlib.pyplot as plt

f = netCDF4.Dataset(f'exp-logmean-PC/2015/01/state_20150116.nc', 'r')
ny, nx = f.dimensions['y'].size, f.dimensions['x'].size
print (ny, nx)
incr = (f['PHD'][1, 0] - f['PHD'][0, 0]).ravel()
i = np.nanargmax(np.abs(incr))
mj = i//nx
mi = i - mj*nx
f.close()

filename = 'eORCA1_1d_20150101_20150331_ptrc1_T_201501-201501.nc'
f = netCDF4.Dataset(f'exp-logmean-PC/2015/01/{filename}', 'r')

print (f['PHD_std'])
plt.figure(1)
plt.pcolormesh(np.ma.masked_where(f['PHD_std'][15, :, :, mi] > 1e10, f['PHD_std'][15, :, :, mi]))
plt.colorbar()
plt.savefig('test.png')
#fens = [netCDF4.Dataset(f'{expname}/{year}/{cycle_month}/ensemble_{i}/{filename}', 'r') for i in range(1, _Ne+1)]
