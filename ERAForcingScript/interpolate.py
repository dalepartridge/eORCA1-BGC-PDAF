import os
import cfgrib
import netCDF4
import numpy as np
import xarray as xr

_varnames = [['u10', 'v10', 't2m', 'msl', 'd2m', 'sp', ],
            ['msr', 'msdwlwrf', 'msdwswrf', 'mtpr', ]]

def getSPH(ens, year):
    f = netCDF4.Dataset(os.path.join(f'ERA5_ens_{ens}', f'ERA5_d2m_y{year}.nc'), 'r')
    d2m = f['d2m'][-1]
    f.close()
    f = netCDF4.Dataset(os.path.join(f'ERA5_ens_{ens}', f'ERA5_sp_y{year}.nc'), 'r')
    sp = f['sp'][-1]
    f.close()
    esat = 611.21 * np.exp( 17.502 * (d2m-273.16) / (d2m-32.19) )
    dyrvap = 287.0597 / 461.5250
    f = netCDF4.Dataset(os.path.join(f'ERA5_ens_{ens}', f'ERA5_SPH_y{year}.nc'), 'r+')
    f['SPH'][-1, :, :] = dyrvap * esat / ( sp - (1-dyrvap) * esat)
    f.close()


def removeLaststep(ens, year):
    # remove last time steps
    for vname in _varnames[1]:
        os.system(f'ncks -d time,0,729 ERA5_ens_{ens}/ERA5_{vname}_y{year}.nc ERA5_ens_{ens}/ERA5_{vname}_y{year}-short.nc')
        os.system(f'mv ERA5_ens_{ens}/ERA5_{vname}_y{year}-short.nc ERA5_ens_{ens}/ERA5_{vname}_y{year}.nc')


def getFilename(dirname, year):
    return os.path.join(dirname, f'download{year}.grib')


def interpolate(i, year, ne):
    """argument i is the index from cfgridb datasets
    i = 0 is usually 3 hourly data
    i = 1 is usually 12 hourly data
    but it may vary with the data, so check it first.
    """
    if i == 1:
        for ens in range(ne):
            getSPH(ens, year)
            removeLaststep(ens, year)

    if year == 2015:
        if i == 0:
            f0 = cfgrib.open_datasets('/work/n01/n01/dapa/NCEO/eORCA1-ERA5/download.grib',
                                      backend_kwargs={
                                      'indexpath':'/work/n01/n01/ymchen/ERA5Forcing/download2015.idx'
                                      })[0]
        elif i == 1:
            f0 = cfgrib.open_datasets('/work/n01/n01/dapa/NCEO/eORCA1-ERA5/download_2.grib',
                                      backend_kwargs={
                                      'indexpath':'/work/n01/n01/ymchen/ERA5Forcing/download2015_2.idx'
                                      })[0]
    else:
        f0 = cfgrib.open_datasets(getFilename(dirname, year))[i]

    f1 = cfgrib.open_datasets(getFilename(dirname, year+1))[i]
    if i == 1:
        f0 = f0.isel(step=1)
        f1 = f1.isel(step=1)

    for ens in range(ne):
        for vname in _varnames[i]:
            f = netCDF4.Dataset(os.path.join(f'ERA5_ens_{ens}', f'ERA5_{vname}_y{year}.nc'), 'r+')
            ds0 = f0.isel(number=ens)[vname]
            ds1 = f1.isel(number=ens)[vname]
            f[vname][-1, :, :] = 0.5*(ds0[-1 - i] + ds1[0])
            print (f'{ens+1}-th ensemble {vname} done')
            f.close()
    f0.close()
    f1.close()


if __name__ == '__main__':
    pass
