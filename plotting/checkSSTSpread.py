import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.colors as mcolors
import matplotlib.gridspec
import cmocean

_Ne = 30
_nx = 362
_ny = 332

def getMeanSTD(it, varname):
    fileformat = '2015/01/ensemble_{i}/eORCA1_1d_20150101_20150111_grid_T_201501-201501.nc'
    val = np.zeros((_Ne, _ny, _nx))
    for i in range(_Ne):
        with netCDF4.Dataset(fileformat.format(i=i+1), 'r') as f:
            val[i] = f[varname][it, 0]
    return  val.mean(axis=0), val.std(axis=0, ddof=1)

def wrap_data(lons):
    fixed_lons = lons.copy()
    for i, start in enumerate(np.argmax(np.abs(np.diff(lons)) > 180, axis=1)):
        fixed_lons[i, start+1:] += 360
    return fixed_lons


def getMask():
    f = netCDF4.Dataset('../../INPUTS/PHYSICS/DOM/eORCA_R1_zps_domcfg.nc', 'r')
    mask = f['top_level'][0] < 1e-3
    f.close()
    return mask


def getCoord():
    f = netCDF4.Dataset('../../INPUTS/PHYSICS/DOM/eORCA_R1_zps_domcfg.nc', 'r')
    lons = wrap_data( f['nav_lon'][:] )
    lats = f['nav_lat'][:]
    f.close()
    return lons, lats


def getModelonObs(exp, year, month):
    var = np.zeros((getMonthDays(year, month), _ny, _nx))
    f = netCDF4.Dataset(f'{expname}/{year}/{month}/eORCA1_1d_20150101_20150131_ptrc1_T_201501-201501.nc', 'r')
    var = f['CHD'][:, 0] + f['CHN'][:, 0]
    f.close()
    return var


def getBackgroundMean(expname, year, month, day):
    ensmean = 0.
    for i in range(1, _nEns+1):
        f = netCDF4.Dataset(f'{expname}/{year}/{month}/state_chlorophyll_{year}{month}{day}_{str(i).zfill(3)}.nc', 'r')
        ensmean += f['chlorophyll'][0]/_nEns
        f.close()
    return ensmean


def getModelSpread(expname, year, month, day):
    ensmean = 0.
    for i in range(1, _nEns+1):
        f = netCDF4.Dataset(f'{expname}/{year}/{month}/state_chlorophyll_{year}{month}{day}_{str(i).zfill(3)}.nc', 'r')
        ensmean += f['chlorophyll'][0]/_nEns
        f.close()
    return ensmean


def getObservation(year, month, day):
    f = netCDF4.Dataset(f'/work/n01/n01/ymchen/obs/chlo_{year}{month}{day}.nc', 'r')
    variable = f['chlorophyll'][:]
    f.close()
    return variable

def getR(year, month, day):
    f = netCDF4.Dataset(f'/work/n01/n01/ymchen/obs/chlo_{year}{month}{day}.nc', 'r')
    variable = f['chlorophyll_unc'][:]
    f.close()
    return variable


def getMonthDays(year, month):
    return calendar.monthrange(int(year), int(month))[1]


def getIncrement(ensMean):
    increment = dict()
    for varname in ensMean:
        increment[varname] = ensMean[varname][1] - ensMean[varname][0]
    return increment


def plotAxes(fig, ax, data, properties):
    ax.set_global()
    pc = ax.pcolormesh(data['x'], data['y'], data['c'], cmap=properties['cmap'],
                       norm=properties['norm'],
                       transform=ccrs.PlateCarree())
    ax.coastlines(color='k', linewidth=.8)
    ax.set_title(properties['title'])
    fig.colorbar(pc, ax=ax, orientation='horizontal', shrink=0.6, pad=0.02)


def plotConc(data, axProperties, figProperty):
    fig = plt.figure()
    naxes = len(data)
    h, w = fig.get_size_inches()
    w = naxes*w
    fig.set_size_inches(w, h)
    gs = matplotlib.gridspec.GridSpec(1, naxes, figure=fig, wspace=0.01, left=0., right=1., bottom=0.0, top=1.)
    fig.clf()
    for ax_gs, ax_data, ax_property in zip(gs, data, axProperties):
        ax = fig.add_subplot(ax_gs, projection=ccrs.Robinson())
        plotAxes(fig, ax, ax_data, ax_property)

    fig.subplots_adjust(left=0.0, bottom=0, right=1., top=1)
    fig.savefig(figProperty['title'], dpi=300)
    plt.close(fig)


def initialisePlot(n):
    data = [dict() for i in range(n)]
    axProperties = [dict() for i in range(n)]
    figProperty = dict()
    return data, axProperties, figProperty

if __name__ == '__main__':
    obsvarnames = ['chlorophyll',]
    varnames =['toce_con', 'soce_abs',]
    year = '2015'
    month = '01'
    ilvl = 0

    lons, lats = getCoord()
    mask = getMask()
    data, axProperties, figProperty = initialisePlot(n=2)

    def f(it, mask):
        
        for i, varname in enumerate(varnames):
            figProperty['title'] = f'spread_SSTSSS_{it}.png'

            # x = np.array([])
            
            _, std = getMeanSTD(it, varname)

            data[i]['x'] = lons
            data[i]['y'] = lats
            data[i]['c'] = np.ma.masked_where(mask, std)
            # x = np.concatenate([np.abs(data[0]['c'])[~mask].ravel(), x])
            axProperties[i]['cmap'] = 'cmo.amp'
            axProperties[i]['title'] = varname
            vmax = np.quantile(np.abs(data[i]['c'])[~mask].ravel(), 0.9)
            axProperties[i]['norm'] = mcolors.Normalize(vmax=vmax, vmin=-vmax)

        plotConc(data, axProperties, figProperty)

    import multiprocessing as mp
    processes = []
    for it in range(0, 11):
        p = mp.Process(target=f, args=(it, mask ))
        p.start()
        processes.append(p)


