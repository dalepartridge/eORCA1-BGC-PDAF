import calendar
import datetime
import os
os.environ['MPLCONFIGDIR'] = "/work/n01/n01/ymchen/config/matplotlib"
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

def getMeanSTD(year, month, it, basefilename):
    val = np.zeros((_Ne, _ny, _nx))
    for i in range(_Ne):
        with netCDF4.Dataset(f'{year}/{month}/ensemble_{i+1}/{basefilename}', 'r') as f:
            val[i] = f['CHD'][it, 0] + f['CHN'][it, 0]
            val[i][val[i] <= 0] = -14
            val[i][val[i] > 0] = np.log10(val[i][val[i] > 0.])
    return  val.mean(axis=0), val.std(axis=0, ddof=1)


def wrap_data(lons):
    fixed_lons = lons.copy()
    for i, start in enumerate(np.argmax(np.abs(np.diff(lons)) > 180, axis=1)):
        fixed_lons[i, start+1:] += 360
    return fixed_lons


def getMask():
    f = netCDF4.Dataset('/work/n01/n01/ymchen/eORCA1-BGC-PDAF/INPUTS/PHYSICS/DOM/eORCA_R1_zps_domcfg.nc', 'r')
    mask = f['top_level'][0] < 1e-3
    f.close()
    return mask


def getCoord():
    f = netCDF4.Dataset('/work/n01/n01/ymchen/eORCA1-BGC-PDAF/INPUTS/PHYSICS/DOM/eORCA_R1_zps_domcfg.nc', 'r')
    lons = wrap_data( f['nav_lon'][:] )
    lats = f['nav_lat'][:]
    f.close()
    return lons, lats


def getEnsMeanSTD(varnames, year, month, day):
    val = dict()
    for varname in varnames:
        val[varname] = np.zeros((_Ne, 2, _ny, _nx))
        for i in range(1, _Ne+1):
            f = netCDF4.Dataset(f'{year}/{month}/state_chlorophyll_{year}{month}{day}_{str(i).zfill(3)}.nc', 'r')
            val[varname][i-1] = f[varname][:]
            f.close()

    mean = dict()
    std = dict()
    for varname in varnames:
        mean[varname] = val[varname].mean(axis=0)
        std[varname] = val[varname].std(axis=0, ddof=1)
    return mean, std


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


# def BackgroundErrorResidual(dob):
#     Edob = dob@dob.T
    

#     - (HBH - R)


def plot2DAxes(fig, ax, data, properties):
    ax.set_global()
    pc = ax.pcolormesh(data['x'], data['y'], data['c'], cmap=properties['cmap'],
                       norm=properties['norm'],
                       transform=ccrs.PlateCarree())
    ax.coastlines(color='k', linewidth=.8)
    ax.set_title(properties['title'])
    fig.colorbar(pc, ax=ax, orientation='horizontal', shrink=0.6, pad=0.02)


def plot2D(data, axProperties, figProperty):
    fig = plt.figure()
    naxes = len(data)
    nx = (naxes - 1)//3 + 1
    ny = min(naxes, 3)
    w, h = fig.get_size_inches()
    w, h = ny*w, nx*h
    fig.set_size_inches(w, h)
    gs = matplotlib.gridspec.GridSpec(nx, ny, figure=fig, wspace=0.01, left=0., right=1., bottom=0.0, top=1.)
    fig.clf()
    for ax_gs, ax_data, ax_property in zip(gs, data, axProperties):
        ax = fig.add_subplot(ax_gs, projection=ccrs.Robinson())
        plot2DAxes(fig, ax, ax_data, ax_property)

    fig.subplots_adjust(left=0.0, bottom=0, right=1., top=1)
    fig.savefig(figProperty['title'], dpi=300)
    plt.close(fig)


def initialisePlot(n):
    data = [dict() for i in range(n)]
    axProperties = [dict() for i in range(n)]
    figProperty = dict()
    return data, axProperties, figProperty


def plotDA(it, mask, title, year, month, obsvarnames):
    lons, lats = getCoord()
    mask = getMask()
    data, axProperties, figProperty = initialisePlot(n=4)

    for varname in obsvarnames:
        figProperty['title'] = title
        obs = getObservation(year, month, str(it+1).zfill(2))
        mean, std = getEnsMeanSTD(obsvarnames, year, month, str(it+1).zfill(2))
        r = getR(year, month, str(it+1).zfill(2))

        d = obs - mean[0]
        mask = mask | (obs > 1e10)
        mask = mask | (mean > 10)

        # innovation
        x = np.array([])
        data[0]['x'] = lons
        data[0]['y'] = lats
        data[0]['c'] = np.ma.masked_where(mask, d)
        x = np.concatenate([np.abs(data[0]['c'])[~mask].ravel(), x])
        axProperties[0]['cmap'] = 'cmo.balance'
        axProperties[0]['title'] = 'innovation'
        vmax = np.quantile(np.abs(data[0]['c'])[~mask].ravel(), 0.9)
        axProperties[0]['norm'] = mcolors.Normalize(vmax=vmax, vmin=-vmax)

        data[1]['x'] = lons
        data[1]['y'] = lats
        data[1]['c'] = np.ma.masked_where(mask, mean[1] - mean[0])
        x = np.concatenate([np.abs(data[1]['c'])[~mask].ravel(), x])
        axProperties[1]['cmap'] = 'cmo.balance'
        axProperties[1]['title'] = 'increments'
        vmax = np.quantile(x, 0.9)
        axProperties[1]['norm'] = mcolors.Normalize(vmax=vmax, vmin=-vmax)

        data[2]['x'] = lons
        data[2]['y'] = lats
        data[2]['c'] = np.ma.masked_where(mask, std)
        axProperties[2]['cmap'] = 'cmo.amp'
        axProperties[2]['title'] = 'spread'
        axProperties[2]['norm'] = mcolors.Normalize(vmax=vmax)

        data[3]['x'] = lons
        data[3]['y'] = lats
        data[3]['c'] = np.ma.masked_where(mask, r)
        axProperties[3]['norm'] =  mcolors.Normalize(vmax=vmax)
        axProperties[3]['cmap'] = 'cmo.amp'
        axProperties[3]['title'] = 'r'

        plot2D(data, axProperties, figProperty)


def saveInnovations(it, mask, year, month, obsvarnames):
    mask = getMask()
    mean, std = getEnsMeanSTD(obsvarnames, year, month, str(it+1).zfill(2))
    for varname in obsvarnames:
        obs = getObservation(year, month, str(it+1).zfill(2))

        d = obs - mean[0]
        mask = mask | (obs > 1e10)
        mask = mask | (mean[0] > 10) | (mean[1] > 10)
        d[mask] = 0.
        dd = d.data
        dd[d.mask] = 0.
        np.save(f'd_{year}{month}{str(it+1).zfill(2)}.npy', dd)


def saveSpread(obsvarnames, year, month, it):
    mask = getMask()
    _, std = getEnsMeanSTD(obsvarnames, year, month, str(it+1).zfill(2))
    spread = dict()
    for varname in obsvarnames:
        spread[varname] = np.zeros(2)
        spread[varname][:] = std[~mask].mean(axis=(-1, -2))
        np.savez(f'spread_{year}{month}{str(it+1).zfill(2)}.npz', *spread)


def plotInnovationTemporalMean(rng):
    lons, lats = getCoord()
    d = 0.
    for it in rng:
        d += np.load(f'd_{year}{month}{str(it+1).zfill(2)}.npy', allow_pickle=True)/len(rng)
    data, axProperties, figProperty = initialisePlot(n=1)
    figProperty['title'] = 'bias.png'
    data[0]['x'] = lons
    data[0]['y'] = lats
    data[0]['c'] = np.ma.masked_where(d == 0, d)
    vmax = np.quantile(np.abs(data[0]['c'][d!=0]).ravel(), 0.9)
    axProperties[0]['norm'] = mcolors.Normalize(vmin=-vmax, vmax=vmax)
    axProperties[0]['cmap'] = 'cmo.balance'
    axProperties[0]['title'] = 'bias'
    plotConc(data, axProperties, figProperty)


def plotSpreadSpatialMean(year, months):
    lons, lats = getCoord()
    std_f = []
    std_a = []
    t = []
    for month in months:
        rng = range(getMonthDays(year, month)) if month != '01' else range(12, getMonthDays(year, month))
        for it in rng:
            t.append(datetime.date(int(year), int(month), it+1))
            std_f.append(
                np.load(f'spread_{year}{month}{str(it+1).zfill(2)}.npz',
                        allow_pickle=True)['chlorophyll'][0]
                )
            std_a.append(
                np.load(f'spread_{year}{month}{str(it+1).zfill(2)}.npz',
                        allow_pickle=True)['chlorophyll'][1]
                )
    fig = plt.figure(1)
    fig.clf()
    ax = fig.add_subplot(111)
    ax.plot(t, std_f, label='forecast spread')
    ax.plot(t, std_a, label='analysis spread')
    ax.legend()
    fig.savefig('spreadSeries.png', dpi=300)
    plt.close(fig)


if __name__ == '__main__':
    obsvarnames = ['chlorophyll',]
    varnames =['CHD', 'CHN',]
    year = '2015'
    month = '01'
    ilvl = 0

    import multiprocessing as mp
    processes = []
    for month in ['01', '02']:
        rng = range(getMonthDays(year, month)) if month != '01' else range(12, getMonthDays(year, month))
        for it in rng:
            p = mp.Process(target=saveSpread, args=(obsvarnames, year, month, it))
            p.start()
            processes.append(p)

    for p in processes:
        p.join()

    # plotSpreadSpatialMean(year, ['01', '02'])