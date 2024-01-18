import calendar
import datetime
import cartopy.crs as ccrs
import cmocean
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.gridspec
import numpy as np
import netCDF4

plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']
plt.rcParams['font.size'] = 12
_ny = 332
_nx = 362
_nEns = 30


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


def ModelOutputTreatmentBGC(val, i, it, *args):
    val[i] = 0.
    for var in args:
        val[i] += var[it, 0]
    val[i][val[i] <= 0] = -14
    val[i][val[i] > 0] = np.log10(val[i][val[i] > 0.])
    return val

def ModelOutputTreatmentPhys(val, i, it, *args):
    val[i] = args[0][it, 0]
    return val


def getModelOutputSTD(year, month, it, varinfo):
    val = np.zeros((_nEns, _ny, _nx))
    for i in range(_nEns):
        with netCDF4.Dataset(f'{year}/{month}/ensemble_{i+1}/{varinfo["filename"]}', 'r') as f:
            val = varinfo['op'](val, i, it, *(f[varname] for varname in varinfo['varname']))
    return  val.mean(axis=0), val.std(axis=0, ddof=1)


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


def getMonthDays(year, month):
    return calendar.monthrange(int(year), int(month))[1]


def getIncrement(ensMean):
    increment = dict()
    for varname in ensMean:
        increment[varname] = ensMean[varname][1] - ensMean[varname][0]
    return increment


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


def initialisePlot(expnames):
    data = [dict() for i in range(len(expnames))]
    axProperties = [dict() for i in range(len(expnames))]
    figProperty = dict()
    return data, axProperties, figProperty


def plotFreeRunSTDOverTime():
    import multiprocessing as mp
    print (mp.cpu_count())
    lons, lats = getCoord()
    mask = getMask()
    varinfo = dict()
    varinfo['chlorophyll'] = dict(varname=['CHN', 'CHD'],
                                  op=ModelOutputTreatmentBGC,
                                  gridType='ptrc1_T')
    varinfo['nitrogen'] = dict(varname=['PHN', 'PHD' ],
                          op=ModelOutputTreatmentBGC,
                          gridType='ptrc1_T')
    varinfo['SST'] = dict(varname=['toce_con', ],
                          op=ModelOutputTreatmentPhys,
                          gridType='grid_T')
    varinfo['SSS'] = dict(varname=['soce_abs', ],
                          op=ModelOutputTreatmentPhys,
                          gridType='grid_T')
    data, axProperties, figProperty = initialisePlot(varinfo)
    def f(it, month, mask, cmonth, filerange):
        figProperty['title'] = f'fcst_{month}{str(it).zfill(2)}.png'
        for i, varname in enumerate(varinfo):
            varinfo[varname]['filename'] = f'eORCA1_1d_{filerange}_{varinfo[varname]["gridType"]}_2015{month}-2015{month}.nc'
            print (month, it, varname, varinfo[varname]['filename'])
            meanVal, stdVal = getModelOutputSTD('2015', cmonth, it, varinfo[varname])
            mask = mask | (meanVal > 1e10)
            data[i]['x'] = lons
            data[i]['y'] = lats
            data[i]['c'] = np.ma.masked_where(mask, stdVal)
            axProperties[i]['cmap'] = 'cmo.amp'
            axProperties[i]['title'] = varname
            axProperties[i]['norm'] = mcolors.Normalize(vmax=np.quantile(data[i]['c'], 0.9))
        plot2D(data, axProperties, figProperty)

    import multiprocessing as mp
    processes = []
    for month, cmonth in zip(['01', '02', '03', '04'], 
                             ['01', '01', '03', '03'],
                             ['20150101_20150228', '20150101_20150228',
                              '20150301_20150430', '20150301_20150430' ]):
        rng = range(getMonthDays(year, month)) if month != '01' else range(12, getMonthDays(year, month))
        for it in rng:
            p = mp.Process(target=f, args=(it, month, mask, cmonth))
            p.start()
            processes.append(p)

        for p in processes:
            p.join()


def saveSpread():
    import multiprocessing as mp
    print (mp.cpu_count())
    lons, lats = getCoord()
    mask = getMask()
    varinfo = dict()
    varinfo['chlorophyll'] = dict(varname=['CHN', 'CHD'],
                                  op=ModelOutputTreatmentBGC,
                                  gridType='ptrc1_T')
    varinfo['nitrogen'] = dict(varname=['PHN', 'PHD' ],
                          op=ModelOutputTreatmentBGC,
                          gridType='ptrc1_T')
    varinfo['SST'] = dict(varname=['toce_con', ],
                          op=ModelOutputTreatmentPhys,
                          gridType='grid_T')
    varinfo['SSS'] = dict(varname=['soce_abs', ],
                          op=ModelOutputTreatmentPhys,
                          gridType='grid_T')
    def f(it, month, mask, cmonth, filerange):
        std = dict()
        for i, varname in enumerate(varinfo):
            varinfo[varname]['filename'] = f'eORCA1_1d_{filerange}_{varinfo[varname]["gridType"]}_2015{month}-2015{month}.nc'
            print (month, it, varname, varinfo[varname]['filename'])
            meanVal, stdVal = getModelOutputSTD('2015', cmonth, it, varinfo[varname])
            mask = mask | (meanVal > 1e10)
            std[varname] = np.mean(np.ma.masked_where(mask, stdVal))
        np.savez(f'meanstd_{month}{str(it).zfill(2)}.npz', **std)

    import multiprocessing as mp
    processes = []
    year = '2015'
    for month, cmonth, filerange in zip(['01', '02', '03', '04'], 
                             ['01', '01', '03', '03'],
                             ['20150101_20150228', '20150101_20150228',
                              '20150301_20150430', '20150301_20150430' ]):
        rng = range(getMonthDays(year, month)) if month != '01' else range(12, getMonthDays(year, month))
        for it in rng:
            p = mp.Process(target=f, args=(it, month, mask, cmonth, filerange))
            p.start()
            processes.append(p)

        for p in processes:
            p.join()


def plotSpreadSpatialMean():
    std = dict()
    std['chlorophyll'] = []
    std['nitrogen'] = []
    std['SSS'] = []
    std['SST'] = []
    t = []
    year = '2015'
    for month in ['01', '02', '03', '04']:
        rng = range(getMonthDays(year, month)) if month != '01' else range(12, getMonthDays(year, month))
        for it in rng:
            t.append(datetime.date(int(year), int(month), it+1))
            for varname in std:
                std[varname].append(
                    np.load(f'meanstd_{month}{str(it).zfill(2)}.npz',
                            allow_pickle=True)[varname]
                    )
    fig = plt.figure(1)
    fig.clf()
    ax = fig.add_subplot(111)
    for varname in std:
        ax.plot(t, std[varname], label=varname)
    for label in ax.get_xticklabels(which='major'):
        label.set(rotation=45, horizontalalignment='right')
    ax.legend()
    fig.savefig('spread.png', dpi=300)
    plt.close(fig)


if __name__ == '__main__':
    # plotFreeRunSTDOverTime()
    # saveSpread()
    plotSpreadSpatialMean()
    # obsvarnames = ['chlorophyll',]
    # varnames =['CHD', 'CHN', 'PHD', 'PHN', 'PDS']
    # year = '2015'
    # month = '01'
    # ilvl = 0

    # lons, lats = getCoord()
    # mask = getMask()
    # data, axProperties, figProperty = initialisePlot(expnames + ['innovation',])

    # def f(it, mask):
    #     for varname in obsvarnames:
    #         print (it, varname)
    #         figProperty['title'] = f'innovation_{varname}_{it}.png'
    #         obs = getObservation(year, month, str(it+1).zfill(2))
    #         background = getBackgroundMean('IAU_bn', year, month, str(it+1).zfill(2))

    #         mask = mask | (obs > 1e10)
    #         mask = mask | (background > 1e10)
    #         x = np.array([])
    #         data[0]['x'] = lons
    #         data[0]['y'] = lats
    #         data[0]['c'] = np.ma.masked_where(mask, background)
    #         x = np.concatenate([data[0]['c'][~mask].ravel(), x])
    #         axProperties[0]['cmap'] = 'cmo.algae'
    #         axProperties[0]['title'] = 'IAU_a'

    #         data[1]['x'] = lons
    #         data[1]['y'] = lats
    #         data[1]['c'] = np.ma.masked_where(mask, obs)
    #         x = np.concatenate([data[1]['c'][~mask].ravel(), x])
    #         axProperties[1]['cmap'] = 'cmo.algae'
    #         axProperties[1]['title'] = 'obs'

    #         print (np.quantile(x, 0.9))
    #         axProperties[0]['norm'] = mcolors.Normalize(vmax=np.quantile(x, 0.9))
    #         axProperties[1]['norm'] = mcolors.Normalize(vmax=np.quantile(x, 0.9))

    #         data[2]['x'] = lons
    #         data[2]['y'] = lats
    #         data[2]['c'] = np.ma.masked_where(mask, obs - background)
    #         vmax = np.quantile(data[2]['c'][~mask].ravel(), 0.9)
    #         print (vmax)
    #         axProperties[2]['norm'] =  mcolors.TwoSlopeNorm(vcenter=0., vmin=-1.*vmax, vmax=vmax)
    #         axProperties[2]['cmap'] = 'cmo.balance'
    #         axProperties[2]['title'] = 'innovation'

    #         plotConc(data, axProperties, figProperty)

    # import multiprocessing as mp
    # processes = []
    # for it in range(1, 15):
    #     p = mp.Process(target=f, args=(it, mask ))
    #     p.start()
    #     processes.append(p)

    # for p in processes:
    #     p.join()

    # for it in range(15, 31):
    #     p = mp.Process(target=f, args=(it, mask ))
    #     p.start()
    #     processes.append(p)

    # for p in processes:
    #     p.join()

