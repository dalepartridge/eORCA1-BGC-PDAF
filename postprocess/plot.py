import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.gridspec
import numpy as np
import cmocean
import cartopy.crs as ccrs


_vmin = {'spread': {'PHD': 0. , 'PHN': 0.,  'CHD': 0., 'CHN':0.0},
         'increment_spread': {'PHD': -0.6, 'PHN': -0.12, 'CHD': -0.1, 'CHN':-0.1},
         'state':  {'PHD': 0. , 'PHN': 0. , 'CHD': 0. , 'CHN': 0. },
         'increment_state' : {'PHD': -1.3, 'PHN': -1.1, 'CHD': -1.5, 'CHN': -3 },
         }
_vmax = {'spread': {'PHD': 0.9, 'PHN': 0.5,  'CHD': 0.1, 'CHN':0.1},
         'increment_spread': {'PHD':  0.6, 'PHN':  0.12, 'CHD': 0.1, 'CHN':0.1},
         'state':  {'PHD': 1.5*2, 'PHN': 1.5*4, 'CHD': 0.7, 'CHN':0.7},
         'increment_state' : {'PHD':  1.3, 'PHN':  1.1*2, 'CHD': 0.7, 'CHN': 0.7},
         }
def plot(f, lons, lats, mask, varname, date, typename):
    var = f[varname]
    cmap = 'cmo.algae'
    if typename == 'variance':
        typename = 'spread'
        var = np.sqrt(var)
        cmap = 'cmo.amp'

    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']
    plt.rcParams['font.size'] = 12

    print (f'start plotting {typename} {varname} on {date}')
    fig = plt.figure(1, figsize=(3*6.4, 4.8))
    gs = matplotlib.gridspec.GridSpec(1, 3, wspace=0.01, left=0., right=1., bottom=0.0, top=1.)
    fig.clf()

    vmax = _vmax[typename][varname] # np.nanquantile(var[:, 0], 0.99)
    vmin = _vmin[typename][varname] # max(np.nanquantile(var[:, 0], 0.2), 0.)

    ax = fig.add_subplot(gs[0], projection=ccrs.Robinson())
    ax.set_global()
    data = var.isel(t=0, z=0)
    pc = ax.pcolormesh(lons, lats, np.ma.masked_where(np.isnan(data) | mask, data),
                      cmap=cmap, vmin=vmin, vmax=vmax , transform=ccrs.PlateCarree())
    ax.coastlines(color='k', linewidth=.8)
    ax.set_title(f'{varname} forecast {typename}')
    fig.colorbar(pc, ax=ax, orientation='horizontal', shrink=0.6, pad=0.02)

    ax = fig.add_subplot(gs[1], projection=ccrs.Robinson())
    ax.set_global()
    data = var.isel(t=1, z=0)
    pc = ax.pcolormesh(lons, lats, np.ma.masked_where(np.isnan(data) | mask, data),
                       cmap=cmap, vmin=vmin, vmax=vmax , transform=ccrs.PlateCarree())
    ax.coastlines(color='k', linewidth=.8)
    ax.set_title(f'{varname} analysis {typename}')
    fig.colorbar(pc, ax=ax, orientation='horizontal', shrink=0.6, pad=0.02)

    ax = fig.add_subplot(gs[2], projection=ccrs.Robinson())
    ax.set_global()
    ax.set_title(f'{typename} increment')
    data = var[1, 0] - var[0, 0]
    vmax = _vmax['increment_'+typename][varname] # np.nanquantile(var[:, 0], 0.99)
    offset = mcolors.TwoSlopeNorm(vcenter=0., vmin=-1.*vmax, vmax=vmax)
    pc = ax.pcolormesh(lons, lats, np.ma.masked_where(np.isnan(data) | mask, data),
                       cmap='cmo.balance', norm=offset ,  transform=ccrs.PlateCarree())
    ax.coastlines(color='k', linewidth=.8)
    fig.colorbar(pc, ax=ax, orientation='horizontal', shrink=0.6, pad=0.02)
    fig.subplots_adjust(left=0.0, bottom=0, right=1., top=1)

    print (f'saving {typename} {varname} on {date}')
    plt.savefig(f'{typename}_{varname}_{date}.png', dpi=300)


def wrap_data(lons):
    fixed_lons = lons.copy()
    for i, start in enumerate(np.argmax(np.abs(np.diff(lons)) > 180, axis=1)):
        fixed_lons[i, start+1:] += 360
    return fixed_lons


def plotObs(varname, lons, lats, mask, obs, forecast, date):
    cmap = 'cmo.algae'

    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']
    plt.rcParams['font.size'] = 12

    print (f'start plotting innovation {varname} on {date}')
    fig = plt.figure(1, figsize=(3*6.4, 4.8))
    gs = matplotlib.gridspec.GridSpec(1, 3, wspace=0.01, left=0., right=1., bottom=0.0, top=1.)
    fig.clf()

    vmax =  _vmax['state'][varname] #  np.nanquantile(obs[varname], 0.99)
    vmin =  _vmin['state'][varname] # max(np.nanquantile(obs[varname], 0.2), 0.)

    ax = fig.add_subplot(gs[0], projection=ccrs.Robinson())
    # pc = ax.tripcolor(lons, lats, obs[varname].ravel()[~isnan],
    data = obs[varname]
    pc = ax.pcolormesh(lons, lats, np.ma.masked_where(np.isnan(data) | mask, data),
                      cmap=cmap, vmin=vmin, vmax=vmax, transform=ccrs.PlateCarree())
    ax.coastlines()
    ax.set_title(f'{varname} observation')
    fig.colorbar(pc, ax=ax, orientation='horizontal', shrink=0.6, pad=0.02)

    ax = fig.add_subplot(gs[1], projection=ccrs.Robinson())
    # pc = ax.tripcolor(lons, lats, forecast[varname].ravel()[~isnan],
    data = forecast[varname]
    pc = ax.pcolormesh(lons, lats, np.ma.masked_where(np.isnan(data) | mask, data),
                      cmap=cmap, vmin=vmin, vmax=vmax, transform=ccrs.PlateCarree())
    ax.coastlines()
    ax.set_title(f'{varname} forecast')
    fig.colorbar(pc, ax=ax, orientation='horizontal', shrink=0.6, pad=0.02)

    ax = fig.add_subplot(gs[2], projection=ccrs.Robinson())
    ax.coastlines()
    ax.set_title(f'{varname} innovation')
    vmax = _vmax['increment_state'][varname]
    offset = mcolors.TwoSlopeNorm(vcenter=0., vmin=-vmax, vmax=vmax)
    data = obs[varname] - forecast[varname]
    pc = ax.pcolormesh(lons, lats, np.ma.masked_where(np.isnan(data) | mask, data),
                      cmap='cmo.balance', norm=offset, transform=ccrs.PlateCarree())
    fig.colorbar(pc, ax=ax, orientation='horizontal', shrink=0.6, pad=0.02)
    fig.subplots_adjust(left=0.0, bottom=0, right=1., top=1)

    print (f'start saving innovation {varname} on {date}')
    fig.savefig(f'obs_{varname}_{date}.png', dpi=300)


def getPdafState(day, month, year, cycle_month, level, expname, varnames):
    f = xr.open_dataset(f'{expname}/{year}/{cycle_month}/state_{year}{month}{day}.nc')
    ny, nx = f.dims['y'], f.dims['x']
    forecast = dict()
    analysis = dict()
    for varname in varnames:
        forecast[varname] = np.zeros((ny, nx))
        analysis[varname] = np.zeros((ny, nx))

        forecast[varname] = np.squeeze(f[varname][0, level].to_numpy())
        analysis[varname] = np.squeeze(f[varname][1, level].to_numpy())

    lons = wrap_data(np.ma.masked_where(f['nav_lon'].data > 1000., f['nav_lon'].data))
    lats = np.ma.masked_where(f['nav_lat'].data > 1000., f['nav_lat'].data)
    f.close()
    return lons, lats, forecast, analysis


def getPdafVariance(day, month, year, cycle_month, level, expname, varnames):
    f = xr.open_dataset(f'{expname}/{year}/{cycle_month}/variance_{year}{month}{day}.nc')
    ny, nx = f.dims['y'], f.dims['x']
    forecast = dict()
    analysis = dict()
    for varname in varnames:
        forecast[varname] = np.zeros((ny, nx))
        analysis[varname] = np.zeros((ny, nx))

        forecast[varname] = np.squeeze(f[varname][0, level].to_numpy())
        analysis[varname] = np.squeeze(f[varname][1, level].to_numpy())

    lons = wrap_data(np.ma.masked_where(f['nav_lon'].data > 1000., f['nav_lon'].data))
    lats = np.ma.masked_where(f['nav_lat'].data > 1000., f['nav_lat'].data)
    f.close()
    return lons, lats, forecast, analysis


def plotInnovation(year, month, day, cycle_month, expname, varnames,
                   obstype):
    lons, lats, forecast, _ = getPdafState(day, month, year, cycle_month, 0, expname, varnames)
    ny, nx = lons.shape

    # f = xr.open_dataset(f'processed_obs/{obstype}_{year}{month}.nc')
    f = xr.open_dataset(f'{obstype}_{year}{month}.nc')
    obs = dict()
    for varname in varnames:
        obs[varname] = f[varname].to_numpy()
        obs[varname][f[varname] > 1e10] = np.nan
    f.close()

    fm = xr.open_dataset('../INPUTS/PHYSICS/DOM/eORCA_R1_zps_domcfg.nc')
    mask = fm['top_level'][0].to_numpy() < 1e-3
    fm.close()

    for varname in varnames:
        plotObs(varname, lons, lats, mask, obs, forecast, f'{year}{month}')


def plotStateIncrement(year, month, day, typename, cycle_month, expname, varnames):
    f = xr.open_dataset(f'{expname}/{year}/{cycle_month}/{typename}_{year}{month}{day}.nc')
    lons = wrap_data(np.ma.masked_where(f['nav_lon'].data > 1000., f['nav_lon'].data))
    lats = np.ma.masked_where(f['nav_lat'].data > 1000., f['nav_lat'].data)
    fm = xr.open_dataset('../INPUTS/PHYSICS/DOM/eORCA_R1_zps_domcfg.nc')
    mask = fm['top_level'][0].to_numpy() < 1e-3
    fm.close()

    for varname in varnames:
        plot(f, lons, lats, mask, varname, f'{year}{month}', typename)


def plotInit(expname, varname):
    f = xr.open_dataset(f'{expname}/2015/01/state_20150101_ini.nc')
    frestart = xr.open_dataset('/work/n01/n01/ymchen/eORCA1-BGC-PDAF/INPUTS/MEDUSA/DOM/restart_trc.nc')

    lons = wrap_data(np.ma.masked_where(f['nav_lon'].data > 1000., f['nav_lon'].data))
    lats = np.ma.masked_where(f['nav_lat'].data > 1000., f['nav_lat'].data)

    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']
    plt.rcParams['font.size'] = 12

    fig = plt.figure(1, figsize=(3*6.4, 4.8))
    gs = matplotlib.gridspec.GridSpec(1, 3, wspace=0.01, left=0., right=1., bottom=0.0, top=1.)
    fig.clf()
    cmap = 'cmo.algae'

    ax = fig.add_subplot(gs[0], projection=ccrs.Robinson())
    var = f[varname]
    vmax = np.nanquantile(var[0, 0], 0.99)
    vmin = max(np.nanquantile(var[0, 0], 0.2), 0.)
    data = var.isel(t=0, z=0)
    pc = ax.pcolormesh(lons, lats, np.ma.masked_where(np.isnan(data) | mask, data),
                       cmap=cmap, vmin=vmin, vmax=vmax , transform=ccrs.PlateCarree())
    ax.coastlines()
    ax.set_title(f'{varname} init')
    fig.colorbar(pc, ax=ax, orientation='horizontal', shrink=0.6, pad=0.02)

    ax = fig.add_subplot(gs[1], projection=ccrs.Robinson())
    var = frestart[f'TRNpelagic_{varname}']
    var = var.where(var > 1e-5)
    vmax = np.nanquantile(var[0, 0], 0.99)
    vmin = max(np.nanquantile(var[0, 0], 0.2), 0.)
    data = var.isel(t=0, z=0)
    pc = ax.pcolormesh(lons, lats, np.ma.masked_where(np.isnan(data) | mask, data),
                       cmap=cmap, vmin=vmin, vmax=vmax , transform=ccrs.PlateCarree())
    ax.coastlines()
    ax.set_title(f'{varname} restart')
    fig.colorbar(pc, ax=ax, orientation='horizontal', shrink=0.6, pad=0.02)

    fobs = xr.open_dataset('../../obs/BICEP_NCEO_PC_ESA-OC-L3S-MERGED-1M_MONTHLY_9km_mapped_201501-fv5.nc')
    var = fobs['C_microphyto']*2/3/6.625
    var = var.where(var > 1e-5)

    lons = fobs['longitude'].to_numpy()
    lats = fobs['latitude'].to_numpy()

    ax = fig.add_subplot(gs[2], projection=ccrs.Robinson())
    vmax = np.nanquantile(var[0], 0.99)
    vmin = max(np.nanquantile(var[0], 0.2), 0.)
    pc = ax.pcolormesh(lons, lats, var[0].to_numpy(),
                       cmap=cmap, vmin=vmin, vmax=vmax, transform=ccrs.PlateCarree())
    ax.coastlines()
    ax.set_title(f'PHD obs')
    fig.colorbar(pc, ax=ax, orientation='horizontal', shrink=0.6, pad=0.02)

    fig.savefig('init.png', dpi=300)


def selectPosition(filename, varname):
    f = xr.open_dataset(filename)
    ny, nx = f.dims['y'], f.dims['x']
    incr = (f[varname][1, 0] - f[varname][0, 0]).to_numpy().ravel()
    i = np.nanargmax(np.abs(incr))
    mj = i//nx
    mi = i - mj*nx
    lons = f['nav_lon'].to_numpy().ravel()
    lats = f['nav_lat'].to_numpy().ravel()
    f.close()
    return mi, mj, lons, lats


def plotTimeSeries(selectionfilename, expname, obstype):
    # choose position
    i, mi, mj, lons, lats = selectPosition(selectionfilename, 'PHD')

    forvals = dict(PHD=[], PHN=[])
    anavals = dict(PHD=[], PHN=[])
    obsvals = dict(PHD=[], PHN=[])
    for month, day, dirname in zip(['01', '02', '03', '04', '05', '06'],
                                   ['16', '14', '15', '15', '15', '15'],
                                   ['01', '01', '01', '04', '04', '04']):
        obsvar = np.load(f'processed_obs/{obstype}_{year}{month}.npz')
        _, _, forecast, analysis = getPdafState(day, month, '2015', dirname, 0, expname)
        for varname in obsvar:
            forvals[varname].append(forecast[varname].ravel()[i])
            anavals[varname].append(analysis[varname].ravel()[i])
            obsvals[varname].append(obsvar[varname].ravel()[i])

    for varname in obsvals:
        fig = plt.figure(1, figsize=(6.4, 4.8))
        gs = matplotlib.gridspec.GridSpec(1, 1)
        fig.clf()
        ax = fig.add_subplot(gs[0])
        ax.plot(range(1, 7), obsvals[varname], color='k', label='observation')
        ax.plot(range(1, 7), forvals[varname], color='r', label='forecast')
        ax.plot(range(1, 7), anavals[varname], color='b', label='analysis')
        ax.legend()
        ax.set_title(f'time series of PHD at ({lons[i]}, {lats[i]})')
        fig.savefig(f'timeseries_Jan_{varname}.png', dpi=300)


def compareSTD(day, month, year, cycle_month, expname,
               obstype):
    lons, lats, forecast, _ = getPdafVariance(day, month, year, cycle_month, 0)
    forecast = {varname:np.sqrt(forecast[varname]).ravel() for varname in forecast}

    fo = np.load(f'processed_obs/{obstype}_{year}{month}.npz')
    obs = {varname:0.3*obs[varname] for varname in obs}
    fo.close()

    fm = xr.open_dataset('../INPUTS/PHYSICS/DOM/eORCA_R1_zps_domcfg.nc')
    mask = fm['top_level'][0].to_numpy() < 1e-3
    fm.close()

    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']
    plt.rcParams['font.size'] = 12

    for varname in forecast:
        fig = plt.figure(1, figsize=(2*6.4, 4.8))
        gs = matplotlib.gridspec.GridSpec(1, 2, wspace=0.01, left=0., right=1., bottom=0.0, top=1.)
        fig.clf()

        cmap = 'cmo.amp'

        ax = fig.add_subplot(gs[0], projection=ccrs.Robinson())
        vmax = np.nanquantile(forecast[varname].ravel(), 0.99)
        vmin = max(np.nanquantile(forecast[varname].ravel(), 0.2), 0.)
        data = forecast[varname]
        pc = ax.pcolormesh(lons, lats, np.ma.masked_where(np.isnan(data) | mask, data),
                        cmap=cmap, vmin=vmin, vmax=vmax, transform=ccrs.PlateCarree())
        ax.coastlines()
        ax.set_title(f'{varname} forecast standard deviation')
        fig.colorbar(pc, ax=ax, orientation='horizontal', shrink=0.6, pad=0.02)

        ax = fig.add_subplot(gs[1], projection=ccrs.Robinson())
        vmax = np.nanquantile(obs[varname].ravel(), 0.99)
        vmin = max(np.nanquantile(obs[varname].ravel(), 0.2), 0.)
        data = obs[varname]
        pc = ax.pcolormesh(lons, lats, np.ma.masked_where(np.isnan(data) | mask, data),
                        cmap=cmap, vmin=vmin, vmax=vmax, transform=ccrs.PlateCarree())
        ax.coastlines()
        ax.set_title(f'{varname} observation standard deviation')
        fig.colorbar(pc, ax=ax, orientation='horizontal', shrink=0.6, pad=0.02)

        plt.savefig(f'std_{varname}_{year}{month}{day}.png', dpi=300)


def plotSTDInnovation(varname, lons, lats, isnan, obs, forecast, forecast_var, date):
    cmap = 'cmo.algae'

    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']
    plt.rcParams['font.size'] = 12

    fig = plt.figure(1, figsize=(3*6.4, 4.8))
    gs = matplotlib.gridspec.GridSpec(1, 3, wspace=0.01, left=0., right=1., bottom=0.0, top=1.)
    fig.clf()
    vmax =  0.3*_vmax['state'][varname] #  np.nanquantile(obs[varname], 0.99)
    vmin =  0.3*_vmin['state'][varname] # max(np.nanquantile(obs[varname], 0.2), 0.)
    ax = fig.add_subplot(gs[0], projection=ccrs.Robinson())
    pc = ax.tripcolor(lons, lats, 0.3*obs[varname].ravel()[~isnan],
                      cmap=cmap, vmin=vmin, vmax=vmax, transform=ccrs.PlateCarree())
    ax.coastlines()
    ax.set_title(f'{varname} observation std')
    fig.colorbar(pc, ax=ax, orientation='horizontal', shrink=0.6, pad=0.02)

    ax = fig.add_subplot(gs[1], projection=ccrs.Robinson())
    pc = ax.tripcolor(lons, lats, np.sqrt(forecast_var[varname].ravel()[~isnan]),
                      cmap=cmap, vmin=vmin, vmax=vmax, transform=ccrs.PlateCarree())
    ax.coastlines()
    ax.set_title(f'{varname} forecast std')
    fig.colorbar(pc, ax=ax, orientation='horizontal', shrink=0.6, pad=0.02)

    ax = fig.add_subplot(gs[2], projection=ccrs.Robinson())
    ax.coastlines()
    ax.set_title(f'{varname} variance difference (y - x)')
    innovation = 0.3*obs[varname] - forecast_var[varname]
    vmax = np.nanquantile(obs[varname], 0.7) # _vmax['increment_state'][varname]
    offset = mcolors.TwoSlopeNorm(vcenter=0., vmin=-vmax, vmax=vmax)
    pc = ax.tripcolor(lons, lats, innovation.ravel()[~isnan],
                      cmap='cmo.balance', norm=offset, transform=ccrs.PlateCarree())
    fig.colorbar(pc, ax=ax, orientation='horizontal', shrink=0.6, pad=0.02)
    fig.subplots_adjust(left=0.0, bottom=0, right=1., top=1)

    fig.savefig(f'var_{varname}_{date}.png', dpi=300)


if __name__ == '__main__':
    varnames =['CHD', 'CHN']
    obstype = 'chlo'
    expname = 'exp-logmean-chlo'
    assim_day = {'01':'16', '02':'14', '03':'15',
                 '04':'15', '05':'15', '06':'15',
                 '07':'16', '08':'15', '09':'15',
                 '10':'16', '11':'15', '12':'15',}
    import multiprocessing as mp
    print (mp.cpu_count())
    processes = []
    # for ic in ('01', '04', '07', '10'):
    #     for month in range(int(ic), int(ic)+3):
    #         print (ic, month)
    #         for typename in ['state', 'variance']:
    #             p = mp.Process(target=plotStateIncrement,
    #                            args=('2015', str(month).zfill(2),
    #                                  assim_day[str(month).zfill(2)],
    #                                  typename, ic, expname, varnames))
    #             p.start()
    #             processes.append(p)

    # for ic in ('01', '04', '07', '10'):
    #     for month in range(int(ic), int(ic)+3):
    ic = '01'
    month = 1
    print (ic, month)
    p = mp.Process(target=plotInnovation,
                   args=('2015', str(month).zfill(2),
                         assim_day[str(month).zfill(2)],
                         ic, expname, varnames, obstype))

    p.start()
    processes.append(p)

    for p in processes:
        p.join()


    ## plotInit('exp-logmean')
    #plotTimeSeries('exp-logmean/2015/01/state_20150116.nc', 'exp-logmean')
    ## plotTimeSeries('exp-logmean/2015/07/state_20150716.nc', 'exp-logmean')
    #days = ['16', '14', '15', '15', '15', '15']
    #months = ['01', '02', '03', '04', '05', '06']
    #years = ['2015']*6
    #cycle_months = ['01', '01', '01', '04', '04', '04']
    #for day, month, year, cycle_month in zip(days, months, years, cycle_month):
    #    print (day, month, year)
    #    compareSTD(day, month, year, cycle_month, 'exp-logmean')
