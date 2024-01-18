import calendar

import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.gridspec
import numpy as np
import cmocean
import cartopy
import cartopy.crs as ccrs


_Ne = 30
_range = {'01': ('0101', '0331'), '04':('0401', '0630'),
          '07': ('0701', '0930'), '10':('1001', '1231')}
_vmin = {'spread': {'PHD': 0. , 'PHN': 0. , 'CHN': 0, 'CHD': 0.},
         'increment_spread': {'PHD': -0.6, 'PHN': -0.12},
         'state':  {'PHD': 0. , 'PHN': 0. , 'CHN': 0, 'CHD': 0.},
         'increment_state' : {'PHD': -0.1, 'PHN': -0.13,
                              'CHD': -0.08, 'CHN': -0.1 }
         }
_vmax = {'spread': {'PHD': 0.9, 'PHN': 0.5},
         'increment_spread': {'PHD':  0.6, 'PHN':  0.12},
         'state':  {'PHD': 0.45 , 'PHN': 0.5 , 'CHN': 0.5, 'CHD': 0.45},
         'increment_state' : {'PHD':  0.1, 'PHN':  0.13,
                              'CHD':  0.08, 'CHN':  0.1 }
         }
class ForecastFile:
    def __init__(self, year, cycle_month, month, expname, plotExpname):
        filename = f'eORCA1_1m_{year}{_range[cycle_month][0]}_{year}' \
                    f'{_range[cycle_month][1]}_ptrc_T_{year}{month}-{year}{month}.nc'
        print (filename)
        self.f = [xr.open_dataset(f'{expname}/{year}/{cycle_month}/ensemble_{i}/{filename}') for i in range(1, _Ne+1)]
        self.lons = self.wrap_lons(np.ma.masked_where(self.f[0]['nav_lon'].data > 1000., self.f[0]['nav_lon'].data))
        self.lats = np.ma.masked_where(self.f[0]['nav_lat'].data > 1000., self.f[0]['nav_lat'].data)
        self.month = int(month)
        self.getMask()
        self.expname = plotExpname

    def wrap_lons(self, lons):
        fixed_lons = lons.copy()
        for i, start in enumerate(np.argmax(np.abs(np.diff(lons)) > 180, axis=1)):
            fixed_lons[i, start+1:] += 360
        return fixed_lons

    def getMask(self):
        f = xr.open_dataset('../INPUTS/PHYSICS/DOM/eORCA_R1_zps_domcfg.nc')
        self.mask = f['top_level'][0].to_numpy() < 1e-3
        f.close()

    def checkBounds(self):
        chd = self.f['CHD'].where(self.f['CHD']>-1e10)
        phd = self.f['PHD'].where(self.f['CHD']>-1e10)
        print ('CHD', np.max(chd))
        print ('CHD_q', np.nanquantile(chd, 0.99))
        print ('PHD', np.max(phd))
        print ('PHD_q', np.nanquantile(phd, 0.99))

        chn = self.f['CHN'].where(self.f['CHN']>-1e10)
        phn = self.f['PHN'].where(self.f['CHN']>-1e10)
        print ('CHN', np.max(chn))
        print ('CHN_q', np.nanquantile(chn, 0.99))
        print ('PHN', np.max(phn))
        print ('PHN_q', np.nanquantile(phn, 0.99))

    def checkDiffBounds(self, f):
        chd = self.f['CHD'].where(self.f['CHD']>-1e10)
        phd = self.f['PHD'].where(self.f['CHD']>-1e10)
        chd_f = f.f['CHD'].where(f.f['CHD']>-1e10)
        phd_f = f.f['PHD'].where(f.f['CHD']>-1e10)
        print ('CHD_q', np.nanquantile(np.abs(chd - chd_f), 0.99))
        print ('PHD_q', np.nanquantile(np.abs(phd - phd_f), 0.99))

        chn = self.f['CHN'].where(self.f['CHN']>-1e10)
        phn = self.f['PHN'].where(self.f['CHN']>-1e10)
        chn_f = f.f['CHN'].where(f.f['CHN']>-1e10)
        phn_f = f.f['PHN'].where(f.f['CHN']>-1e10)
        print ('CHN_q', np.nanquantile(np.abs(chn-chn_f), 0.99))
        print ('PHN_q', np.nanquantile(np.abs(phn-phn_f), 0.99))

    def getEnsMean(self, varname):
        print (f'...Computing the mean of {varname}...')
        print ((varname, self.f[0]))
        return np.mean([self.f[i][varname][0, 0] for i in range(_Ne)], axis=0)

    def plot(self, varname1, varname2, typename):
        var1 = self.getEnsMean(varname1)
        var2 = self.getEnsMean(varname2)
        var1 = np.ma.masked_where(np.isnan(var1) | self.mask, var1)
        var2 = np.ma.masked_where(np.isnan(var2) | self.mask, var2)

        lons = self.lons
        lats = self.lats

        cmap = 'cmo.amp' if typename == 'std' else 'cmo.algae'

        plt.rcParams['font.family'] = 'serif'
        plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']
        plt.rcParams['font.size'] = 12

        fig = plt.figure(1)
        w, h = fig.get_size_inches()
        w = 3*w
        h = h
        fig.set_size_inches(w, h)
        gs = matplotlib.gridspec.GridSpec(1, 3, figure=fig, wspace=0.01, left=0., right=1., bottom=0.0, top=1.)
        fig.clf()

        ax = fig.add_subplot(gs[0], projection=ccrs.Robinson())
        ax.set_title(f'{varname1} ensemble mean at {calendar.month_abbr[self.month]}')
        ax.coastlines()
        ax.add_feature(cartopy.feature.LAND, zorder=100)
        pc = ax.pcolormesh(lons, lats, var1, cmap=cmap, transform=ccrs.PlateCarree(),
                          vmin=_vmin[typename][varname1], vmax=_vmax[typename][varname1])
        # ax.set_global()
        fig.colorbar(pc, ax=ax, orientation='horizontal', shrink=0.6, pad=0.02)

        ax = fig.add_subplot(gs[1], projection=ccrs.Robinson())
        ax.set_title(f'{varname2} ensemble mean at {calendar.month_abbr[self.month]}')
        ax.coastlines()
        ax.add_feature(cartopy.feature.LAND, zorder=100)
        pc = ax.pcolormesh(lons, lats, var2, cmap=cmap, transform=ccrs.PlateCarree(),
                          vmin=_vmin[typename][varname2], vmax=_vmax[typename][varname2])
        ax.set_global()
        fig.colorbar(pc, ax=ax, orientation='horizontal', shrink=0.6, pad=0.02)

        ax = fig.add_subplot(gs[2], projection=ccrs.Robinson())
        ax.set_title(f'{varname1}/{varname2} ratio')
        # offset = mcolors.TwoSlopeNorm(vcenter=1., vmin=-1.*vmax, vmax=vmax)
        offset = mcolors.TwoSlopeNorm(vcenter=1.)
        ax.coastlines()
        ax.add_feature(cartopy.feature.LAND, zorder=100)
        pc = ax.pcolormesh(lons, lats, var1/var2, cmap='cmo.balance',
                           transform=ccrs.PlateCarree(), norm=offset)
        ax.set_global()
        fig.colorbar(pc, ax=ax, orientation='horizontal', shrink=0.6, pad=0.02)

        plt.savefig(f'{varname1}_{varname2}_ratio_{self.month}.png', dpi=300)
        plt.close(fig)

    def diff(self, f, varname, typename):
        var1 = self.getEnsMean(varname)
        var2 = f.getEnsMean(varname)
        var1 = np.ma.masked_where(np.isnan(var1) | self.mask, var1)
        var2 = np.ma.masked_where(np.isnan(var2) | self.mask, var2)
        lons = self.lons
        lats = self.lats

        cmap = 'cmo.amp' if typename == 'std' else 'cmo.algae'
        plt.rcParams['font.family'] = 'serif'
        plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']
        plt.rcParams['font.size'] = 12

        fig = plt.figure(1)
        w, h = fig.get_size_inches()
        w = 3*w
        h = h
        fig.set_size_inches(w, h)
        gs = matplotlib.gridspec.GridSpec(1, 3, figure=fig, wspace=0.01, left=0., right=1., bottom=0.0, top=1.)
        fig.clf()

        ax = fig.add_subplot(gs[0], projection=ccrs.Robinson())
        ax.set_title(f'{varname} ensemble mean at {calendar.month_abbr[self.month]} in {self.expname}')
        ax.coastlines()
        ax.add_feature(cartopy.feature.LAND, zorder=100)
        pc = ax.pcolormesh(lons, lats, var1, cmap=cmap, transform=ccrs.PlateCarree(),
                          vmin=_vmin[typename][varname], vmax=_vmax[typename][varname])
        # ax.set_global()
        fig.colorbar(pc, ax=ax, orientation='horizontal', shrink=0.6, pad=0.02)

        ax = fig.add_subplot(gs[1], projection=ccrs.Robinson())
        ax.set_title(f'{varname} ensemble mean at {calendar.month_abbr[self.month]} in {f.expname}')
        ax.coastlines()
        ax.add_feature(cartopy.feature.LAND, zorder=100)
        pc = ax.pcolormesh(lons, lats, var2, cmap=cmap, transform=ccrs.PlateCarree(),
                          vmin=_vmin[typename][varname], vmax=_vmax[typename][varname])
        ax.set_global()
        fig.colorbar(pc, ax=ax, orientation='horizontal', shrink=0.6, pad=0.02)

        ax = fig.add_subplot(gs[2], projection=ccrs.Robinson())
        ax.set_title(f'{self.expname} - {f.expname}')
        # offset = mcolors.TwoSlopeNorm(vcenter=1., vmin=-1.*vmax, vmax=vmax)
        offset = mcolors.TwoSlopeNorm(vcenter=0.,
                                      vmin=_vmin['increment_'+typename][varname],
                                      vmax=_vmax['increment_'+typename][varname])
        ax.coastlines()
        ax.add_feature(cartopy.feature.LAND, zorder=100)
        pc = ax.pcolormesh(lons, lats, var1 - var2, cmap='cmo.balance', transform=ccrs.PlateCarree(), norm=offset)
        ax.set_global()
        fig.colorbar(pc, ax=ax, orientation='horizontal', shrink=0.6, pad=0.02)

        plt.savefig(f'diff_{varname}_ratio_{self.month}.png', dpi=300)
        plt.close(fig)


if __name__ == '__main__':
    suffix1 = '-free'
    year = '2015'
    suffix2 = ''
    expname = 'exp-logmean-chlo'
    def target(ic, month, varname):
        f1 = ForecastFile('2015', ic, str(month).zfill(2), 'exp-free', 'free')
        f2 = ForecastFile('2015', ic, str(month).zfill(2), expname, 'assim Chlorophyll')
        f2.diff(f1, varname, 'state')

    import multiprocessing as mp
    processes = []
    for ic in ['01', '04', '07', '10']:
        for month in range(int(ic), int(ic)+3):
            for varname in [ 'CHD', 'CHN', 'PHD', 'PHN',]:
                p = mp.Process(target=target,
                                    args=(ic, month, varname))
                p.start()
                processes.append(p)

    for p in processes:
        p.join()



