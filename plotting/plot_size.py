"""Comaprison of ratio between phytoplankton constituents.
"""
import itertools
import os

import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature  # type: ignore
import cmocean  # type: ignore # pylint: disable=unused-import
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker
import matplotlib.gridspec as mgs
import matplotlib.pyplot as plt
import numpy as np

import config
import utils


plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman'] + \
    plt.rcParams['font.serif']
plt.rcParams['font.size'] = 18


def get_series(vname: str) -> None:
    """Get monthly series for CO2 and O2 fluxes.
    """
    years: list[str]
    months: list[str]
    years = ['2015', '2016']
    months = [str(i).zfill(2) for i in range(1, 13)]

    n_months: int = len(years) * len(months) - 1
    data: np.ndarray = np.zeros(n_months)

    yearmonths = itertools.product(years, months)
    for i, (year, month) in enumerate(yearmonths):
        if year == '2015' and month == '01':
            continue
        f = np.load(
            os.path.join(
                'data', config.exp,
                'ensmean_'f'{vname}_{year}{month}01.npz'
            )
        )

        data[i-1] = np.ma.masked_array(f['model'], f['model_mask']).sum()

    np.savez(os.path.join('data', config.exp, f'flux_{vname}_series.npz'),
             data=data)


def get_map(vname: str) -> None:
    """Get total concentration of phytoplankton nitrogen
    and chlorophyll over time.
    """
    years: list[str]
    months: list[str]
    years = ['2015', '2016']
    months = [str(i).zfill(2) for i in range(1, 13)]

    n_months: int = len(years) * len(months) - 1
    data: np.ndarray = np.zeros((config.ny, config.nx))

    yearmonths = itertools.product(years, months)
    for i, (year, month) in enumerate(yearmonths):
        if year == '2015' and month == '01':
            continue

        f = np.load(
            os.path.join(
                'data', config.exp,
                'ensmean_'f'{vname}_{year}{month}01.npz'
            )
        )
        data = data + np.ma.masked_array(f['model'], f['model_mask'])

    data = data/n_months

    np.savez(os.path.join('data', config.exp, f'flux_{vname}_map.npz'),
             data=data)


def plot_map(vnames: list[str]) -> None:
    """plot the differences between DA experiment and free run"""
    os.makedirs('figs', exist_ok=True)
    lons, lats = utils.get_coord()
    exps: list[str] = ['free', 'chlo', 'chlo-monthly',
                       'chlo-monthly-update', 'PC',
                       'PC-update', 'chlo-pc',
                       ]

    vmin: dict[str, float] = {
        'PHD': 0.,
        'PHN': 0.,
        'CHD': 0.,
        'CHN': 0.,
        'ZMI': 0.,
        'ZME': 0.
    }
    vmax: dict[str, float] = {
        'PHD': 0.5,
        'PHN': 0.5,
        'CHD': 0.5,
        'CHN': 0.5,
        'ZMI': 0.5,
        'ZME': 0.5
    }
    range_change: dict[str, float] = {
        'PHD': 0.03,
        'PHN': 0.03,
        'CHD': 0.015,
        'CHN': 0.015,
        'ZMI': 0.02,
        'ZME': 0.02
    }

    fig: plt.Figure = plt.figure()
    w, h = fig.get_size_inches()
    fig.set_size_inches(w*2, h*2)
    gs: mgs.GridSpec = mgs.GridSpec(8, 4,
                                    figure=fig,
                                    wspace=0.0,
                                    hspace=0.5,
                                    left=0., right=1.,
                                    bottom=0.01, top=0.94)
    ax: cartopy.mpl.geoaxes.GeoAxes
    for j, vname in enumerate(vnames):
        for i, exp in enumerate(exps):
            if exp == 'free':
                ax = fig.add_subplot(
                    gs[1 + 4*j:3 + 4*j, 0], projection=ccrs.Robinson())
            else:
                m = slice(0 + 4*j, 2 + 4*j) if i < 4 else slice(2 + 4*j, 4 + 4*j)
                n = i if i < 4 else i - 3
                ax = fig.add_subplot(
                    gs[m, n], projection=ccrs.Robinson())
            ax.set_global()
            f = np.load(os.path.join('data', exp, f'flux_{vname}_map.npz'))
            r = f['data']
            norm = mcolors.Normalize(vmin=vmin[vname], vmax=vmax[vname])
            cmap = 'cmo.algae' if vname in [
                'PHD', 'PHN', 'CHD', 'CHN'] else 'cmo.amp'
            if exp != 'free':
                norm = mcolors.TwoSlopeNorm(
                    vcenter=0., vmin=-range_change[vname],
                    vmax=range_change[vname])
                if vname in ['PHD', 'PHN'] and exp == 'chlo-monthly':
                    norm = mcolors.TwoSlopeNorm(
                        vcenter=0., vmin=-0.002, vmax=0.002)
                if vname in ['CHD', 'CHN'] and exp == 'chlo':
                    norm = mcolors.TwoSlopeNorm(
                        vcenter=0., vmin=-0.1, vmax=0.1)
                if vname in ['ZMI', 'ZME'] and exp == 'chlo-monthly':
                    norm = mcolors.TwoSlopeNorm(
                        vcenter=0., vmin=-0.002, vmax=0.002)
                cmap = 'cmo.balance'
                f = np.load(os.path.join('data', 'free',
                                         f'flux_{vname}_map.npz'))
                r = r - f['data']
            r[r < -1e+20] = np.nan
            print(np.nanmin(r), np.nanmax(r))
            pc = ax.pcolormesh(
                lons, lats, np.squeeze(r), transform=ccrs.PlateCarree(),
                cmap=cmap, norm=norm)
            if exp == 'free':
                ax.set_title(f'{vname}\n'+config.exp_labels[exp])
            else:
                ax.set_title(config.exp_labels[exp])
            ax.coastlines(color='k', linewidth=.8)
            ax.add_feature(cfeature.LAND, zorder=3)
            cb = fig.colorbar(pc, ax=ax, orientation='horizontal', pad=0.05,
                              shrink=0.5, norm=norm)
            cb.ax.tick_params(labelsize=12)

    fig.savefig(f'figs/{'_'.join(vnames)}_map.png', dpi=300)


if __name__ == '__main__':
    expnames: list[str] = ['free', 'chlo', 'chlo-monthly',
                           'chlo-monthly-update', 'PC',
                           'PC-update', 'chlo-pc',]
    variables = ['PHD', 'PHN', 'CHD', 'CHN']
    variables = ['ZMI', 'ZME']
    # for expname in expnames:
    #     config.exp = expname
    #     config.output_path = config.output_path_format.format(exp=expname)
    #     for varname in variables:
    #         get_map(varname)

    # plot_map(['CHD', 'CHN'])
    # plot_map(['PHD', 'PHN'])
    plot_map(['ZMI', 'ZME'])
