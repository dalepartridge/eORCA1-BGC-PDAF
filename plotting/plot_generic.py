"""Comaprison of ratio between phytoplankton constituents.
"""
import itertools
import os

import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature  # type: ignore
import cmocean  # type: ignore # pylint: disable=unused-import
import matplotlib.colors as mcolors
import matplotlib.dates as mdates
import matplotlib.gridspec as mgs
import matplotlib.lines as mlines
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


def plot_timeseries(vname: str) -> None:
    """Plot the timeseries for multiple years and months.

    Parameters
    ----------
    vname : str
        The variable name.
    """
    # plotting experiments

    exps: list[str] = ['free', 'chlo', 'chlo-monthly',
                       'chlo-monthly-update', 'PC',
                       'PC-update', 'chlo-pc',
                       ]
    linestyles: list[str] = ['-', ':', '-', ':', '-', ':', ':', ]
    colours: list[str] = ['gray', '#1E88E5', '#FFC107', '#FFC107',
                          '#2CF8BA', '#2CF8BA', '#D81B1B']
    # time array
    t: np.ndarray
    t = np.arange('2015-02', '2017-01', dtype='datetime64[M]')
    print(t)
    # making the plot
    fig: plt.Figure = plt.figure()
    # increase the width of the figure as we will have two subplots
    w: float
    h: float
    w, h = fig.get_size_inches()
    fig.set_size_inches(w, h*1)
    fig.clf()
    gs: mgs.GridSpec = mgs.GridSpec(
        1, 1, figure=fig,
        wspace=0.4, hspace=0., left=0.2, right=0.85,
        bottom=0.33, top=0.93)

    locator: mdates.MonthLocator
    locator = mdates.MonthLocator(bymonth=range(2, 13, 3))

    # loop over the subplots
    ax: plt.Axes = fig.add_subplot(gs[0])
    ax1 = ax.twinx()
    lines: list[mlines.Line2D] = []
    for exp, linestyle, colour in zip(
            exps[1:],
            linestyles[1:],
            colours[1:]):
        f_free = np.load(os.path.join(
            'data', 'free', f'flux_{vname}_series.npz'))
        f = np.load(os.path.join('data', exp, f'flux_{vname}_series.npz'))
        r: np.ndarray = f['data'] - f_free['data']
        line: mlines.Line2D
        line, = ax.plot(t, r, color=colour,
                        linestyle=linestyle,
                        label=config.exp_labels[exp], alpha=1)
        lines.append(line)

    # Set the major locator to be every day
    ax.xaxis.set_major_locator(locator)
    # Set the major formatter to display the date in 'Month-Day' format
    ax.xaxis.set_major_formatter(mdates.AutoDateFormatter(locator))
    ax.tick_params(axis='x', rotation=30)
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    ax.set_xlabel('Time')
    ax.set_ylabel('Differences from freerun')
    ax.set_title(vname)
    line, = ax1.plot(t, f_free['data'], color=colours[0],
                     linestyle=linestyles[0],
                     label=config.exp_labels['free'], alpha=1)
    ax1.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

    fig.legend(loc='outside lower center',
               handles=[mlines.Line2D(
                   [0],
                   [0],
                   linewidth=1,  linestyle=ls, color=colour)
                   for ls, colour in zip(linestyles, colours)],
               labels=[config.exp_labels[exp] for exp in exps],
               ncols=4, fontsize=12, markerscale=0.5)
    fig.savefig(f'flux_{vname}.png', dpi=300)
    plt.close(fig)


def plot_map(vname: str) -> None:
    """plot the differences between DA experiment and free run"""
    os.makedirs('figs', exist_ok=True)
    lons, lats = utils.get_coord()
    exps: list[str] = ['free', 'chlo', 'chlo-monthly',
                       'chlo-monthly-update', 'PC',
                       'PC-update', 'chlo-pc',
                       ]

    fig: plt.Figure = plt.figure()
    w, h = fig.get_size_inches()
    fig.set_size_inches(w*2.6, h*1.2)
    gs: mgs.GridSpec = mgs.GridSpec(4, 4,
                                    figure=fig,
                                    wspace=0.0,
                                    hspace=0.4,
                                    left=0., right=1.,
                                    bottom=0.01, top=0.94)
    ax: cartopy.mpl.geoaxes.GeoAxes
    for i, exp in enumerate(exps):
        if exp == 'free':
            ax = fig.add_subplot(
                gs[1:3, 0], projection=ccrs.Robinson())
        else:
            m = slice(0, 2) if i < 4 else slice(2, 4)
            n = i if i < 4 else i - 3
            ax = fig.add_subplot(
                gs[m, n], projection=ccrs.Robinson())

        ax.set_global()
        f = np.load(os.path.join('data', exp, f'flux_{vname}_map.npz'))
        r = f['data']
        if vname == 'TCO2':
            norm = mcolors.Normalize(vmin=1500., vmax=2200)
        elif vname == 'CO2FLUX':
            norm = mcolors.TwoSlopeNorm(vmin=-100., vcenter=0., vmax=10)
        elif vname == 'OCN_PCO2':
            norm = mcolors.Normalize(vmin=100., vmax=600)
        elif vname == 'O2FLUX':
            norm = mcolors.TwoSlopeNorm(vmin=-50., vcenter=0., vmax=50)
        elif vname == 'OXY':
            norm = mcolors.Normalize(vmin=50., vmax=400)
        elif vname == 'PHD':
            norm = mcolors.Normalize(vmin=0., vmax=0.3)
        elif vname == 'PHN':
            norm = mcolors.Normalize(vmin=0., vmax=0.3)
        elif vname == 'CHD':
            norm = mcolors.Normalize(vmin=0., vmax=3)
        elif vname == 'CHN':
            norm = mcolors.Normalize(vmin=0., vmax=0.5)
        elif vname == 'PRD':
            norm = mcolors.Normalize(vmin=0., vmax=2)
        elif vname == 'PRN':
            norm = mcolors.Normalize(vmin=0., vmax=5)
        else:
            norm = mcolors.Normalize(vmin=0., vmax=1.)
        cmap = 'cmo.amp' if isinstance(
            norm, mcolors.Normalize) else 'cmo.balance'
        if exp != 'free':
            if vname == 'TCO2':
                norm = mcolors.TwoSlopeNorm(
                    vcenter=0., vmin=-1, vmax=1.)
                if exp == 'chlo-monthly':
                    norm = mcolors.TwoSlopeNorm(
                        vcenter=0., vmin=-0.1, vmax=0.1)
            elif vname == 'CO2FLUX':
                norm = mcolors.TwoSlopeNorm(
                    vcenter=0., vmin=-0.6, vmax=0.6)
            elif vname == 'OCN_PCO2':
                norm = mcolors.TwoSlopeNorm(
                    vcenter=0., vmin=-3, vmax=3)
            elif vname == 'O2FLUX':
                norm = mcolors.TwoSlopeNorm(
                    vcenter=0., vmin=-1, vmax=1)
            elif vname == 'OXY':
                norm = mcolors.TwoSlopeNorm(
                    vcenter=0., vmin=-0.3, vmax=0.3)
            elif vname == 'PHD':
                norm = mcolors.TwoSlopeNorm(vcenter=0., vmin=-0.01, vmax=0.01)
            elif vname == 'PHN':
                norm = mcolors.TwoSlopeNorm(vcenter=0., vmin=-0.01, vmax=0.01)
            elif vname == 'CHD':
                norm = mcolors.TwoSlopeNorm(vcenter=0., vmin=-0.01, vmax=0.01)
            elif vname == 'CHN':
                norm = mcolors.TwoSlopeNorm(vcenter=0., vmin=-0.01, vmax=0.01)
            elif vname == 'PRD':
                norm = mcolors.TwoSlopeNorm(vcenter=0., vmin=-0.2, vmax=0.2)
            elif vname == 'PRN':
                norm = mcolors.TwoSlopeNorm(vcenter=0., vmin=-0.5, vmax=0.5)
            elif vname == 'DET':
                norm = mcolors.TwoSlopeNorm(vcenter=0., vmin=-0.01, vmax=0.01)
            elif vname == 'MPN':
                norm = mcolors.TwoSlopeNorm(vcenter=0., vmin=-0.02, vmax=0.02)
            elif vname == 'MPD':
                norm = mcolors.TwoSlopeNorm(vcenter=0., vmin=-0.1, vmax=0.1)
            elif vname == 'GMIPn':
                norm = mcolors.TwoSlopeNorm(vcenter=0., vmin=-0.7, vmax=0.7)
            elif vname == 'GMEPN':
                norm = mcolors.TwoSlopeNorm(vcenter=0., vmin=-0.1, vmax=0.1)
            elif vname == 'GMEPD':
                norm = mcolors.TwoSlopeNorm(vcenter=0., vmin=-0.1, vmax=0.1)
            elif vname == 'PD_JLIM':
                norm = mcolors.TwoSlopeNorm(vcenter=0., vmin=-0.1, vmax=0.1)
            elif vname == 'PD_NLIM':
                norm = mcolors.TwoSlopeNorm(vcenter=0., vmin=-0.1, vmax=0.1)
            elif vname == 'PD_FELIM':
                norm = mcolors.TwoSlopeNorm(vcenter=0., vmin=-0.1, vmax=0.1)
            elif vname == 'PD_JLIM':
                norm = mcolors.TwoSlopeNorm(vcenter=0., vmin=-0.1, vmax=0.1)
            elif vname == 'PD_NLIM':
                norm = mcolors.TwoSlopeNorm(vcenter=0., vmin=-0.1, vmax=0.1)
            elif vname == 'PD_FELIM':
                norm = mcolors.TwoSlopeNorm(vcenter=0., vmin=-0.1, vmax=0.1)
            elif vname == 'PD_SILIM':
                norm = mcolors.TwoSlopeNorm(vcenter=0., vmin=-20, vmax=20)
            elif vname == 'PDSILIM2':
                norm = mcolors.TwoSlopeNorm(vcenter=0., vmin=-0.1, vmax=0.1)
            elif vname == 'PN_JLIM':
                norm = mcolors.TwoSlopeNorm(vcenter=0., vmin=-0.1, vmax=0.1)
            elif vname == 'PN_NLIM':
                norm = mcolors.TwoSlopeNorm(vcenter=0., vmin=-0.9, vmax=0.9)
            elif vname == 'PN_FELIM':
                norm = mcolors.TwoSlopeNorm(vcenter=0., vmin=-0.9, vmax=0.9)
            elif vname == 'PD_LLOSS':
                norm = mcolors.TwoSlopeNorm(vcenter=0., vmin=-0.01, vmax=0.01)
            elif vname == 'PN_LLOSS':
                norm = mcolors.TwoSlopeNorm(vcenter=0., vmin=-0.01, vmax=0.01)
            cmap = 'cmo.balance'
            f = np.load(os.path.join('data', 'free',
                                     f'flux_{vname}_map.npz'))
            r = r - f['data']
        r[r < -1e+20] = np.nan
        print(np.nanmin(r), np.nanmax(r))
        pc = ax.pcolormesh(
            lons, lats, np.squeeze(r), transform=ccrs.PlateCarree(),
            cmap=cmap, norm=norm)
        ax.set_title(config.exp_labels[exp])
        ax.coastlines(color='k', linewidth=.8)
        ax.add_feature(cfeature.LAND, zorder=3)
        fig.colorbar(pc, ax=ax, orientation='horizontal', pad=0.05,
                     shrink=0.5, norm=norm)

    fig.savefig(f'figs/flux_{vname}_map.png', dpi=300)


if __name__ == '__main__':
    expnames: list[str] = ['free', 'chlo', 'chlo-monthly',
                           'chlo-monthly-update', 'PC',
                           'PC-update', 'chlo-pc',]
    variables = ['TCO2', 'CO2FLUX', 'OCN_PCO2', 'O2FLUX', 'OXY', ]
    # variables = ['PHD', 'PHN', 'CHD', 'CHN']
    variables = ['PRD', 'PRN']
    variables = ['DET', 'MPN', 'MPD', 'GMIPn', 'GMEPN', 'GMEPD', 'PD_JLIM',
                 'PD_NLIM', 'PD_FELIM', 'PD_SILIM', 'PDSILIM2', 'PN_JLIM',
                 'PN_NLIM', 'PN_FELIM']
    variables = ['PD_LLOSS', 'PN_LLOSS']
    for expname in expnames:
        config.exp = expname
        config.output_path = config.output_path_format.format(exp=expname)
        for varname in variables:
            get_series(varname)
            get_map(varname)

    for varname in variables:
        # plot_timeseries(vname)
        print(varname)
        plot_map(varname)
