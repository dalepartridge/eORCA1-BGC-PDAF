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
import diag
import file_info
import utils

plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman'] + \
    plt.rcParams['font.serif']
plt.rcParams['font.size'] = 18


def get_n_chl_series(is_obs: bool) -> None:
    """Get monthly series for ratio of phytoplankton nitrogen
    and chlorophyll.
    """
    years: list[str]
    months: list[str]
    years = ['2015', '2016']
    months = [str(i).zfill(2) for i in range(1, 13)]

    n_months: int = len(years) * len(months) - 1
    ratio: np.ndarray = np.zeros(n_months)

    yearmonths = itertools.product(years, months)
    for i, (year, month) in enumerate(yearmonths):
        if year == '2015' and month == '01':
            continue
        # Read the physical observation
        f_info: file_info.FileInfo
        f_info = file_info.FileInfo(year, month, '01', '')
        obs_n, _ = diag.read_physical_observation(
            f_info.get_obsfilename('pc'), 'pc')
        obs_chl, _ = diag.read_physical_observation(
            f_info.get_obsfilename('chlo-monthly'), 'chlo-monthly')
        mask = np.logical_or(obs_n.mask, obs_chl.mask)
        if is_obs:
            ratio[i-1] = (12*np.ma.masked_array(obs_n.data, mask) /
                          np.ma.masked_array(obs_chl.data, mask)).mean()
        else:
            f_chl = np.load(
                os.path.join(
                    'data', config.exp,
                    'ensmean_chlo-monthly'f'_{year}{month}01.npz'
                )
            )
            f_n = np.load(
                os.path.join(
                    'data', config.exp,
                    'ensmean_nitrogen-monthly'f'_{year}{month}01.npz'
                )
            )
            ratio[i-1] = (12*np.ma.masked_array(f_n['model'], mask) /
                          np.ma.masked_array(f_chl['model'], mask)).mean()

    if is_obs:
        os.makedirs(os.path.join('data', 'obs'), exist_ok=True)
        np.savez(os.path.join('data', 'obs', 'mor_n_chl_series.npz'),
                 ratio=ratio)
    else:
        np.savez(os.path.join('data', config.exp, 'mor_n_chl_series.npz'),
                 ratio=ratio)


def get_chl_n_series(is_obs: bool) -> None:
    """Get monthly series for ratio of phytoplankton nitrogen
    and chlorophyll.
    """
    years: list[str]
    months: list[str]
    years = ['2015', '2016']
    months = [str(i).zfill(2) for i in range(1, 13)]

    n_months: int = len(years) * len(months) - 1
    ratio: np.ndarray = np.zeros(n_months)

    yearmonths = itertools.product(years, months)
    for i, (year, month) in enumerate(yearmonths):
        if year == '2015' and month == '01':
            continue
        # Read the physical observation
        f_info: file_info.FileInfo
        f_info = file_info.FileInfo(year, month, '01', '')
        obs_n, _ = diag.read_physical_observation(
            f_info.get_obsfilename('pc'), 'pc')
        obs_chl, _ = diag.read_physical_observation(
            f_info.get_obsfilename('chlo-monthly'), 'chlo-monthly')
        mask = np.logical_or(obs_n.mask, obs_chl.mask)
        if is_obs:
            ratio[i-1] = (np.ma.masked_array(obs_chl.data, mask) /
                          12*np.ma.masked_array(obs_n.data, mask)).mean()
        else:
            f_chl = np.load(
                os.path.join(
                    'data', config.exp,
                    'ensmean_chlo-monthly'f'_{year}{month}01.npz'
                )
            )
            f_n = np.load(
                os.path.join(
                    'data', config.exp,
                    'ensmean_nitrogen-monthly'f'_{year}{month}01.npz'
                )
            )
            ratio[i-1] = (np.ma.masked_array(f_chl['model'], mask) /
                          12*np.ma.masked_array(f_n['model'], mask)).mean()

    if is_obs:
        os.makedirs(os.path.join('data', 'obs'), exist_ok=True)
        np.savez(os.path.join('data', 'obs', 'mor_chl_n_series.npz'),
                 ratio=ratio)
    else:
        np.savez(os.path.join('data', config.exp, 'mor_chl_n_series.npz'),
                 ratio=ratio)


def get_n_chl_map(is_obs: bool) -> None:
    """Get map for ratio of phytoplankton nitrogen
    and chlorophyll.
    """
    years: list[str]
    months: list[str]
    years = ['2015', '2016']
    months = [str(i).zfill(2) for i in range(1, 13)]

    ratio: np.ndarray = np.zeros((config.ny, config.nx))
    n_values: np.ndarray = np.zeros((config.ny, config.nx))

    yearmonths = itertools.product(years, months)
    for i, (year, month) in enumerate(yearmonths):
        if year == '2015' and month == '01':
            continue
        # Read the physical observation
        f_info: file_info.FileInfo
        f_info = file_info.FileInfo(year, month, '01', '')
        obs_n, _ = diag.read_physical_observation(
            f_info.get_obsfilename('pc'), 'pc')
        obs_chl, _ = diag.read_physical_observation(
            f_info.get_obsfilename('chlo-monthly'), 'chlo-monthly')
        mask = np.logical_or(obs_n.mask, obs_chl.mask)
        if is_obs:
            r = 12*np.ma.masked_array(obs_n.data, mask) / \
                np.ma.masked_array(obs_chl.data, mask)
            mask[np.isnan(r) | np.isinf(r)] = True
            r[np.isnan(r) | np.isinf(r)] = 0
            ratio = ratio + r.filled(0)
            n_values = n_values + np.logical_not(mask)
        else:
            f_chl = np.load(
                os.path.join(
                    'data', config.exp,
                    'ensmean_chlo-monthly'f'_{year}{month}01.npz'
                )
            )
            f_n = np.load(
                os.path.join(
                    'data', config.exp,
                    'ensmean_nitrogen-monthly'f'_{year}{month}01.npz'
                )
            )
            r = 12*np.ma.masked_array(f_n['model'], f_n['model_mask']) / \
                np.ma.masked_array(f_chl['model'], f_chl['model_mask'])
            n_values = n_values + np.logical_not(
                f_n['model_mask'] & f_chl['model_mask'] & np.isnan(r) & np.isinf(
                    r)).astype(int)
            r[np.isnan(r) | np.isinf(r)] = 0
            ratio = ratio + r.filled(0)

    if is_obs:
        os.makedirs(os.path.join('data', 'obs'), exist_ok=True)
        np.savez(os.path.join('data', 'obs', 'mor_n_chl_map.npz'),
                 ratio=ratio/n_values)
    else:
        np.savez(os.path.join('data', config.exp, 'mor_n_chl_map.npz'),
                 ratio=ratio/n_values)


def get_phytoplankton_series() -> None:
    """Get monthly series for ratio of diatom and non-diatom
    phytoplankton nitrogen.
    """
    years: list[str]
    months: list[str]
    years = ['2015', '2016']
    months = [str(i).zfill(2) for i in range(1, 13)]

    n_months: int = len(years) * len(months) - 1
    data: dict[str, np.ndarray] = {}
    data['r_chl'] = np.zeros(n_months)
    data['r_n'] = np.zeros(n_months)

    yearmonths = itertools.product(years, months)
    for i, (year, month) in enumerate(yearmonths):
        if year == '2015' and month == '01':
            continue

        f_d = np.load(
            os.path.join(
                'data', config.exp,
                f'ensmean_CHD_{year}{month}01.npz'
            )
        )
        f_nd = np.load(
            os.path.join(
                'data', config.exp,
                f'ensmean_CHN_{year}{month}01.npz'
            )
        )
        data['r_chl'][i-1] = (np.ma.masked_array(f_d['model'],
                                                 f_d['model_mask']) /
                              np.ma.masked_array(f_nd['model'],
                                                 f_nd['model_mask'])).mean()

        f_d = np.load(
            os.path.join(
                'data', config.exp,
                f'ensmean_PHD_{year}{month}01.npz'
            )
        )
        f_nd = np.load(
            os.path.join(
                'data', config.exp,
                f'ensmean_PHN_{year}{month}01.npz'
            )
        )
        data['r_n'][i-1] = (np.ma.masked_array(f_d['model'],
                                               f_d['model_mask']) /
                            np.ma.masked_array(f_nd['model'],
                                               f_nd['model_mask'])).mean()

    np.savez(
        os.path.join(
            'data', config.exp, 'mor_time_series_phytoplankton.npz'),
        **data)  # type: ignore


def get_phytoplankton_map() -> None:
    """Get monthly series for ratio of diatom and non-diatom
    phytoplankton nitrogen.
    """
    years: list[str]
    months: list[str]
    years = ['2015', '2016']
    months = [str(i).zfill(2) for i in range(1, 13)]

    data: dict[str, np.ndarray] = {}
    data['r_chl'] = np.zeros((config.ny, config.nx))
    data['r_n'] = np.zeros((config.ny, config.nx))
    n_values_chl = np.zeros((config.ny, config.nx))
    n_values_n = np.zeros((config.ny, config.nx))

    yearmonths = itertools.product(years, months)
    for i, (year, month) in enumerate(yearmonths):
        if year == '2015' and month == '01':
            continue

        f_d = np.load(
            os.path.join(
                'data', config.exp,
                f'ensmean_CHD_{year}{month}01.npz'
            )
        )
        f_nd = np.load(
            os.path.join(
                'data', config.exp,
                f'ensmean_CHN_{year}{month}01.npz'
            )
        )
        r = np.squeeze(np.ma.masked_array(f_d['model'],
                                          f_d['model_mask']) /
                       np.ma.masked_array(f_nd['model'],
                                          f_nd['model_mask']))
        n_values_chl = n_values_chl + np.logical_not(np.squeeze(f_d['model_mask']) & np.squeeze(
            f_nd['model_mask']) & np.isnan(r) & np.isinf(r)).astype(int)
        r[np.isnan(r) | np.isinf(r)] = 0
        data['r_chl'] += r

        f_d = np.load(
            os.path.join(
                'data', config.exp,
                f'ensmean_PHD_{year}{month}01.npz'
            )
        )
        f_nd = np.load(
            os.path.join(
                'data', config.exp,
                f'ensmean_PHN_{year}{month}01.npz'
            )
        )
        r = np.squeeze(np.ma.masked_array(f_d['model'],
                                          f_d['model_mask']) /
                       np.ma.masked_array(f_nd['model'],
                                          f_nd['model_mask']))
        n_values_n = n_values_n + np.logical_not(np.squeeze(f_d['model_mask']) & np.squeeze(
            f_nd['model_mask']) & np.isnan(r) & np.isinf(r)).astype(int)
        r[np.isnan(r) | np.isinf(r)] = 0
        data['r_n'] += r

    data['r_chl'] /= n_values_chl
    data['r_n'] /= n_values_n

    np.savez(
        os.path.join(
            'data', config.exp, 'mor_map_phytoplankton.npz'),
        **data)  # type: ignore


def plot_timeseries() -> None:
    """Plot the timeseries of ratios for multiple years and months.

    Parameters
    ----------
    varname : str
        The variable name.
    """
    # plotting experiments

    exps: list[str] = ['obs', 'free', 'chlo', 'chlo-monthly',
                       'chlo-monthly-update', 'PC',
                       'PC-update', 'chlo-pc',
                       ]
    linestyles: list[str] = ['-', '-', ':', '-', ':', '-', ':', ':', ]
    colours: list[str] = ['k', 'gray', '#1E88E5', '#FFC107', '#FFC107',
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
    fig.set_size_inches(3*w, h*1)
    fig.clf()
    gs: mgs.GridSpec = mgs.GridSpec(
        1, 3, figure=fig,
        wspace=0.1, hspace=0., left=0.04, right=1.,
        bottom=0.36, top=0.92)

    locator: mdates.MonthLocator
    locator = mdates.MonthLocator(bymonth=range(2, 13, 3))

    # loop over the subplots
    ax: plt.Axes = fig.add_subplot(gs[0])
    lines: list[mlines.Line2D] = []
    for exp, linestyle, colour in zip(
            exps,
            linestyles,
            colours):
        # f = np.load(os.path.join('data', exp, 'mor_n_chl_series.npz'))
        f = np.load(os.path.join('data', exp, 'mor_chl_n_series.npz'))
        r: np.ndarray = f['ratio']
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
    ax.set_xlabel('Time')
    ax.set_title(r'$Chl/N$ ratio')

    ax = fig.add_subplot(gs[1])
    for exp, linestyle, colour in zip(
            exps[1:],
            linestyles[1:],
            colours[1:]):
        f = np.load(os.path.join(
            'data', exp, 'mor_time_series_phytoplankton.npz'))
        r = f['r_chl']
        line, = ax.plot(t, r, color=colour,
                        linestyle=linestyle,
                        label=config.exp_labels[exp], alpha=1)

    # Set the major locator to be every day
    ax.xaxis.set_major_locator(locator)
    # Set the major formatter to display the date in 'Month-Day' format
    ax.xaxis.set_major_formatter(mdates.AutoDateFormatter(locator))
    ax.tick_params(axis='x', rotation=30)
    ax.set_xlabel('Time')
    ax.set_title(r'$Chl_{diatom}/Chl_{non-diatom}$ ratio')

    ax = fig.add_subplot(gs[2])
    for exp, linestyle, colour in zip(
            exps[1:],
            linestyles[1:],
            colours[1:]):
        f = np.load(os.path.join(
            'data', exp, 'mor_time_series_phytoplankton.npz'))
        r = f['r_n']
        line, = ax.plot(t, r, color=colour,
                        linestyle=linestyle,
                        label=config.exp_labels[exp], alpha=1)

    # Set the major locator to be every day
    ax.xaxis.set_major_locator(locator)
    # Set the major formatter to display the date in 'Month-Day' format
    ax.xaxis.set_major_formatter(mdates.AutoDateFormatter(locator))
    ax.tick_params(axis='x', rotation=30)
    ax.set_xlabel('Time')
    ax.set_title(r'$N_{diatom}/N_{non-diatom}$ ratio')

    fig.legend(loc='outside lower center',
               handles=[mlines.Line2D(
                   [0],
                   [0],
                   linewidth=1,  linestyle=ls, color=colour)
                   for ls, colour in zip(linestyles, colours)],
               labels=[config.exp_labels[exp] for exp in exps],
               ncols=4, fontsize=12, markerscale=0.5)
    fig.savefig('mean_of_ratio_chl_n_ts.pdf', dpi=300)
    plt.close(fig)


def plot_axes_n_chl(fig, ax, exp: str) -> None:
    """plot the axes"""
    lons, lats = utils.get_coord()
    ax.set_global()
    f = np.load(os.path.join('data', exp, 'mor_n_chl_map.npz'))
    r = np.squeeze(f['ratio'])
    norm = mcolors.Normalize(vmin=3., vmax=50)
    cmap = 'cmo.amp'
    if exp == 'free':
        norm = mcolors.SymLogNorm(  # pylint: disable=unexpected-keyword-arg
            linthresh=12, vmin=-90, vmax=90)
        cmap = 'cmo.balance'
        f = np.load(os.path.join('data', 'obs', 'mor_n_chl_map.npz'))
        r = r - f['ratio']
    elif exp != 'obs':
        norm = mcolors.SymLogNorm(  # pylint: disable=unexpected-keyword-arg
            linthresh=0.5, vmin=-80., vmax=80)
        cmap = 'cmo.balance'
        f = np.load(os.path.join('data', 'free', 'mor_n_chl_map.npz'))
        r = r - f['ratio']
    print(np.nanmin(r), np.nanmax(r))
    pc = ax.pcolormesh(
        lons, lats, np.squeeze(r), transform=ccrs.PlateCarree(),
        cmap=cmap, norm=norm)
    ax.set_title(config.exp_labels[exp])
    ax.coastlines(color='k', linewidth=.8)
    ax.add_feature(cfeature.LAND, zorder=3)
    fig.colorbar(pc, ax=ax, orientation='horizontal', pad=0.05,
                 shrink=0.5, norm=norm)


def plot_axes_phytoplanktons(
        fig, ax, exp: str, vname: str) -> None:
    """plot the axes"""
    lons, lats = utils.get_coord()
    ax.set_global()
    f = np.load(os.path.join('data', exp, 'mor_map_phytoplankton.npz'))
    r = np.squeeze(f[f'r_{vname}'])
    norm = mcolors.Normalize(vmin=0., vmax=1.5)
    cmap = 'cmo.amp'
    if exp != 'free':
        if vname == 'n':
            norm = mcolors.TwoSlopeNorm(vmin=-0.2, vcenter=0., vmax=0.2)
        else:
            norm = mcolors.TwoSlopeNorm(vmin=-0.1, vcenter=0., vmax=0.1)
        cmap = 'cmo.balance'
        f = np.load(os.path.join('data', 'free', 'mor_map_phytoplankton.npz'))
        r = r - f[f'r_{vname}']
    print(np.nanmin(r), np.nanmax(r))
    pc = ax.pcolormesh(
        lons, lats, np.squeeze(r), transform=ccrs.PlateCarree(),
        cmap=cmap, norm=norm)
    ax.set_title(config.exp_labels[exp])
    ax.coastlines(color='k', linewidth=.8)
    ax.add_feature(cfeature.LAND, zorder=3)
    fig.colorbar(pc, ax=ax, orientation='horizontal', pad=0.05,
                 shrink=0.5, norm=norm)


def plot_map(vname: str) -> None:
    """plot the differences between DA experiment and free run"""
    os.makedirs('figs', exist_ok=True)
    exps: list[str] = ['obs', 'free', 'chlo', 'chlo-monthly',
                       'chlo-monthly-update', 'PC',
                       'PC-update', 'chlo-pc',
                       ]

    fig: plt.Figure = plt.figure()
    w, h = fig.get_size_inches()
    fig.set_size_inches(w*2.6, h*1.2)
    n_subplots = 2 if vname == 'n_chl' else 4
    hspace = 0.15 if vname == 'n_chl' else 0.4
    gs: mgs.GridSpec = mgs.GridSpec(n_subplots, 4,
                                    figure=fig,
                                    wspace=0.0,
                                    hspace=hspace,
                                    left=0., right=1.,
                                    bottom=0.01, top=0.94)
    for i, exp in enumerate(exps):
        if vname == 'n_chl':
            ax = fig.add_subplot(
                gs[i], projection=ccrs.Robinson())
            plot_axes_n_chl(fig, ax, exp)
        else:
            if exp == 'obs':
                continue
            if exp == 'free':
                ax = fig.add_subplot(
                    gs[1:3, 0], projection=ccrs.Robinson())
            else:
                m = slice(0, 2) if i <= 4 else slice(2, 4)
                n = i - 1 if i <= 4 else i - 4
                ax = fig.add_subplot(
                    gs[m, n], projection=ccrs.Robinson())
            plot_axes_phytoplanktons(fig, ax, exp, vname)

    fig.savefig(f'figs/mean_of_ratio_map_{vname}.png', dpi=300)


def plot_monthly_map() -> None:
    """Plot the monthly map for the ratio of N/Chl."""
    os.makedirs('figs/monthly', exist_ok=True)
    lons, lats = utils.get_coord()

    exps: list[str] = ['obs', 'free', 'chlo', 'chlo-monthly',
                       'chlo-monthly-update', 'PC',
                       'PC-update', 'chlo-pc',
                       ]
    years: list[str]
    months: list[str]
    years = ['2015', '2016']
    months = [str(i).zfill(2) for i in range(1, 13)]

    yearmonths = itertools.product(years, months)
    for i, (year, month) in enumerate(yearmonths):
        if year == '2015' and month == '01':
            continue

        fig: plt.Figure = plt.figure()
        w, h = fig.get_size_inches()
        fig.set_size_inches(w*3, h*1.2)
        gs: mgs.GridSpec = mgs.GridSpec(2, 4,
                                        figure=fig,
                                        wspace=0.01,
                                        hspace=0.01,
                                        left=0., right=1.,
                                        bottom=0.08, top=0.98)
        for i, exp in enumerate(exps):
            config.exp = exp
            config.output_path = config.output_path_format.format(exp=exp)
            ax: cartopy.mpl.geoaxes.GeoAxes = fig.add_subplot(
                gs[i], projection=ccrs.Robinson())
            ax.set_global()
            # Read the physical observation
            f_info: file_info.FileInfo
            f_info = file_info.FileInfo(year, month, '01', '')
            obs_n, _ = diag.read_physical_observation(
                f_info.get_obsfilename('pc'), 'pc')
            obs_chl, _ = diag.read_physical_observation(
                f_info.get_obsfilename('chlo-monthly'), 'chlo-monthly')
            mask = np.logical_or(obs_n.mask, obs_chl.mask)
            if exp == 'obs':
                r = 12*np.ma.masked_array(obs_n.data, mask) / \
                    np.ma.masked_array(obs_chl.data, mask)
            else:
                f_chl = np.load(
                    os.path.join(
                        'data', config.exp,
                        'ensmean_chlo-monthly'f'_{year}{month}01.npz'
                    )
                )
                f_n = np.load(
                    os.path.join(
                        'data', config.exp,
                        'ensmean_nitrogen-monthly'f'_{year}{month}01.npz'
                    )
                )
                r = 12*np.ma.masked_array(f_n['model'], f_n['model_mask']) / \
                    np.ma.masked_array(f_chl['model'], f_chl['model_mask'])
            print(r.min(), r.max())
            ax.pcolormesh(
                lons, lats, np.squeeze(r), transform=ccrs.PlateCarree(),
                cmap='cmo.amp', norm=mcolors.Normalize(vmin=9., vmax=60))
            ax.set_title('N/Chl ratio in '+config.exp_labels[exp])
            ax.coastlines(color='k', linewidth=.8)
            ax.add_feature(cfeature.LAND, zorder=3)
        fig.savefig(
            f'figs/monthly/mean_of_ratio_map_{year}{month}.png', dpi=300)


if __name__ == '__main__':
    expnames: list[str] = ['free', 'chlo', 'chlo-monthly',
                           'chlo-monthly-update', 'PC',
                           'PC-update', 'chlo-pc',]
    for expname in expnames:
        config.exp = expname
        config.output_path = config.output_path_format.format(exp=expname)
        # get_n_chl_series(False)
        get_chl_n_series(False)
        # get_phytoplankton_series()
    #     get_n_chl_map(False)
    #     get_phytoplankton_map()

    config.exp = ''
    config.output_path = config.output_path_format.format(exp='')
    get_chl_n_series(True)
    # # get_n_chl_series(True)
    # get_n_chl_map(True)

    # plot_timeseries()
    # plot_map('n_chl')
    # plot_map('n')
    plot_map('chl')
    # plot_monthly_map()
