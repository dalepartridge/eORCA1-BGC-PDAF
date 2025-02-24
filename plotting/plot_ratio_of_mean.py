"""Comaprison of ratio between phytoplankton constituents.
"""
import itertools
import os

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
    chl: np.ndarray = np.zeros(n_months)
    n: np.ndarray = np.zeros(n_months)

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
            n[i-1] = 12*np.ma.masked_array(obs_n.data, mask).sum()
            chl[i-1] = np.ma.masked_array(obs_chl.data, mask).sum()
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
            n[i-1] = 12*np.ma.masked_array(f_n['model'], mask).sum()
            chl[i-1] = np.ma.masked_array(f_chl['model'], mask).sum()

    if is_obs:
        os.makedirs(os.path.join('data', 'obs'), exist_ok=True)
        np.savez(os.path.join('data', 'obs', 'n_chl_series.npz'),
                 n=n, chl=chl)
    else:
        np.savez(os.path.join('data', config.exp, 'n_chl_series.npz'),
                 n=n, chl=chl)


def get_n_chl_map(is_obs: bool) -> None:
    """Get total concentration of phytoplankton nitrogen
    and chlorophyll over time.
    """
    years: list[str]
    months: list[str]
    years = ['2015', '2016']
    months = [str(i).zfill(2) for i in range(1, 13)]

    chl: np.ndarray = np.zeros((config.ny, config.nx))
    n: np.ndarray = np.zeros((config.ny, config.nx))
    n_values: np.ndarray = np.zeros((config.ny, config.nx), dtype=int)

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
            n = n + 12*np.ma.masked_array(obs_n.data, mask).filled(0.)
            chl = chl + np.ma.masked_array(obs_chl.data, mask).filled(0.)
            n_values = n_values + np.logical_not(mask).astype(int)
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
            n = n + 12*np.ma.masked_array(f_n['model'], f_n['model_mask'])
            chl = chl + np.ma.masked_array(f_chl['model'], f_chl['model_mask'])
            n_values = n_values + np.logical_not(
                f_n['model_mask'] & f_chl['model_mask']).astype(int)

    if is_obs:
        os.makedirs(os.path.join('data', 'obs'), exist_ok=True)
        np.savez(os.path.join('data', 'obs', 'rom_n_chl_map.npz'),
                 n=n/n_values, chl=chl/n_values)
    else:
        np.savez(os.path.join('data', config.exp, 'rom_n_chl_map.npz'),
                 n=n/n_values, chl=chl/n_values)


def get_phytoplankton_series() -> None:
    """Get monthly series for ratio of diatom and non-diatom
    phytoplankton constituents.
    """
    years: list[str]
    months: list[str]
    years = ['2015', '2016']
    months = [str(i).zfill(2) for i in range(1, 13)]

    n_months: int = len(years) * len(months) - 1
    data: dict[str, np.ndarray] = {}
    data['CHD'] = np.zeros(n_months)
    data['CHN'] = np.zeros(n_months)
    data['PHD'] = np.zeros(n_months)
    data['PHN'] = np.zeros(n_months)

    yearmonths = itertools.product(years, months)
    for i, (year, month) in enumerate(yearmonths):
        if year == '2015' and month == '01':
            continue

        for varname in data:  # pylint: disable=consider-using-dict-items
            f = np.load(
                os.path.join(
                    'data', config.exp,
                    f'ensmean_{varname}_{year}{month}01.npz'
                )
            )
            data[varname][i-1] = np.ma.masked_array(f['model'],
                                                    f['model_mask']).sum()

    np.savez(os.path.join('data', config.exp, 'time_series_phytoplankton.npz'),
             **data)  # type: ignore


def get_phytoplankton_map() -> None:
    """Get monthly map for diatom and non-diatom
    phytoplankton constituents.
    """
    years: list[str]
    months: list[str]
    years = ['2015', '2016']
    months = [str(i).zfill(2) for i in range(1, 13)]

    n_months: int = len(years) * len(months) - 1
    data: dict[str, np.ndarray] = {}
    data['CHD'] = np.zeros((config.ny, config.nx))
    data['CHN'] = np.zeros((config.ny, config.nx))
    data['PHD'] = np.zeros((config.ny, config.nx))
    data['PHN'] = np.zeros((config.ny, config.nx))

    yearmonths = itertools.product(years, months)
    for i, (year, month) in enumerate(yearmonths):
        if year == '2015' and month == '01':
            continue

        for varname in data:  # pylint: disable=consider-using-dict-items
            f = np.load(
                os.path.join(
                    'data', config.exp,
                    f'ensmean_{varname}_{year}{month}01.npz'
                )
            )
            data[varname] += np.squeeze(np.ma.masked_array(f['model'],
                                                           f['model_mask']))

    for varname in data:  # pylint: disable=consider-using-dict-items
        data[varname] = data[varname]/n_months

    np.savez(os.path.join('data', config.exp, 'rom_map_phytoplankton.npz'),
             **data)  # type: ignore


def get_zooplankton_series() -> None:
    """Get monthly series for ratio of meso- and micro-
    zooplanktons.
    """
    years: list[str]
    months: list[str]
    years = ['2015', '2016']
    months = [str(i).zfill(2) for i in range(1, 13)]

    n_months: int = len(years) * len(months) - 1
    data: dict[str, np.ndarray] = {}
    data['ZMI'] = np.zeros(n_months)
    data['ZME'] = np.zeros(n_months)

    yearmonths = itertools.product(years, months)
    for i, (year, month) in enumerate(yearmonths):
        if year == '2015' and month == '01':
            continue

        for varname in data:  # pylint: disable=consider-using-dict-items
            f = np.load(
                os.path.join(
                    'data', config.exp,
                    f'ensmean_{varname}_{year}{month}01.npz'
                )
            )
            data[varname][i-1] = np.ma.masked_array(f['model'],
                                                    f['model_mask']).sum()

    np.savez(os.path.join('data', config.exp, 'time_series_zooplankton.npz'),
             **data)  # type: ignore


def get_zooplankton_map() -> None:
    """Get monthly map for meso- and micro-zooplankton.
    """
    years: list[str]
    months: list[str]
    years = ['2015', '2016']
    months = [str(i).zfill(2) for i in range(1, 13)]

    n_months: int = len(years) * len(months) - 1
    data: dict[str, np.ndarray] = {}
    data['ZMI'] = np.zeros((config.ny, config.nx))
    data['ZME'] = np.zeros((config.ny, config.nx))

    yearmonths = itertools.product(years, months)
    for i, (year, month) in enumerate(yearmonths):
        if year == '2015' and month == '01':
            continue

        for varname in data:  # pylint: disable=consider-using-dict-items
            f = np.load(
                os.path.join(
                    'data', config.exp,
                    f'ensmean_{varname}_{year}{month}01.npz'
                )
            )
            data[varname] += np.squeeze(np.ma.masked_array(f['model'],
                                                           f['model_mask']))

    for varname in data:  # pylint: disable=consider-using-dict-items
        data[varname] = data[varname]/n_months

    np.savez(os.path.join('data', config.exp, 'rom_map_zooplankton.npz'),
             **data)  # type: ignore


def plot_phyto_timeseries() -> None:
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
        f = np.load(os.path.join('data', exp, 'n_chl_series.npz'))
        r: np.ndarray = f['n']/f['chl']
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
    ax.set_title(r'$N/Chl$ ratio')

    ax = fig.add_subplot(gs[1])
    for exp, linestyle, colour in zip(
            exps[1:],
            linestyles[1:],
            colours[1:]):
        f = np.load(os.path.join('data', exp, 'time_series_phytoplankton.npz'))
        r = f['CHD']/f['CHN']
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
        f = np.load(os.path.join('data', exp, 'time_series_phytoplankton.npz'))
        r = f['PHD']/f['PHN']
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
    fig.savefig('ratio_of_mean_n_chl_ts.pdf', dpi=300)
    plt.close(fig)


def plot_zoo_timeseries() -> None:
    """Plot the timeseries of ratios for multiple years and months.

    Parameters
    ----------
    varname : str
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
    fig.set_size_inches(3*w, h*1)
    fig.clf()
    gs: mgs.GridSpec = mgs.GridSpec(
        1, 3, figure=fig,
        wspace=0.4, hspace=0., left=0.065, right=0.935,
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
            'data', 'free', 'time_series_zooplankton.npz'))
        f = np.load(os.path.join('data', exp, 'time_series_zooplankton.npz'))
        r: np.ndarray = f['ZME']/f['ZMI'] - f_free['ZME']/f_free['ZMI']
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
    ax.set_ylabel('Differences from freerun')
    ax.set_title(r'$ZME/ZMI$ ratio')
    line, = ax1.plot(t, f_free['ZME']/f_free['ZMI'], color=colours[0],
                     linestyle=linestyles[0],
                     label=config.exp_labels['free'], alpha=1)

    ax = fig.add_subplot(gs[1])
    ax1 = ax.twinx()
    for exp, linestyle, colour in zip(
            exps[1:],
            linestyles[1:],
            colours[1:]):
        f_free = np.load(os.path.join(
            'data', 'free', 'time_series_zooplankton.npz'))
        f = np.load(os.path.join('data', exp, 'time_series_zooplankton.npz'))
        r = f['ZME'] - f_free['ZME']
        line, = ax.plot(t, r, color=colour,
                        linestyle=linestyle,
                        label=config.exp_labels[exp], alpha=1)

    # Set the major locator to be every day
    ax.xaxis.set_major_locator(locator)
    # Set the major formatter to display the date in 'Month-Day' format
    ax.xaxis.set_major_formatter(mdates.AutoDateFormatter(locator))
    ax.tick_params(axis='x', rotation=30)
    ax.set_xlabel('Time')
    ax.set_title('ZME')
    line, = ax1.plot(t, f_free['ZME'], color=colours[0],
                     linestyle=linestyles[0],
                     label=config.exp_labels['free'], alpha=1)

    ax = fig.add_subplot(gs[2])
    ax1 = ax.twinx()
    for exp, linestyle, colour in zip(
            exps[1:],
            linestyles[1:],
            colours[1:]):
        f_free = np.load(os.path.join(
            'data', 'free', 'time_series_zooplankton.npz'))
        f = np.load(os.path.join('data', exp, 'time_series_zooplankton.npz'))
        r = f['ZMI'] - f_free['ZMI']
        line, = ax.plot(t, r, color=colour,
                        linestyle=linestyle,
                        label=config.exp_labels[exp], alpha=1)

    # Set the major locator to be every day
    ax.xaxis.set_major_locator(locator)
    # Set the major formatter to display the date in 'Month-Day' format
    ax.xaxis.set_major_formatter(mdates.AutoDateFormatter(locator))
    ax.tick_params(axis='x', rotation=30)
    ax.set_xlabel('Time')
    ax.set_title('ZMI')
    line, = ax1.plot(t, f_free['ZMI'], color=colours[0],
                     linestyle=linestyles[0],
                     label=config.exp_labels['free'], alpha=1)
    ax1.set_ylabel('Freerun')

    fig.legend(loc='outside lower center',
               handles=[mlines.Line2D(
                   [0],
                   [0],
                   linewidth=1,  linestyle=ls, color=colour)
                   for ls, colour in zip(linestyles, colours)],
               labels=[config.exp_labels[exp] for exp in exps],
               ncols=7, fontsize=12, markerscale=0.5)
    fig.savefig('ratio_of_mean_zoo_ts.pdf', dpi=300)
    plt.close(fig)


def plot_axes_n_chl(fig, ax, expname: str) -> None:
    """plot the axes"""
    lons, lats = utils.get_coord()
    ax.set_global()
    f = np.load(os.path.join('data', expname, 'rom_n_chl_map.npz'))
    r = f['n']/f['chl']
    norm = mcolors.Normalize(vmin=3., vmax=30)
    cmap = 'cmo.amp'
    if expname == 'free':
        norm = mcolors.SymLogNorm(  # pylint: disable=unexpected-keyword-arg
            linthresh=30., vmin=-80, vmax=80)
        cmap = 'cmo.balance'
        f = np.load(os.path.join('data', 'obs', 'rom_n_chl_map.npz'))
        r = r - f['n']/f['chl']
    elif expname != 'obs':
        norm = mcolors.SymLogNorm(  # pylint: disable=unexpected-keyword-arg
            linthresh=2., vmin=-30, vmax=30)
        cmap = 'cmo.balance'
        f = np.load(os.path.join('data', 'free', 'rom_n_chl_map.npz'))
        r = r - f['n']/f['chl']
    print(np.nanmin(r), np.nanmax(r))
    pc = ax.pcolormesh(
        lons, lats, np.squeeze(r), transform=ccrs.PlateCarree(),
        cmap=cmap, norm=norm)
    ax.set_title(config.exp_labels[expname])
    ax.coastlines(color='k', linewidth=.8)
    ax.add_feature(cfeature.LAND, zorder=3)
    fig.colorbar(pc, ax=ax, orientation='horizontal', pad=0.05,
                 shrink=0.5, norm=norm)


def plot_axes_phytoplanktons(
        fig, ax, expname: str, vname: str) -> None:
    """plot the axes"""
    lons, lats = utils.get_coord()
    ax.set_global()
    f = np.load(os.path.join('data', expname, 'rom_map_phytoplankton.npz'))
    r = f[f'{vname}HD']/f[f'{vname}HN']
    norm = mcolors.TwoSlopeNorm(vmin=0., vcenter=1., vmax=2)
    cmap = 'cmo.balance'
    if expname != 'free':
        if vname == 'P':
            norm = mcolors.SymLogNorm(  # pylint: disable=unexpected-keyword-arg
                linthresh=0.1, vmin=-0.3, vmax=0.3)
        else:
            norm = mcolors.SymLogNorm(  # pylint: disable=unexpected-keyword-arg
                linthresh=0.01, vmin=-2, vmax=2)
        cmap = 'cmo.balance'
        f = np.load(os.path.join('data', 'free', 'rom_map_phytoplankton.npz'))
        r = r - f[f'{vname}HD']/f[f'{vname}HN']
    print(np.nanmin(r), np.nanmax(r))
    pc = ax.pcolormesh(
        lons, lats, np.squeeze(r), transform=ccrs.PlateCarree(),
        cmap=cmap, norm=norm)
    ax.set_title(config.exp_labels[expname])
    ax.coastlines(color='k', linewidth=.8)
    ax.add_feature(cfeature.LAND, zorder=3)
    fig.colorbar(pc, ax=ax, orientation='horizontal', pad=0.05,
                 shrink=0.5, norm=norm)


def plot_axes_zooplanktons_ratio(
        fig, ax, expname: str) -> None:
    """plot the axes"""
    lons, lats = utils.get_coord()
    ax.set_global()
    f = np.load(os.path.join('data', expname, 'rom_map_zooplankton.npz'))
    r = f['ZME']/f['ZMI']
    norm = mcolors.TwoSlopeNorm(vcenter=1,  vmin=0., vmax=3.)
    cmap = 'cmo.balance'
    if expname != 'free':
        norm = mcolors.SymLogNorm(  # pylint: disable=unexpected-keyword-arg
            linthresh=0.03, vmin=-0.3, vmax=0.3)
        cmap = 'cmo.balance'
        f = np.load(os.path.join('data', 'free', 'rom_map_zooplankton.npz'))
        r = r - f['ZME']/f['ZMI']
    print(np.nanmin(r), np.nanmax(r))
    pc = ax.pcolormesh(
        lons, lats, np.squeeze(r), transform=ccrs.PlateCarree(),
        cmap=cmap, norm=norm)
    ax.set_title(config.exp_labels[expname])
    ax.coastlines(color='k', linewidth=.8)
    ax.add_feature(cfeature.LAND, zorder=3)
    fig.colorbar(pc, ax=ax, orientation='horizontal', pad=0.05,
                 shrink=0.5, norm=norm)


def plot_axes_zooplanktons_ZME(
        fig, ax, expname: str) -> None:
    """plot the axes"""
    lons, lats = utils.get_coord()
    ax.set_global()
    f = np.load(os.path.join('data', expname, 'rom_map_zooplankton.npz'))
    r = f['ZME']
    norm = mcolors.Normalize(vmin=0., vmax=1.)
    cmap = 'cmo.amp'
    if expname != 'free':
        norm = mcolors.SymLogNorm(  # pylint: disable=unexpected-keyword-arg
            linthresh=0.005, vmin=-0.04, vmax=0.04)
        cmap = 'cmo.balance'
        f = np.load(os.path.join('data', 'free', 'rom_map_zooplankton.npz'))
        r = r - f['ZME']
    print(np.nanmin(r), np.nanmax(r))
    pc = ax.pcolormesh(
        lons, lats, np.squeeze(r), transform=ccrs.PlateCarree(),
        cmap=cmap, norm=norm)
    ax.set_title(config.exp_labels[expname])
    ax.coastlines(color='k', linewidth=.8)
    ax.add_feature(cfeature.LAND, zorder=3)
    fig.colorbar(pc, ax=ax, orientation='horizontal', pad=0.05,
                 shrink=0.5, norm=norm)


def plot_axes_zooplanktons_ZMI(
        fig, ax, expname: str) -> None:
    """plot the axes"""
    lons, lats = utils.get_coord()
    ax.set_global()
    f = np.load(os.path.join('data', expname, 'rom_map_zooplankton.npz'))
    r = f['ZMI']
    norm = mcolors.Normalize(vmin=0., vmax=0.8)
    cmap = 'cmo.amp'
    if expname != 'free':
        norm = mcolors.SymLogNorm(  # pylint: disable=unexpected-keyword-arg
            linthresh=0.005, vmin=-0.1, vmax=0.1)
        cmap = 'cmo.balance'
        f = np.load(os.path.join('data', 'free', 'rom_map_zooplankton.npz'))
        r = r - f['ZMI']
    print(np.nanmin(r), np.nanmax(r))
    pc = ax.pcolormesh(
        lons, lats, np.squeeze(r), transform=ccrs.PlateCarree(),
        cmap=cmap, norm=norm)
    ax.set_title(config.exp_labels[expname])
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
            if vname == 'zoo':
                plot_axes_zooplanktons_ratio(fig, ax, exp)
            elif vname == 'ZME':
                plot_axes_zooplanktons_ZME(fig, ax, exp)
            elif vname == 'ZMI':
                plot_axes_zooplanktons_ZMI(fig, ax, exp)
            else:
                plot_axes_phytoplanktons(fig, ax, exp, vname)

    fig.savefig(f'figs/ratio_of_mean_map_{vname}.png', dpi=300)


if __name__ == '__main__':
    expnames: list[str] = ['free', 'chlo', 'chlo-monthly',
                           'chlo-monthly-update', 'PC',
                           'PC-update', 'chlo-pc',]
    # for expname in expnames:
    #     config.exp = expname
    #     config.output_path = config.output_path_format.format(exp=expname)
    #     get_n_chl_series(False)
    #     get_phytoplankton_series()
    #     get_n_chl_map(False)
    #     get_phytoplankton_map()
    #     get_zooplankton_series()
    #     get_zooplankton_map()

    # config.exp = ''
    # config.output_path = config.output_path_format.format(exp='')
    # get_n_chl_series(True)
    # get_n_chl_map(True)

    # plot_phyto_timeseries()

    # plot_map('n_chl')
    # plot_map('P')
    # plot_map('C')

    # # plot_zoo_timeseries()
    plot_map('zoo')
    # plot_map('ZME')
    # plot_map('ZMI')
