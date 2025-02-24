"""Diagnosing the spread in DA step and model output
"""
import itertools
import multiprocessing as mp
import os

import numpy as np
import matplotlib.colors as mcolors
import matplotlib.dates as mdates
import matplotlib.gridspec as mgs
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import cartopy  # type: ignore
import cartopy.crs as ccrs  # type: ignore
import cmocean  # type: ignore # pylint: disable=unused-import

import config
import diag
import file_info
import utils


plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']
plt.rcParams['font.size'] = 16


def get_model_timeseries(varname: str, exp: str) -> np.ndarray:
    """get the monthly model spread timeseries
    """
    config.exp = exp
    config.output_path = config.output_path_format.format(exp=exp)
    if varname == 'chlorophyll':
        varname = 'chlo-monthly'
    if varname == 'nitrogen':
        varname = 'nitrogen-monthly'

    years: list[str] = ['2015', '2016']
    months: list[str] = [str(month).zfill(2)
                         for month in range(1, 13)]

    years = ['2015', ]

    yearmonth = itertools.product(years, months)
    spread: np.ndarray = np.zeros(len(years)*len(months) - 1)
    i: int = 0
    for year, month in yearmonth:
        if year == '2015' and month == '01':
            continue

        # read the monthly model ensemble mean
        f_info: file_info.FileInfo
        f_info = file_info.FileInfo(year, month, '01', str(
            config.varinfo[varname]['file_format']))

        spread[i] = diag.get_model_output_std(f_info.get_nemo_filename(),
                                              False, 1, np.array([0]),
                                              config.varinfo[varname]
                                              ).mean()
        i += 1

    return spread


def get_pdaf(f_info: file_info.FileInfo, exp: str,
             varname: str) -> np.ma.MaskedArray:
    """get pdaf spread data
    """
    vname_file: str = 'carbon' if 'PC' in exp else varname
    fname = f_info.get_pdaf_filename(vname_file)
    return diag.get_pdaf_output_std(fname, vname_file,
                                    f_info.stats_pdaf)


def get_pdaf_val(f_info: file_info.FileInfo, exp: str,
                 varname: str) -> np.ma.MaskedArray:
    """get pdaf ensemble mean data
    """
    vname_file: str = 'carbon' if 'PC' in exp else varname
    fname = f_info.get_pdaf_filename(vname_file)
    return diag.get_ensemble_mean_from_pdaf(fname, vname_file,
                                            f_info.stats_pdaf)


def get_pdaf_physic(f_info: file_info.FileInfo, exp: str,
                    varname: str) -> np.ma.MaskedArray:
    """get spread in physical space
    """
    vname_file: str = 'carbon' if 'PC' in exp else varname
    fname = f_info.get_pdaf_filename(vname_file)
    return diag.get_pdaf_output_physical_std(fname, vname_file,
                                             f_info.stats_pdaf)


def get_pdaf_physicval(f_info: file_info.FileInfo, exp: str,
                       varname: str) -> np.ma.MaskedArray:
    """get ensemble mean in physical space
    """
    vname_file: str = 'carbon' if 'PC' in exp else varname
    fname = f_info.get_pdaf_filename(vname_file)
    return diag.get_physical_ensemble_mean_from_pdaf(fname, vname_file,
                                                     f_info.stats_pdaf)


def get_obs(f_info: file_info.FileInfo, obs_type: str
            ) -> np.ma.MaskedArray:
    """get obs err data
    """
    fname = f_info.get_obsfilename(obs_type)
    std_obs: np.ma.MaskedArray
    _, std_obs = diag.read_observation(fname, obs_type)
    return std_obs


def plot_forecast_obs(f_info_pdaf: file_info.FileInfo,
                      exp: str, varname: str,
                      obs_type: str) -> None:
    """Plot the forecast and observation spread
    """
    lons, lats = utils.get_coord()
    fig: plt.Figure = plt.figure()
    fig.clf()
    w, h = fig.get_size_inches()
    fig.set_size_inches(w * 3, h * 1)
    gs: mgs.GridSpec = mgs.GridSpec(1, 3, figure=fig)
    ax: cartopy.mpl.geoaxes.GeoAxes = fig.add_subplot(
        gs[0], projection=ccrs.Robinson())
    ax.set_global()
    pc = ax.pcolormesh(
        lons, lats, get_pdaf(f_info_pdaf, exp, varname)[0],
        cmap='cmo.amp',
        norm=mcolors.Normalize(vmin=0, vmax=0.2),
        transform=ccrs.PlateCarree())
    ax.coastlines(color='k', linewidth=.8)
    fig.colorbar(pc, ax=ax, orientation='horizontal',
                 shrink=0.6, pad=0.02)
    ax.set_title('forecast spread')

    ax = fig.add_subplot(
        gs[1], projection=ccrs.Robinson())
    ax.set_global()
    pc = ax.pcolormesh(
        lons, lats, get_obs(f_info_pdaf, obs_type),
        cmap='cmo.amp',
        norm=mcolors.Normalize(vmin=0, vmax=0.5),
        transform=ccrs.PlateCarree())
    ax.coastlines(color='k', linewidth=.8)
    fig.colorbar(pc, ax=ax, orientation='horizontal',
                 shrink=0.6, pad=0.02)
    ax.set_title('obs err')

    ax = fig.add_subplot(
        gs[2], projection=ccrs.Robinson())
    ax.set_global()
    pc = ax.pcolormesh(
        lons, lats, get_pdaf(f_info_pdaf, exp, varname)[0] -
        get_obs(f_info_pdaf, obs_type),
        cmap='cmo.balance', norm=mcolors.CenteredNorm(
            vcenter=0, halfrange=0.5),
        transform=ccrs.PlateCarree())
    ax.coastlines(color='k', linewidth=.8)
    fig.colorbar(pc, ax=ax, orientation='horizontal',
                 shrink=0.6, pad=0.02)
    ax.set_title('obs err')

    fname: str = f'{varname}_{f_info_pdaf.year}{f_info_pdaf.month}.png'
    fig.savefig(os.path.join('figs', 'spread_DA', exp, fname), dpi=300)
    plt.close(fig)


def plot_increment(f_info: file_info.FileInfo,
                   exp: str, varname: str) -> None:
    """Plot the increment of analysis and forecast spread
    """
    lons, lats = utils.get_coord()
    fig: plt.Figure = plt.figure()
    fig.clf()
    w, h = fig.get_size_inches()
    fig.set_size_inches(w * 2, h * 1)
    gs: mgs.GridSpec = mgs.GridSpec(1, 2, figure=fig,
                                    wspace=0,
                                    left=0., right=1,
                                    bottom=0., top=0.93)
    ax: cartopy.mpl.geoaxes.GeoAxes = fig.add_subplot(
        gs[0], projection=ccrs.Robinson())
    ax.set_global()
    spread = get_pdaf_physic(f_info, exp, varname)
    dspread = spread[1] - spread[0]
    print(np.abs(dspread).min(), np.abs(dspread).max(),
          np.nanpercentile(np.abs(dspread.filled(np.nan)), 95))
    pc = ax.pcolormesh(
        lons, lats, dspread,
        cmap='cmo.balance',
        norm=mcolors.CenteredNorm(vcenter=0, halfrange=0.05),
        transform=ccrs.PlateCarree())
    val = get_pdaf_physicval(f_info, exp, varname)
    mask = val[1] - val[0] > 0
    ax.plot(lons[mask], lats[mask],
            'o', markersize=0.1, color='k', alpha=0.5,
            transform=ccrs.PlateCarree())
    ax.coastlines(color='k', linewidth=.8)
    fig.colorbar(pc, ax=ax, orientation='horizontal',
                 shrink=0.6, pad=0.02)
    ax.set_title('physical space')

    ax = fig.add_subplot(
        gs[1], projection=ccrs.Robinson())
    ax.set_global()
    spread = get_pdaf(f_info, exp, varname)
    dspread = spread[1] - spread[0]
    pc = ax.pcolormesh(
        lons, lats, dspread,
        cmap='cmo.balance',
        norm=mcolors.CenteredNorm(vcenter=0, halfrange=0.05),
        transform=ccrs.PlateCarree())
    val = get_pdaf_physicval(f_info, exp, varname)
    mask = val[1] - val[0] > 0
    ax.plot(lons[mask], lats[mask],
            'o', markersize=0.1, color='k', alpha=0.5,
            transform=ccrs.PlateCarree())
    ax.coastlines(color='k', linewidth=.8)
    fig.colorbar(pc, ax=ax, orientation='horizontal',
                 shrink=0.6, pad=0.02)
    ax.set_title('lognormal parameter space')

    fname = f'{varname}_inc_{f_info.year}{f_info.month}{f_info.day}.png'
    if config.exp != 'chlo':
        fname = f'{varname}_inc_{f_info.year}{f_info.month}.png'
    fig.savefig(
        os.path.join('figs', 'spread_DA', exp, fname),
        dpi=300)
    plt.close(fig)


def plot_timeseries() -> None:
    """Plot the spread timeseries for multiple years and months.

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
    linestyles: list[str] = [':', '-', '-', '-', '-', '-', '-', ]
    colours: list[str] = ['r', '#1E88E5', '#FFC107', '#48B03E',
                          'k', '#2CF8BA', '#D81B1B']
    # time array
    t: np.ndarray
    t = np.arange('2015-02', '2016-01', dtype='datetime64[M]')
    # making the plot
    fig: plt.Figure = plt.figure()
    # increase the width of the figure as we will have two subplots
    w: float
    h: float
    w, h = fig.get_size_inches()
    fig.set_size_inches(2*w, h)
    fig.clf()
    gs: mgs.GridSpec = mgs.GridSpec(1, 2, figure=fig,
                                    wspace=0.36,
                                    left=0.09, right=0.93,
                                    bottom=0.3, top=0.93)
    # loop over the subplots
    for i, varname in enumerate(['chlorophyll', 'nitrogen']):
        ax: plt.Axes = fig.add_subplot(gs[i])
        lines: list[mlines.Line2D] = []

        for exp, linestyle, colour in zip(exps, linestyles, colours):
            spread = get_model_timeseries(varname, exp)
            line: mlines.Line2D
            line, = ax.plot(t, spread, color=colour,
                            linestyle=linestyle,
                            label=config.exp_labels[exp], alpha=1)
            lines.append(line)

        # Set the major locator to be every day
        locator: mdates.MonthLocator
        locator = mdates.MonthLocator(bymonth=range(2, 13, 3))
        ax.xaxis.set_major_locator(locator)
        # Set the major formatter to display the date in 'Month-Day' format
        ax.xaxis.set_major_formatter(mdates.AutoDateFormatter(locator))
        ax.tick_params(axis='x', rotation=20)
        if i == 0:
            ax.set_ylabel('normalised RMSD difference')
        ax.set_xlabel('Time')
        if varname == 'chlo':
            varname = 'chlorophyll'
        if i == 0:
            ax.set_title(f'a) phytoplankton {varname}')
        if i == 1:
            ax.set_title(f'b) phytoplankton {varname}')
        lines.append(line)

    fig.legend(loc='outside lower center',
               handles=lines,
               labels=[config.exp_labels[exp] for exp in exps] +
               ['Freerun', ],
               ncols=7, fontsize=12, markerscale=0.5)
    fig.savefig('spread_timeseries.pdf', dpi=300)
    plt.close(fig)


def do_time_loop(exp: str, varname: str, obs_type: str) -> None:
    """Loop over time to plot the spread
    """
    config.exp = exp
    config.output_path = config.output_path_format.format(exp=exp)

    years: list[str]
    months: list[str]
    years = ['2015', '2016']
    months = [str(month).zfill(2) for month in range(1, 13)]
    yearmonths = itertools.product(years, months)

    os.makedirs(os.path.join('spread_DA', exp), exist_ok=True)

    process: list[mp.Process] = []

    for year, month in yearmonths:
        if year == '2015' and month == '01':
            continue

        # Get the number of days in the month
        ndays: int = utils.get_month_days(int(year), int(month))
        days: range = range(ndays//2, ndays//2 + 1)
        if exp == 'chlo':
            days = range(1, ndays + 1)

        for day in days:
            day_str = str(day).zfill(2)
            f_info_pdaf: file_info.FileInfo
            f_info_pdaf = file_info.FileInfo(
                year, month, day_str, exp)

            # p = mp.Process(
            #     target=plot_forecast_obs,
            #     args=(f_info_pdaf, f_info_obs, exp, varname, obs_type)
            # )
            p = mp.Process(
                target=plot_increment,
                args=(f_info_pdaf, exp, varname)
            )
            process.append(p)

    batched_processes = itertools.batched(process, 16)
    for processes in batched_processes:
        for p in processes:
            p.start()
        for p in processes:
            p.join()


def get_onestep_physc(spread, inc, i, f_info: file_info.FileInfo, exp: str,
                      varname: str) -> None:
    """Get the one-step spread and ensemble mean increment in physical space"""
    data = get_pdaf_physic(f_info, exp, varname)
    print ()
    spread[i*config.ny*config.nx:(i+1)*config.ny*config.nx] = \
        (data[1] - data[0]).ravel()
    data = get_pdaf_physicval(f_info, exp, varname)
    inc[i*config.ny*config.nx:(i+1)*config.ny*config.nx] = \
        (data[1] - data[0]).ravel()


def get_onestep(spread, inc, i, f_info: file_info.FileInfo, exp: str,
                varname: str) -> None:
    """Get the one-step spread and ensemble mean increment"""
    data = get_pdaf(f_info, exp, varname)
    print (data[0].shape, config.nx, config.nx)
    spread[i*config.ny*config.nx:(i+1)*config.ny*config.nx] = \
        (data[1] - data[0]).ravel()
    data = get_pdaf_val(f_info, exp, varname)
    inc[i*config.ny*config.nx:(i+1)*config.ny*config.nx] = \
        (data[1] - data[0]).ravel()


def get_all_pdaf(exp: str, varname: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    config.exp = exp
    config.output_path = config.output_path_format.format(exp=exp)

    years: list[str]
    months: list[str]
    years = ['2015', '2016']
    months = [str(month).zfill(2) for month in range(1, 13)]
    yearmonths = itertools.product(years, months)

    nmonths = len(years)*len(months) - 1
    spread = mp.Array('d', np.ma.zeros(nmonths*config.ny*config.nx))
    inc = mp.Array('d', np.ma.zeros(nmonths*config.ny*config.nx))
    spread_physc = mp.Array('d', np.ma.zeros(nmonths*config.ny*config.nx))
    inc_physc = mp.Array('d', np.ma.zeros(nmonths*config.ny*config.nx))

    process: list[mp.Process] = []

    i = 0
    for year, month in yearmonths:
        if year == '2015' and month == '01':
            continue

        # Get the number of days in the month
        ndays: int = utils.get_month_days(int(year), int(month))
        days: range = range(ndays//2, ndays//2 + 1)
        if exp == 'chlo':
            days = range(1, ndays + 1)

        for day in days:
            day_str = str(day).zfill(2)
            f_info: file_info.FileInfo
            f_info = file_info.FileInfo(
                year, month, day_str, exp)

            p = mp.Process(
                target=get_onestep_physc,
                args=(spread_physc, inc_physc, i, f_info, exp, varname)
            )
            process.append(p)

            p = mp.Process(
                target=get_onestep,
                args=(spread, inc, i, f_info, exp, varname)
            )
            process.append(p)
            i = i + 1

    batched_processes = itertools.batched(process, config.n_process)
    for processes in batched_processes:
        for p in processes:
            p.start()
        for p in processes:
            p.join()

    np.savez(
        os.path.join('data', exp, f'spread_inc_scatter_{varname}.npz'),
        spread=np.asarray(spread),
        inc=np.asarray(inc),
        spread_physc=np.asarray(spread_physc),
        inc_physc=np.asarray(inc_physc))

    return np.asarray(spread), np.asarray(inc), np.asarray(spread_physc), np.asarray(
        inc_physc)


def plot_scatter() -> None:
    """Plot the scatter plot of the spread against increments.
    """

    exps: list[str] = ['chlo-monthly',
                       'chlo-monthly-update', 'PC',
                       'PC-update', 'chlo-pc',
                       ]
    colours: list[str] = ['#FFC107', '#48B03E',
                          'k', '#2CF8BA', '#D81B1B']

    fig: plt.Figure = plt.figure()
    fig.clf()
    w, h = fig.get_size_inches()
    fig.set_size_inches(w * 4, h * 2)
    gs: mgs.GridSpec = mgs.GridSpec(2, 4, figure=fig,
                                    wspace=0.25,
                                    hspace=0.27,
                                    left=0.04, right=0.99,
                                    bottom=0.12, top=0.96)

    axes = [fig.add_subplot(gs[i]) for i in range(8)]
    for ax in axes:
        ax.axvline(0, color='gray', linestyle='--')
        ax.axhline(0, color='gray', linestyle='--')
        ax.set_yscale('symlog', linthresh=0.3)
        ax.set_xscale('symlog', linthresh=0.3)
        ax.set_xlabel('Spread increment')
        ax.set_ylabel('Ensemble mean increment')

    for i, varname in enumerate(['chlorophyll', 'nitrogen']):
        vname = ''
        if varname == 'chlorophyll':
            vname = 'Chl'
        if varname == 'nitrogen':
            vname = 'N'
        j = 0
        for exp, colour in zip(exps, colours):
            if exp in ['PC', 'PC-update', ] and varname == 'chlorophyll':
                continue
            if exp in ['chlo-monthly', 'chlo-monthly-update', ] and varname == 'nitrogen':
                continue
            f = np.load(
                os.path.join('data', exp,
                             f'spread_inc_scatter_{varname}.npz')
            )
            print (varname, exp, j)
            axes[4*i + 2*j].scatter(f['spread'], f['inc'], s=0.2, c=colour, alpha=0.1)
            axes[4*i + 2*j].set_title(f'Lognormal parameters ({config.exp_labels[exp]}: {vname})')
            axes[4*i + 2*j + 1].set_title(f'Physical values ({config.exp_labels[exp]}: {vname})')
            axes[4*i + 2*j + 1].scatter(f['spread_physc'],
                                        f['inc_physc'], s=0.2,
                                        c=colour, alpha=0.1)
            if exp in ['chlo-monthly-update', 'PC-update', ]:
                j += 1

    fig.legend(
        loc='outside lower center',
        handles=[mlines.Line2D(
            [0],
            [0],
            linewidth=0, marker='.', markersize=12, color=colour)
            for colour in colours],
        labels=[config.exp_labels[exp] for exp in exps],
        ncols=5, fontsize=12, markerscale=0.5)

    fig.savefig('figs/spread_scatter.png', dpi=300)


if __name__ == '__main__':
    # for expname in ['chlo-monthly', 'chlo-monthly-update', 'chlo-pc',]:
    #     do_time_loop(expname, 'chlorophyll', 'chlo-monthly')

    # for expname in ['PC', 'chlo-pc', 'PC-update',]:
    #     do_time_loop(expname, 'nitrogen', 'pc')

    # for expname in ['chlo', ]:
    #     do_time_loop(expname, 'chlorophyll', 'chlo')

    # plot_timeseries()
    expnames: list[str] = ['chlo-monthly', 'chlo-monthly-update',
                           'PC', 'PC-update', 'chlo-pc',
                          ]
    for expname in expnames:
        config.output_path = config.output_path_format.format(exp=expname)
        if expname in ['chlo-monthly','chlo-monthly-update', 'chlo-pc']:
            get_all_pdaf(expname, 'chlorophyll')
        if expname in ['PC', 'PC-update', 'chlo-pc']:
            get_all_pdaf(expname, 'nitrogen')
    plot_scatter()
