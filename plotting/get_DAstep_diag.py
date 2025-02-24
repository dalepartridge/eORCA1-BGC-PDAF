import itertools
import os

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as mgs
import matplotlib.colors as mcolors
import matplotlib.lines as mlines
import cartopy  # type: ignore
import cartopy.crs as ccrs  # type: ignore
import cmocean  # type: ignore # pylint: disable=unused-import

import config
import diag
import file_info
import utils


def plot_increment(varname: str, exp: str, label: str) -> None:
    """Plot the increment map."""

    os.makedirs(f'increment/{exp}', exist_ok=True)

    years: list[str] = ['2015', ]
    months: list[str] = [str(month).zfill(2)
                         for month in range(1, 13)]

    lons: np.ma.MaskedArray
    lats: np.ma.MaskedArray
    lons, lats = utils.get_coord()

    fname: str = os.path.join(
        'data',
        exp,
        f'pdaf_ensmean_{varname}_''{year}{month}{day}.npz'
    )

    yearmonths = itertools.product(years, months)

    for year, month in yearmonths:
        if year == '2015' and month == '01':
            continue

        # Get the number of days in the month
        ndays: int = utils.get_month_days(int(year), int(month))
        days: range = range(ndays//2, ndays//2 + 1)
        if exp == 'chlo':
            days = range(1, ndays + 1)

        # Loop over days
        for day in days:
            day_str = str(day).zfill(2)
            fig = plt.figure()
            fig.clf()
            w, h = fig.get_size_inches()
            fig.set_size_inches(w, h)
            gs = mgs.GridSpec(1, 1, figure=fig, wspace=0.01,
                              left=0., right=1., bottom=0.0, top=1.)
            ax: cartopy.mpl.geoaxes.GeoAxes = fig.add_subplot(
                gs[0], projection=ccrs.Robinson())
            ax.set_global()
            # Read the innovation

            f = np.load(fname.format(year=year, month=month,
                                     day=day_str))

            pdaf_data: np.ma.MaskedArray = np.ma.masked_array(
                f['model'], f['model_mask'])

            pc = ax.pcolormesh(
                lons, lats, pdaf_data[1] - pdaf_data[0],
                cmap='cmo.balance',
                norm=mcolors.CenteredNorm(vcenter=0, halfrange=0.1),
                transform=ccrs.PlateCarree())
            ax.coastlines(color='k', linewidth=.8)
            fig.colorbar(pc, ax=ax, orientation='horizontal',
                         shrink=0.6, pad=0.02)
            ax.set_title(f'{label} observation - forecast')

            fig.savefig(
                os.path.join(
                    'increment',
                    exp,
                    f'{varname}_{year}_{month}{day_str}.png'),
                dpi=300)
            plt.close(fig)


def plot_individual_histogram(varnames: list[str]) -> None:
    """Plot the histogram of variables in PDAF output.

    Parameters
    ----------
    varname: str
        The variable name.
    """
    assert len(varnames) == 2, 'Only two variables are plotted.'

    years: list[str] = ['2015', '2016']
    months: list[str] = [str(month).zfill(
        2) for month in range(1, 13)]

    years = ['2015', ]

    exps = ['PC', 'chlo-monthly',
            'chlo-monthly-update', 'chlo-pc']
    linestyles = ['-', '-', ':', ':']
    colours = ['k', '#FFC107', '#004D40', 'r']

    # create the figure
    fig = plt.figure()
    fig.clf()
    w, h = fig.get_size_inches()
    fig.set_size_inches(w*2, h*2)

    gs = mgs.GridSpec(2, 2, figure=fig, wspace=0.17,
                      left=0.08, right=0.99, bottom=0.25, top=0.94)

    varname_data = itertools.product(varnames, ['obs', 'forecasts'])

    for i, (varname, data_type) in enumerate(varname_data):
        ax = fig.add_subplot(gs[i])
        data: dict[str, np.ndarray] = dict()
        exp_yearmonths = itertools.product(exps, years, months)

        # initialise the bias for each exp for varname
        for exp in exps:
            data[exp] = np.array([])

        # get an array of o - b
        for exp, year, month in exp_yearmonths:

            config.exp = exp
            config.output_path = config.output_path_format.format(exp=exp)
            if year == '2015' and month == '01':
                continue

            if varname == 'nitrogen' and exp in ['chlo',
                                                 'chlo-monthly',
                                                 'chlo-monthly-update'
                                                 ]:
                continue

            if varname == 'chlorophyll' and exp in ['PC',]:
                continue

            obs_type = 'pc'
            if varname == 'chlorophyll':
                obs_type = 'chlo-monthly'
                if exp == 'chlo':
                    obs_type = 'chlo'

            # Get the number of days in the month
            ndays: int = utils.get_month_days(int(year), int(month))
            days: range = range(ndays//2, ndays//2 + 1)
            if exp == 'chlo':
                days = range(1, ndays + 1)

            # Loop over days
            for day in days:
                day_str = str(day).zfill(2)
                if data_type == 'obs':
                    f_info = file_info.FileInfo(year, month,
                                                day_str, '')
                    obs, _ = diag.read_observation(
                        f_info.get_obsfilename(obs_type), obs_type)

                    data[exp] = np.append(
                        data[exp], obs.filled(np.nan).ravel())
                else:
                    f_info = file_info.FileInfo(year, month,
                                                day_str, '')
                    vname_file = 'carbon' if exp == 'PC' else varname
                    fname = f_info.get_pdaf_filename(vname_file)
                    for i in range(1, config.ne + 1):
                        model_data = diag.read_pdaf_output_3d(
                            fname.format(str(i).zfill(3)),
                            vname_file) / config.ne
                        data[exp] = np.append(
                            data[exp],
                            model_data.filled(np.nan).ravel()
                        )

        # plot the histogram for current varname
        for exp, linestyle, colour in zip(exps, linestyles, colours):
            print(varname, exp, data[exp].shape)
            ax.hist(data[exp], bins=100,
                    density=True,
                    histtype='step', color=colour,
                    linestyle=linestyle,
                    linewidth=3, label=exp)

        if varname == 'nitrogen':
            ax.set_xlabel(data_type + r' ($mmol/m^3$)')
        if varname == 'chlorophyll':
            ax.set_xlabel(data_type + r' ($mg/m^3$)')

        ax.set_title(varname)
        ax.set_ylabel('Frequency')

    lines = [mlines.Line2D(
        [], [], color=color, linestyle=linestyle)
        for linestyle, color in zip(
        linestyles, colours
    )
    ]

    fig.legend(
        loc='outside lower center', handles=lines,
        labels=[config.exp_labels[exp] for exp in exps],
        ncols=5)

    fig.savefig('hist.pdf', dpi=300)
    plt.close(fig)


def plot_distance_physical(varname: str, exp: str, label: str) -> None:
    """Plot the distance between pdaf output and obs.
    in physical space.
    """
    config.exp = exp
    config.output_path = \
        config.output_path_format.format(exp=exp)

    if varname == 'nitrogen' and exp in ['chlo',
                                         'chlo-monthly',
                                         'chlo-monthly-update'
                                         ]:
        return

    if varname == 'chlorophyll' and exp in ['PC',]:
        return

    obs_type = 'pc'
    if varname == 'chlorophyll':
        obs_type = 'chlo-monthly'
        if exp == 'chlo':
            obs_type = 'chlo'

    os.makedirs(f'monthly_physical_distance/{exp}', exist_ok=True)

    years: list[str] = ['2015', ]
    months: list[str] = [str(month).zfill(2) for month in range(1, 13)]

    lons: np.ma.MaskedArray
    lats: np.ma.MaskedArray
    lons, lats = utils.get_coord()

    yearmonths = itertools.product(years, months)
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

            fig = plt.figure()
            fig.clf()
            w, h = fig.get_size_inches()
            fig.set_size_inches(w, h)
            gs = mgs.GridSpec(1, 1, figure=fig, wspace=0.01,
                              left=0., right=1., bottom=0.0,
                              top=1.)
            ax: cartopy.mpl.geoaxes.GeoAxes = fig.add_subplot(
                gs[0], projection=ccrs.Robinson())
            ax.set_global()
            # Read the physical observation
            f_info: file_info.FileInfo
            f_info = file_info.FileInfo(year, month, '01', '')
            obs_type = 'pc' if varname == 'nitrogen' \
                else 'chlo-monthly'
            obs, _ = diag.read_physical_observation(
                f_info.get_obsfilename(obs_type), obs_type)

            # read the monthly model ensemble mean
            f_info = file_info.FileInfo(year, month,
                                        day_str, '')
            obs, _ = diag.read_physical_observation(
                f_info.get_obsfilename(obs_type), obs_type)

            f_info = file_info.FileInfo(year, month,
                                        day_str, '')
            vname_file = 'carbon' if exp == 'PC' else varname
            fname = f_info.get_pdaf_filename(vname_file)
            model_data = diag.get_physical_ensemble_mean_from_pdaf(
                fname, vname_file, f_info.stats_pdaf)

            pc = ax.pcolormesh(
                lons, lats, np.squeeze(obs - model_data[0]),
                cmap='cmo.balance',
                norm=mcolors.CenteredNorm(vcenter=0, halfrange=0.3),
                transform=ccrs.PlateCarree())
            ax.coastlines(color='k', linewidth=.8)
            fig.colorbar(pc, ax=ax, orientation='horizontal',
                         shrink=0.6, pad=0.02)
            ax.set_title(f'{label} observation - forecast')
            fig.savefig(
                os.path.join(
                    'monthly_physical_distance',
                    exp,
                    f'{varname}_{year}_{month}.png'),
                dpi=300)
            plt.close(fig)


def plot_error_correlation(
        year: str, month: str, day: str, varname_0: str,
        varname_1: str):
    lons, lats = utils.get_coord()

    f_info = file_info.FileInfo(year, month, str(day).zfill(2), '')
    ens_anomaly_0 = diag.get_ensemble_anomaly_from_pdaf(
        f_info.get_pdaf_filename(varname_0),
        varname_0, f_info.stats_pdaf)[:, 0]
    ens_anomaly_1 = diag.get_ensemble_anomaly_from_pdaf(
        f_info.get_pdaf_filename(varname_1),
        varname_1, f_info.stats_pdaf)[:, 0]

    # compute the correlation matrix
    corr = np.sum(ens_anomaly_0*ens_anomaly_1, axis=0)/(config.ne-1)
    corr = corr/np.std(ens_anomaly_0, axis=0,
                       ddof=1)/np.std(ens_anomaly_1, axis=0, ddof=1)

    print(corr.min(), corr.max())
    fig = plt.figure()
    fig.clf()
    gs = mgs.GridSpec(1, 1, figure=fig, wspace=0.01,
                      left=0., right=1., bottom=0.0, top=1.)
    ax: cartopy.mpl.geoaxes.GeoAxes = fig.add_subplot(
        gs[0], projection=ccrs.Robinson())
    ax.set_global()
    pc = ax.pcolormesh(lons, lats, corr, cmap='cmo.balance',
                       norm=mcolors.CenteredNorm(vcenter=0, halfrange=1),
                       transform=ccrs.PlateCarree())
    ax.coastlines(color='k', linewidth=.8)
    fig.colorbar(pc, ax=ax, orientation='horizontal',
                 shrink=0.6, pad=0.02)

    ax.set_title('ensemble cross-correlation')
    fig.savefig(f'cross_corr_{config.exp}_{year}{month}.png', dpi=300)
    plt.close(fig)


def plot_obs_err(year: str, month: str, day: str):
    lons, lats = utils.get_coord()
    f_info = file_info.FileInfo(year, month, day, '')
    obs, obs_unc = diag.read_observation(
        f_info.get_obsfilename('chlo-monthly'), 'chlo-monthly')

    fig = plt.figure()
    fig.clf()
    w, h = fig.get_size_inches()
    fig.set_size_inches(2*w, h)
    gs = mgs.GridSpec(1, 2, figure=fig, wspace=0.01,
                      left=0., right=1., bottom=0.0, top=1.)
    ax: cartopy.mpl.geoaxes.GeoAxes = fig.add_subplot(
        gs[0], projection=ccrs.Robinson())
    ax.set_global()
    pc = ax.pcolormesh(lons, lats, obs, cmap='cmo.matter',
                       # norm = mcolors.Normalize(vmin=0, vmax=0.5),
                       transform=ccrs.PlateCarree())
    ax.coastlines(color='k', linewidth=.8)
    fig.colorbar(pc, ax=ax, orientation='horizontal',
                 shrink=0.6, pad=0.02)
    ax.set_title('Obs. Error of monthly chlorophyll')

    f_info = file_info.FileInfo(year, month, day, '')
    obs, obs_unc = diag.read_observation(
        f_info.get_obsfilename('pc'), 'pc')
    ax = fig.add_subplot(gs[1], projection=ccrs.Robinson())
    ax.set_global()
    pc = ax.pcolormesh(lons, lats, obs, cmap='cmo.matter',
                       # norm = mcolors.Normalize(vmin=0, vmax=0.5),
                       transform=ccrs.PlateCarree())
    ax.coastlines(color='k', linewidth=.8)
    fig.colorbar(pc, ax=ax, orientation='horizontal',
                 shrink=0.6, pad=0.02)

    print(obs_unc[obs < 1e14].mean())
    ax.set_title(f'Obs. Error of nitrogen {obs_unc.mean()}')
    fig.savefig(f'test_obs_{year}{month}.png', dpi=300)
    plt.close(fig)


def check_obs_mask() -> None:
    lons, lats = utils.get_coord()
    months = ['01', '02', '03', '04', '05',
              '06', '07', '08', '09', '10', '11', '12']

    yearmonths = itertools.product(['2015', '2016'], months)
    for year, month in yearmonths:
        if year == '2015' and month == '01':
            continue
        f_info = file_info.FileInfo(year, month, '01', '')
        obs_physc, obs_physc_unc = diag.read_physical_observation(
            f_info.get_obsfilename('pc'), 'pc')
        obs, obs_unc = diag.read_observation(
            f_info.get_obsfilename('pc'), 'pc')

        fig = plt.figure()
        fig.clf()
        w, h = fig.get_size_inches()
        fig.set_size_inches(2*w, 2*h)
        gs = mgs.GridSpec(2, 2, figure=fig, wspace=0.01,
                          left=0., right=1., bottom=0.0, top=1.)
        ax: cartopy.mpl.geoaxes.GeoAxes = fig.add_subplot(
            gs[0], projection=ccrs.Robinson())
        ax.set_global()
        pc = ax.pcolormesh(lons, lats, obs, cmap='cmo.matter',
                           transform=ccrs.PlateCarree())
        ax.coastlines(color='k', linewidth=.8)
        fig.colorbar(pc, ax=ax, orientation='horizontal',
                     shrink=0.6, pad=0.02)
        ax.set_title('obs log10')

        ax = fig.add_subplot(
            gs[1], projection=ccrs.Robinson())
        ax.set_global()
        pc = ax.pcolormesh(lons, lats, obs_unc, cmap='cmo.matter',
                           transform=ccrs.PlateCarree())
        ax.coastlines(color='k', linewidth=.8)
        fig.colorbar(pc, ax=ax, orientation='horizontal',
                     shrink=0.6, pad=0.02)
        ax.set_title('Obs. Error of monthly nitrogen')

        ax = fig.add_subplot(
            gs[2], projection=ccrs.Robinson())
        ax.set_global()
        pc = ax.pcolormesh(lons, lats, obs_physc, cmap='cmo.matter',
                           transform=ccrs.PlateCarree())
        ax.coastlines(color='k', linewidth=.8)
        fig.colorbar(pc, ax=ax, orientation='horizontal',
                     shrink=0.6, pad=0.02)
        ax.set_title('physical obs')

        ax = fig.add_subplot(
            gs[3], projection=ccrs.Robinson())
        ax.set_global()
        pc = ax.pcolormesh(lons, lats, obs_physc_unc, cmap='cmo.matter',
                           transform=ccrs.PlateCarree())
        ax.coastlines(color='k', linewidth=.8)
        fig.colorbar(pc, ax=ax, orientation='horizontal',
                     shrink=0.6, pad=0.02)
        ax.set_title('physical obs err')

        fig.savefig(f'test_obs_new_{year}{month}.png', dpi=300)
        plt.close(fig)


if __name__ == '__main__':
    for exp in ['PC-update', ]:
        config.exp = exp
        plot_increment('nitrogen', exp, config.exp_labels[exp])

    # for exp in ['chlo', 'chlo-monthly', 'chlo-monthly-update']:
    #     config.exp = exp
    #     plot_increment('chlorophyll', exp, config.exp_labels[exp])
    # plot_individual_histogram(['chlorophyll', 'nitrogen',])

    # plot_distance_physical('nitrogen', 'PC', config.exp_labels['PC'])
