"""obtain ensemble mean from model and pdaf output files.
"""
import itertools
import multiprocessing as mp
import os
import typing

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.gridspec as mgs
import matplotlib.colors as mcolors
import cartopy  # type: ignore
import cartopy.crs as ccrs  # type: ignore
import cmocean  # type: ignore # pylint: disable=unused-import

import config
import utils
import diag
from file_info import FileInfo
import operators


def save_ensemble_mean_model(
        f_info: FileInfo, fname_output: str,
        level: np.ndarray, varname: str) -> None:
    """Save ensemble mean to an .npz file that is of given `varname`
    at given vertical `level` from NEMO output.

    Parameters
    ----------
    f_info : file_info.FileInfo
        file information.
    fname_output : str
        output filename.
    day : int
        day of the model data.
    level : np.ndarray
        vertical levels of the model data.
    varname : str
        name of the variable.
    """
    # Calculate/read the ensemble mean from the model data file
    model_data: np.ma.MaskedArray
    model_data = diag.get_model_output_ensemble_mean(
        f_info.get_nemo_filename(),
        f_info.stats, int(f_info.day), level, config.varinfo[varname])

    # Save the monthly climatology to a numpy archive
    print(fname_output.format(year=f_info.year, month=f_info.month,
                              day=f_info.day))
    np.savez(
        fname_output.format(year=f_info.year, month=f_info.month,
                            day=f_info.day),
        model=model_data.data,
        model_mask=model_data.mask
    )
    print(f'saved {fname_output.format(
        year=f_info.year, month=f_info.month, day=f_info.day)}')


def save_ensemble_mode_model(
        f_info: FileInfo, fname_output: str, day: int,
        level: np.ndarray, varname: str) -> None:
    """Save mode of ensemble distribution to an .npz file that is of
    given `varname` at given vertical `level` from NEMO output.

    This is calculated based on the assumption that the ensemble
    follows a log-normal distribution.

    Parameters
    ----------
    f_info : file_info.FileInfo
        file information.
    fname_output : str
        output filename.
    month : str
        month of the model data.
    day : int
        day of the model data.
    level : np.ndarray
        vertical levels of the model data.
    varname : str
        name of the variable.
    """
    # varname must contain -log10 because the mode of
    # the lognormal distribution is calculated by its
    # parameters.
    assert '-log10' in varname, "varname must contain -log10"
    # Calculate the ensemble mean and standard deviation
    # from the model data
    model_data: np.ma.MaskedArray
    model_unc: np.ma.MaskedArray
    model_data = diag.get_model_output_ensemble_mean(
        f_info.get_nemo_filename(), f_info.stats,
        day, level, config.varinfo[varname])
    model_unc = diag.get_model_output_std(f_info.get_nemo_filename(),
                                          f_info.stats, day, level,
                                          config.varinfo[varname])
    # get the lognormal parameters
    mu: np.ma.MaskedArray
    std: np.ma.MaskedArray
    mu, std = operators.get_lognormal_params(model_data, model_unc)
    mode: np.ma.MaskedArray = operators.get_lognormal_mode(mu, std)

    # Save the monthly climatology to a numpy archive
    # We don't want to keep -log10 in the filename because only physical
    np.savez(
        fname_output,
        model=mode.data,
        model_mask=mode.mask
    )


def save_ensemble_mean_multiprocessed(
        varname: str, level: np.ndarray,
        target_func: typing.Callable) -> None:
    """A generic wrapper function used to perform multi-processing
    with daily data.

    Parameters
    ----------
    varname : str
        Name of the variable.
    level : np.ndarray
        Vertical levels
    target_func : typing.Callable
        The target function to be executed.
        This can be :func:`save_ensemble_mean_model` or
        :func:`save_ensemble_mode_model`.
    do_monthly : bool
        Whether to process monthly data.
    """

    # Define the years to process
    years: list[str] = ['2015', '2016']

    # Define the cycling months to process
    months: list[str] = [str(month).zfill(2)
                         for month in range(1, 13, 1)]

    # Create the output directory if it doesn't exist
    os.makedirs(os.path.join('data', config.exp), exist_ok=True)

    # pylint: disable=comparison-with-callable
    if target_func == save_ensemble_mean_model:
        enstype = 'ensmean'
    else:
        enstype = 'ensmode'
    # Specify output filename
    fname_output: str = os.path.join(
        'data', config.exp,
        f'{enstype}_{varname}''_{year}{month}{day}.npz'
    )

    # Determine the file format
    file_format: str = str(config.varinfo[varname]['file_format'])
    do_monthly: bool = '_1m_' in file_format

    yearmonths = itertools.product(years, months)
    processes: list[mp.Process] = []

    # Loop over years
    for year, month in yearmonths:
        if year == '2015' and month == '01':
            continue
        # Get the number of days in the month
        days: range
        if do_monthly:
            days = range(1, 2)
        else:
            days = range(1, utils.get_month_days(
                int(year),
                int(month)) + 1)
        # Loop over days
        for day in days:
            f_info = FileInfo(year, month, str(day).zfill(2),
                              file_format
                              )

            # Create a process for saving spread
            process: mp.Process = mp.Process(
                target=target_func,
                args=(f_info, fname_output, level, varname))
            processes.append(process)

    # start and join the processes
    processes_batch = itertools.batched(processes, config.n_process)
    for batch in processes_batch:
        for p in batch:
            p.start()
        for p in batch:
            p.join()


def save_ensemble_mean_pdaf(
        f_info: FileInfo, fname_output: str, varname: str) -> None:
    """Save ensemble mean of pdaf variable.

    Parameters
    ----------
    f_info : file_info.FileInfo
        file information.
    fname_output : str
        output filename.
    varname : str
        Name of the variable.
    """
   # Calculate the ensemble mean for the model data
    model_data: np.ma.MaskedArray = diag.get_ensemble_mean_from_pdaf(
        f_info.get_pdaf_filename(varname), varname, f_info.stats_pdaf)

    # Save the monthly climatology to a numpy archive
    np.savez(
        fname_output.format(year=f_info.year, month=f_info.month,
                            day=f_info.day),
        model=model_data.data,
        model_mask=model_data.mask
    )


def save_ensemble_mean_pdaf_multiprocessed(varname: str) -> None:
    """Save ensemble mean of PDAF output through
    multiple years and months using multiprocessing.

    Parameters
    ----------
    varname : str
        Name of the variable.
    """

    # Define the years to process
    years: list[str] = ['2015', '2016']

    # Define the cycling months to process
    months: list[str] = [str(month).zfill(2)
                         for month in range(1, 13)]

    # Create the output directory if it doesn't exist
    os.makedirs(os.path.join('data',
                             config.exp), exist_ok=True)

    varname_output: str = 'nitrogen' if varname == 'carbon' \
        else varname
    fname_output: str = os.path.join(
        'data',
        config.exp,
        f'pdaf_ensmean_{varname_output}_''{year}{month}{day}.npz'
    )

    yearmonths = itertools.product(years, months)
    processes: list[mp.Process] = []
    # Loop over years
    for year, month in yearmonths:
        if year == '2015' and month == '01':
            continue

        # Get the number of days in the month
        ndays = utils.get_month_days(int(year), int(month))
        days: range = range(ndays//2, ndays//2 + 1)
        if config.exp == 'chlo':
            days = range(1, ndays + 1)

        # Loop over days
        for day in days:
            f_info = FileInfo(year, month, str(day).zfill(2), '')
            # Create a process for saving spread
            process: mp.Process = mp.Process(
                target=save_ensemble_mean_pdaf,
                args=(f_info, fname_output, varname))
            processes.append(process)

    # start and join the processes
    processes_batch = itertools.batched(processes, config.n_process)
    for batch in processes_batch:
        for p in batch:
            p.start()
        for p in batch:
            p.join()


def plot_monthly_ensmean_from_free(
        exp: str, label: str, varnames: list[str]) -> None:
    """Plot the monthly bias for multiple months in a single year.

    Args:
        year: PDAF output year.
        month: PDAF output month.
    """
    lons: np.ma.MaskedArray
    lats: np.ma.MaskedArray
    mask: np.ma.MaskedArray = utils.get_land_mask()
    lons, lats = utils.get_coord()

    years: list[str] = ['2015', ]
    months: list[str] = [str(month).zfill(2) for month in range(1, 13)]
    # months = ['07',]
    vmax = {'FER': 9e-6, 'DIN': 0.3, 'DIC': 1.5, 'ALK': 0.15, 'SIL': 0.2,
            'DET': 0.15, 'ZMI': 0.01, 'ZME': 0.02, 'OXY': 1.,
            'CHN': 0.1, 'CHD': 0.1, 'PHN': 0.1, 'PHD': 0.1,
            'PRN': 1., 'PRD': 0.3, 'PD_JLIM': 1., 'PN_JLIM': 1.2,
            'GMEPN': 0.45, 'GMEPD': 1.0,
            'PN_LLOSS': 0.15, 'PD_LLOSS': 0.3,
            'MPN': 1.0, 'MPD': 1.0, 'GMIPn': 1.0, 'MED_XPAR': 1.0,
            'Rn': 1.5, 'Rd': 1.5}
    os.makedirs(f'phyto_{exp}_diff', exist_ok=True)
    for year in years:
        for month in months:
            if year == '2015' and month == '01':
                continue
            print(year, month)
            fig: matplotlib.figure.Figure = plt.figure()
            fig.clf()
            w, h = fig.get_size_inches()
            fig.set_size_inches(w*2, h*2)
            gs = mgs.GridSpec(2, 2, figure=fig, wspace=0.01,
                              left=0., right=1., bottom=0.0, top=1.)
            for i, varname in enumerate(varnames):
                ax: cartopy.mpl.geoaxes.GeoAxes = fig.add_subplot(
                    gs[i], projection=ccrs.Robinson())
                f = np.load(
                    f'{exp} /ensmean_{varname} _{year} {month} 01.npz',
                    allow_pickle=True)
                var: np.ma.MaskedArray = np.ma.masked_array(
                    f['model'], f['model_mask'] | mask)
                f_free = np.load(
                    f'free/ensmean_{varname} _{year} {month} 01.npz',
                    allow_pickle=True)
                var_free: np.ma.MaskedArray = np.ma.masked_array(
                    f_free['model'], f_free['model_mask'] | mask)
                diff = np.squeeze(var - var_free)
                diff = np.ma.masked_array(diff, diff.mask | (diff > 0))
                ax.set_global()
                print(varname, np.nanpercentile(
                    np.abs(diff.filled(np.nan)), 90))
                pc = ax.pcolormesh(
                    lons, lats, diff, cmap="cmo.balance",
                    norm=mcolors.CenteredNorm(
                        vcenter=0, halfrange=vmax[varname]),
                    transform=ccrs.PlateCarree())
                ax.coastlines(color='k', linewidth=.8)
                ax.set_title(f'{label}: {varname} in {year}-{month}')
                fig.colorbar(pc, ax=ax, orientation='horizontal',
                             shrink=0.6, pad=0.02)

            fig.savefig(f'phyto_{exp}_diff/phyto_{year}{month}.png', dpi=300)
            plt.close(fig)


def plot_monthly_ensmean_free(
        exp: str, label: str, varnames: list[str]) -> None:
    """Plot the monthly bias for multiple months in a single year.

    Args:
        year: PDAF output year.
        month: PDAF output month.
    """
    lons: np.ma.MaskedArray
    lats: np.ma.MaskedArray
    mask: np.ma.MaskedArray = utils.get_land_mask()
    lons, lats = utils.get_coord()

    years: list[str] = ['2015', ]
    months: list[str] = [str(month).zfill(2) for month in range(1, 13)]
    vmin = {'DIN': 0.3, 'DIC': 1970, 'SIL': 0.4, 'FER': 2e-4,
            'ALK': 2200, 'DET': 2e-3, 'ZMI': 1e-3, 'ZME': 0.,
            'CHN': 0., 'CHD': 0., 'PHN': 0., 'PHD': 0.,
            'PRN': 1., 'PRD': 0.3, 'PD_JLIM': 1., 'PN_JLIM': 1.2,
            'GMEPN': 0.45, 'GMEPD': 1.0,
            'PN_LLOSS': 0.15, 'PD_LLOSS': 0.3,
            'MPN': 1.0, 'MPD': 1.0, 'GMIPn': 1.0, 'MED_XPAR': 1.0,
            'Rn': 0., 'Rd': 0.}
    vmax = {'DIN': 30, 'DIC': 2200, 'SIL': 70, 'FER': 1e-3,
            'ALK': 2400, 'DET': 0.25, 'ZMI': 0.25, 'ZME': 0.45,
            'CHN': 5., 'CHD': 5., 'PHN': 5., 'PHD': 5.,
            'PRN': 1., 'PRD': 0.3, 'PD_JLIM': 1., 'PN_JLIM': 1.2,
            'GMEPN': 0.45, 'GMEPD': 1.0,
            'PN_LLOSS': 0.15, 'PD_LLOSS': 0.3,
            'MPN': 1.0, 'MPD': 1.0, 'GMIPn': 1.0, 'MED_XPAR': 1.0,
            'Rn': 1.2, 'Rd': 1.2}
    os.makedirs(f'phyto_{exp}', exist_ok=True)
    for year in years:
        for month in months:
            if year == '2015' and month == '01':
                continue
            print(year, month)
            fig: matplotlib.figure.Figure = plt.figure()
            fig.clf()
            w, h = fig.get_size_inches()
            fig.set_size_inches(w*2, h*1)
            gs = mgs.GridSpec(1, 2, figure=fig, wspace=0.01,
                              left=0., right=1., bottom=0.0, top=1.)
            for i, varname in enumerate(varnames):
                ax: cartopy.mpl.geoaxes.GeoAxes = fig.add_subplot(
                    gs[i], projection=ccrs.Robinson())
                f = np.load(
                    f'{exp} /ensmean_{varname} _{year} {month} 01.npz',
                    allow_pickle=True)
                clim: np.ma.MaskedArray = np.ma.masked_array(
                    np.squeeze(f['model']),
                    np.squeeze(f['model_mask']) | mask |
                    (np.squeeze(f['model']) >= 1.2))
                ax.set_global()
                # print (varname, np.nanpercentile(np.abs(clim.filled(np.nan)), 10),
                # np.nanpercentile(np.abs(clim.filled(np.nan)), 90))
                print(varname, np.nanpercentile(clim.filled(np.nan), 10),
                      np.nanpercentile(clim.filled(np.nan), 90))
                pc = ax.pcolormesh(lons, lats, clim, cmap="cmo.matter",
                                   norm=mcolors.Normalize(vmin=vmin[varname],
                                                          vmax=vmax[varname]
                                                          ),
                                   transform=ccrs.PlateCarree())
                ax.coastlines(color='k', linewidth=.8)
                ax.set_title(f'{label}: {varname} in {year}-{month}')
                fig.colorbar(pc, ax=ax, orientation='horizontal',
                             shrink=0.6, pad=0.02)

            fig.savefig(f'phyto_{exp}/R_{year}{month}.png', dpi=300)
            plt.close(fig)


if __name__ == '__main__':
    # directory name of experiments
    expnames: list[str] = ['free', 'chlo', 'chlo-monthly',
                           'chlo-monthly-update', 'PC',
                           'PC-update', 'chlo-pc',
                           ]
    # get ensemble mean in model daily output files
    variables: list[str]
    variables = ['chlorophyll', 'nitrogen', ]
    variables = ['chlo-monthly', 'nitrogen-monthly', ]
    variables = ['CHN', 'CHD', 'PHN', 'PHD', 'ZMI', 'ZME']

    # variables = ['PRN', 'PRD']
    variables = ['TCO2', 'CO2FLUX', 'ATM_PCO2', 'OCN_PCO2', 'O2FLUX', 'OXY', ]
    variables = ['PRN', 'PRD', ]  # ['PD_JLIM', 'PN_JLIM']
    variables = ['DET', ]
    # variables = ['PD_LLOSS', 'PN_LLOSS']
    # variables = ['MPN', 'MPD', 'GMIPn', 'GMEPN', 'GMEPD', 'PD_JLIM',
    #              'PD_NLIM', 'PD_FELIM', 'PD_SILIM', 'PDSILIM2', 'PN_JLIM',
    #              'PN_NLIM', 'PN_FELIM']

    for expname in expnames:
        config.exp = expname
        config.output_path = config.output_path_format.format(exp=expname)
        for vname in variables:
            save_ensemble_mean_multiprocessed(vname, np.array([0]),
                                              save_ensemble_mean_model
                                              )

    # # get ensemble mean in PDAF output files
    # for expname in expnames[1:]:
    #     config.exp = expname
    #     config.output_path = config.output_path_format.format(exp=expname)
    #     if expname in ['PC', 'PC-update']:
    #         save_ensemble_mean_pdaf_multiprocessed('carbon')
    #     if expname in ['chlo-pc', ]:
    #         save_ensemble_mean_pdaf_multiprocessed('nitrogen')
    #     if expname in ['chlo', 'chlo-monthly', 'chlo-pc',
    #                    'chlo-monthly-update']:
    #         save_ensemble_mean_pdaf_multiprocessed('chlorophyll')
