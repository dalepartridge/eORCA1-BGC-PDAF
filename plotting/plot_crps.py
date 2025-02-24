"""computing CRPS
"""
import multiprocessing as mp
import itertools
import os

import numpy as np
import matplotlib.ticker as mticker
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import matplotlib.gridspec as mgs
import pyPDAF.PDAF as PDAF

import config
import diag
import file_info

plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']
plt.rcParams['font.size'] = 16


def diag_crps_model(year: str, month: str, varname: str) -> None:
    """Compute CRPS from monthly model output and observations.

    Parameters
    ----------
    year : str
        The year of the model output.
    month : str
        The month of the model output.
    varname : str
        The variable name of the model output.
    """
    model_data: np.ma.MaskedArray = np.ma.zeros(
        (config.ne, config.ny, config.nx), dtype=float)

    f_info: file_info.FileInfo
    file_format: str = str(config.varinfo[varname]['file_format'])
    f_info = file_info.FileInfo(year, month, '01',
                                file_format
                                )
    for i in range(config.ne):
        fname = f_info.get_nemo_filename().format(str(i+1))
        model_data[i] = diag.read_model_output_3d(
            fname, 1, np.array([0]), config.varinfo[varname])

    model_data = model_data.reshape(config.ne, config.ny * config.nx)
    model_data = model_data.T

    mask = np.any(model_data.mask, axis=1)

    # get observations
    obs_type: str
    if varname == 'chlo-monthly':
        obs_type = 'chlo-monthly'
    elif varname == 'nitrogen-monthly':
        obs_type = 'pc'
    else:
        raise ValueError(f'Unknown observation type name: {obs_type}')

    obs: np.ma.MaskedArray
    f_info = file_info.FileInfo(year, month, '01', '')
    fname = f_info.get_obsfilename(obs_type)
    obs, _ = diag.read_physical_observation(fname, obs_type)

    obs = obs.reshape(config.ny * config.nx)
    mask = np.logical_or(mask, obs.mask)
    mask = np.logical_or(mask, np.isnan(obs))
    mask = np.logical_or(mask, np.isinf(obs))
    obs = obs[~mask]
    model_data = model_data[~mask, :]

    crps, reli, resol, uncert, status = PDAF.diag_CRPS_nompi(
        0, model_data, obs)

    assert status == 0, 'CRPS calculation failed'

    np.savez(
        os.path.join(
            'crps',
            config.exp,
            f'{varname}_{year}{month}_model.npz'
        ),
        CRPS=crps, reli=reli, resol=resol, uncert=uncert)


def save_crps_multiprocessed(varname: str) -> None:
    """Save the CRPS of the model output
    """

    os.makedirs(os.path.join('crps', config.exp), exist_ok=True)

    years: list[str]
    months: list[str]
    years = ['2015', '2016']
    months = [str(month).zfill(2) for month in range(1, 13)]

    yearmonths = itertools.product(years, months)

    processes: list[mp.Process] = []
    process: mp.Process

    for year, month in yearmonths:
        if year == '2015' and month == '01':
            continue

        process = mp.Process(
            target=diag_crps_model,
            args=(year, month, varname))
        processes.append(process)

    # start and join the processes
    processes_batch = itertools.batched(processes, config.n_process)
    for batch in processes_batch:
        for p in batch:
            p.start()
        for p in batch:
            p.join()


def get_crps_series(varname: str) -> dict[str, list[float]]:
    """get crps time series for plotting
    """
    metric: dict[str, list[float]] = {}
    exps: list[str] = ['free', 'chlo', 'chlo-monthly',
                       'chlo-monthly-update', 'PC',
                       'PC-update', 'chlo-pc',
                       ]
    years: list[str] = ['2015', '2016']
    months: list[str] = [str(month).zfill(2) for month in range(1, 13)]

    for exp in exps:
        metric[f'{exp}_CRPS'] = []
        metric[f'{exp}_reli'] = []
        metric[f'{exp}_resol'] = []
        metric[f'{exp}_uncert'] = []

    iterator = itertools.product(exps, years, months)
    for exp, year, month in iterator:
        if year == '2015' and month == '01':
            continue
        f = np.load(
            os.path.join(
                'crps',
                exp,
                f'{varname}_{year}{month}_model.npz')
        )
        metric[f'{exp}_CRPS'].append(f['CRPS'])
        metric[f'{exp}_reli'].append(f['reli'])
        metric[f'{exp}_resol'].append(f['resol'])
        metric[f'{exp}_uncert'].append(f['uncert'])

    return metric


def plot_crps_series(varname: str) -> None:
    """Plot the time series of crps.
    """

    metric = get_crps_series(varname)

    exps: list[str] = ['free', 'chlo', 'chlo-monthly',
                       'chlo-monthly-update', 'PC',
                       'PC-update', 'chlo-pc',
                       ]
    linestyles: list[str] = [':', '-', '-', '-', '-', '-', '-', ]
    colours: list[str] = ['r', '#1E88E5', '#FFC107', '#48B03E',
                          'k', '#2CF8BA', '#D81B1B']
    metric_suffixes: list[str] = ['CRPS', 'reli', 'resol', 'uncert']
    metric_names: list[str] = ['CRPS', 'reliability',
                               'potential CRPS', 'uncertainty']

    # time array
    t: np.ndarray
    t = np.arange('2015-02', '2016-01', dtype='datetime64[M]')

    fig: plt.Figure = plt.figure()
    gs: mgs.GridSpec = mgs.GridSpec(
        2, 2, figure=fig,
        wspace=0.12, hspace=0.1, left=0.07, right=0.98,
        bottom=0.13, top=0.96)
    w, h = fig.get_size_inches()
    fig.set_size_inches(w * 2, h * 2)

    for i, (suffix, label) in enumerate(zip(metric_suffixes,
                                            metric_names)):
        ax: plt.Axes = fig.add_subplot(gs[i])
        lines = []
        for exp, linestyle, colour in zip(exps, linestyles, colours):
            line, = ax.plot(
                t, metric[f'{exp}_{suffix}'],
                linestyle=linestyle, marker='.',
                color=colour, label=label)
            lines.append(line)
            ax.set_title(f'{label} of {varname}')

        if i >= 2:
            # Set the major locator to be every day
            locator = mdates.MonthLocator(bymonth=range(2, 13, 2))
            ax.xaxis.set_major_locator(locator)

            # Set the major formatter to display the date in 'Month-Day' format
            ax.xaxis.set_major_formatter(
                mdates.AutoDateFormatter(locator))
            ax.tick_params(axis='x', rotation=15)
        else:
            ax.xaxis.set_major_locator(mticker.NullLocator())

    fig.legend(loc='outside lower center',
               handles=lines,
               labels=[config.exp_labels[exp] for exp in exps],
               ncols=7, fontsize=12)
    fig.savefig(f'CRPS_{varname}.png', dpi=300)
    plt.close(fig)


if __name__ == '__main__':
    expnames: list[str] = ['free', 'chlo', 'chlo-monthly',
                           'chlo-monthly-update', 'PC',
                           'PC-update', 'chlo-pc',]
    for expname in expnames:
        config.exp = expname
        config.output_path = config.output_path_format.format(exp=expname)
        for vname in ['chlo-monthly', 'nitrogen-monthly']:
            save_crps_multiprocessed(vname)

    for vname in ['chlo-monthly', 'nitrogen-monthly']:
        plot_crps_series(vname)
