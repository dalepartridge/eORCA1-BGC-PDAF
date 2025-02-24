import datetime
import multiprocessing as mp
import itertools
import os

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as mgs
import matplotlib.ticker as mticker

import config
import diag
import file_info


plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']
plt.rcParams['font.size'] = 16


def save_acc(year: str, month: str, varname: str) -> None:
    """Calculate anomaly correlation coefficient.
    """
    # Load the climatology data
    f = np.load('plotting/5year_clim.npz')
    clim: np.ma.MaskedArray
    if varname == 'chlo-monthly':
        clim = np.ma.array(
            f['clim_CHD'] + f['clim_CHN'],
            mask=(f['CHD_mask']) & (f['CHN_mask']))
        obs_type = 'chlo-monthly'

    elif varname == 'nitrogen-monthly':
        clim = np.ma.array(
            f['clim_PHD'] + f['clim_PHN'],
            mask=(f['PHD_mask']) & (f['PHN_mask']))
        obs_type = 'pc'
    else:
        raise ValueError(f'Unknown varname: {varname}')

    # get observation anomaly
    y: np.ma.MaskedArray
    f_info = file_info.FileInfo(year, month, '01', '')
    fname = f_info.get_obsfilename(obs_type)
    y, _ = diag.read_physical_observation(fname, obs_type)
    y_anom = y - clim
    # get model anomaly
    file_format: str = str(config.varinfo[varname]['file_format'])
    f_info = file_info.FileInfo(year, month, '01',
                                file_format
                                )

    acc = np.zeros(config.ne)
    for i in range(config.ne):
        fname = f_info.get_nemo_filename().format(str(i + 1))
        x = diag.read_model_output_3d(
            fname, 1, np.array([0]), config.varinfo[varname])
        x_anom = np.ma.masked_where(x < 0., x)
        x_anom = x_anom - clim
        # Calculate the anomaly correlation coefficient
        acc[i] = np.ma.corrcoef(x_anom.ravel(), y_anom.ravel())[0, 1]
    # Save monthly bias to a numpy archive
    np.savez(f'acc/{config.exp}/{varname}_{year}{month}.npz', acc=acc)


def save_acc_multiprocessed(varname: str) -> None:
    """Save monthly bias for multiple years and months using multiprocessing.

    Args:
        None

    Returns:
        None
    """
    os.makedirs(f'acc/{config.exp}', exist_ok=True)

    years: list[str] = ['2015', ]
    months: list[str] = [str(month).zfill(2) for month in range(1, 13)]

    yearmonths = itertools.product(years, months)

    processes: list[mp.Process] = []

    for year, month in yearmonths:
        if year == '2015' and month == '01':
            continue

        # Calculate the number of days in the month
        process: mp.Process = mp.Process(
            target=save_acc, args=(year, month, varname))
        processes.append(process)

    for process in processes:
        process.start()

    for process in processes:
        process.join()


def plot_acc_series(varnames: list[str]) -> None:
    """Plot time series of the acc
    """
    exps: list[str] = ['PC', 'chlo', 'chlo-monthly',
                       'chlo-monthly-update', 'chlo-pc']
    linestyles: list[str] = ['-', ':', '-', '-', '-']
    colours: list[str] = ['k', 'k', '#1E88E5', '#FFC107',
                          '#004D40', 'r']
    varnames_fig: list[str] = [
        'nitrogen' if varname == 'nitrogen-monthly'
        else 'chlorophyll' for varname in varnames]
    # time array
    t: np.ndarray
    t = np.arange('2015-02', '2016-01', dtype='datetime64[M]')
    # time
    years: list[str]
    months: list[str]
    years = ['2015', ]
    months = [str(month).zfill(2) for month in range(1, 13)]
    # ACC data
    acc: dict[str, list[float]] = {}

    for exp in exps:
        for varname in varnames:
            acc[f'{varname}_{exp}'] = []
            yearmonths = itertools.product(years, months)
            for year, month in yearmonths:
                if year == '2015' and month == '01':
                    continue
                x = np.load(f'acc/{exp}/{varname}_{year}{month}.npz',
                            )['acc'][:]
                acc[f'{varname}_{exp}'].append(x)

    # get the figure
    fig: plt.Figure = plt.figure(1)
    fig.clf()
    # expand the figure width
    w, h = fig.get_size_inches()
    fig.set_size_inches(w * 2, h)
    # set the grid for each subplot
    gs = mgs.GridSpec(1, 2, figure=fig, wspace=0.08,
                      left=0.03, right=0.99, bottom=0.19, top=0.95)
    ax: plt.Axes

    for i, (varname, varname_fig) in enumerate(
            zip(varnames, varnames_fig)):
        ax = fig.add_subplot(gs[i])
        # plot the ACC time series for each experiment
        lines: list = []
        for exp, linestyle, colour in zip(exps, linestyles, colours):
            print(exp, varname, np.array(acc[f'{varname}_{exp}']).shape)
            line, = ax.plot(
                t, np.array(acc[f'{varname}_{exp}']).mean(1),
                color=colour, linestyle=linestyle)  # type: ignore
            lines.append(line)
            ax.fill_between(
                t, np.array(acc[f'{varname}_{exp}']).mean(1) -
                0.5*np.array(acc[f'{varname}_{exp}']).std(1),
                np.array(acc[f'{varname}_{exp}']).mean(1) +
                0.5*np.array(acc[f'{varname}_{exp}']).std(1),
                color=colour, alpha=0.1)  # type: ignore
        ax.set_title(varname_fig)
        # ax.set_ylim((-0.05, 0.7))

    fig.autofmt_xdate()
    fig.legend(
        loc='outside lower center', handles=lines,
        labels=[config.exp_labels[exp] for exp in exps],
        ncols=5)
    fig.savefig('acc_ts.png', dpi=300)
    plt.close(fig)


def plot_acc_boxplots(years: list[str], varnames: list[str]) -> None:
    """Plot ACC boxplots
    """
    exps: list[str] = ['chlo', 'PC', 'chlo-monthly',
                       'chlo-monthly-update', 'free', 'chlo-pc']

    acc: dict = dict()
    for exp in exps:
        for varname in varnames:
            acc[f'{varname}_{exp}'] = []

    t: list[datetime.date] = []  # Timestamps

    for year in years:
        for EXP in ['chlo', 'PC', 'chlo-monthly', 'free', 'chlo-pc']:
            for varname in varnames:
                x = np.load(f'{EXP}/acc_ens_{varname}_{year}{month}.npz',
                            )['acc'][:, 0, 1]
                acc[f'{varname}_{EXP}'].append(x)

    fig: plt.Figure = plt.figure(1)
    fig.clf()
    w, h = fig.get_size_inches()
    fig.set_size_inches(w * 2, h)
    gs = mgs.GridSpec(1, 2, figure=fig, wspace=0.1, left=0.04,
                      right=0.99, bottom=0.2, top=0.93)
    ax: plt.Axes

    # set fixed locators for x ticks
    loc = mticker.FixedLocator([-0.2, 2.3, 4.8])
    # set fixed formatter for x ticks
    fmt = mticker.FixedFormatter(t)  # type: ignore

    for i, varname in enumerate(varnames):
        ax = fig.add_subplot(gs[i])
        pos = 2.5*np.arange(3) - 0.8
        for j, (EXP, linestyle, color) in enumerate(
            zip(['PC', 'free', 'chlo', 'chlo-monthly', 'chlo-pc'],
                ['--', '-', ':', ':', '-'],
                ['gray', 'r', 'r', 'k', 'k']
                )
        ):
            medianprops = dict(linestyle=linestyle, linewidth=2.5, color=color)
            positions = pos + j*0.3
            ax.boxplot(np.array(acc[f'{varname}_{EXP}']).T,
                       medianprops=medianprops,
                       widths=0.2,
                       positions=np.round(positions, 3))  # type: ignore
        ax.xaxis.set_major_locator(loc)
        ax.xaxis.set_major_formatter(fmt)
        ax.set_title(f'EXP: ACC {varname}')

    # set legends
    fig.legend(loc='outside lower center',
               labels=['Carbon', 'Freerun', 'daily Chlo',
                       'monthly Chlo', 'monthly Chlo & Carbon'],
               handles=[plt.Line2D([0], [0], color='gray', linestyle='--'),
                        plt.Line2D([0], [0], color='r', linestyle='-'),
                        plt.Line2D([0], [0], color='r', linestyle=':'),
                        plt.Line2D([0], [0], color='k', linestyle=':'),
                        plt.Line2D([0], [0], color='k', linestyle='-')],
               ncols=5)

    fig.savefig('ACC_box.png', dpi=300)
    plt.close(fig)


if __name__ == '__main__':
    # for exp in ['chlo', 'PC', 'chlo-monthly', 'chlo-pc', 'free',
    #             'chlo-monthly-update']:
    #     config.exp = exp
    #     config.output_path = config.output_path_format.format(exp=exp)
    #     save_acc_multiprocessed('nitrogen-monthly')
    #     save_acc_multiprocessed('chlo-monthly')

    plot_acc_series(['chlo-monthly', 'nitrogen-monthly', ])
    # plot_acc_boxplots(['2015',], ['chlorophyll', 'nitrogen',])
