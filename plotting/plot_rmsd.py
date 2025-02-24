"""Comaprison of experiments in RMSD sense.
"""
import os
import typing

import numpy as np
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import matplotlib.gridspec as mgs
import matplotlib.dates as mdates

import config
import diag
import file_info
import plot_crps

plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman'] + \
    plt.rcParams['font.serif']
plt.rcParams['font.size'] = 16


def save_monthly(filename_model: str,
                 obs_type: typing.Literal['chlo-monthly', 'pc']
                 ) -> None:
    """Save monthly RMSE over experiment period

    Parameters
    ----------
    filename_model : str
        filename of the monthly climatology in .npz file.
    obs_type : Literal['chlo-monthly', 'pc']
        name of the variable used to determine the observations.
    """
    assert obs_type in ['pc', 'chlo-monthly'], \
        'Currently, we only calculate RMSEs ' \
        'using `pc` or `chlo-monthly` observations'

    years: list[str] = ['2015', '2016']
    months: list[str] = [str(month).zfill(2)
                         for month in range(1, 13)]

    rmse: list[float] = []
    for year in years:
        for month in months:
            if year == '2015' and month == '01':
                continue

            # Read the physical observation
            f_info: file_info.FileInfo
            f_info = file_info.FileInfo(year, month, '01', '')
            obs, _ = diag.read_physical_observation(
                f_info.get_obsfilename(obs_type), obs_type)

            # read the monthly model ensemble mean
            varname = 'nitrogen' if obs_type == 'pc' \
                else 'chlo'
            f: np.lib.npyio.NpzFile = np.load(
                filename_model.format(year=year, month=month,
                                      varname=varname
                                      )
            )
            model: np.ma.MaskedArray = np.ma.masked_array(
                f['model'], f['model_mask']
            )

            # Calculate RMSE
            err2 = (obs - model)*(obs - model)
            rmse.append(np.sqrt(err2.mean()))

    # Save monthly bias to a numpy archive
    np.savez(
        os.path.join('data', config.exp,
                     f'rmse_{varname}.npz'
                     ),
        RMSE=np.array(rmse)
    )


def plot_timeseries() -> None:
    """Plot the RMSD timeseries for multiple years and months.

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
    t = np.arange('2015-02', '2017-01', dtype='datetime64[M]')
    # making the plot
    fig: plt.Figure = plt.figure()
    # increase the width of the figure as we will have two subplots
    w: float
    h: float
    w, h = fig.get_size_inches()
    fig.set_size_inches(2*w, h*2)
    fig.clf()
    gs: mgs.GridSpec = mgs.GridSpec(
        2, 2, figure=fig,
        wspace=0.43, hspace=0.4, left=0.09, right=0.92,
        bottom=0.14, top=0.96)

    locator: mdates.MonthLocator
    locator = mdates.MonthLocator(bymonth=range(2, 13, 3))

    # loop over the subplots
    for i, varname in enumerate(['chlo', 'nitrogen']):
        ax: plt.Axes = fig.add_subplot(gs[i])
        # providing a second y-axis
        ax1: plt.Axes = ax.twinx()
        # plot difference between  for each experiment
        ax.axhline(y=0, color='lightgrey')
        lines: list[mlines.Line2D] = []
        rmsd_free: np.ndarray = np.load(
            f'data/free/rmse_{varname}.npz')['RMSE']

        for exp, linestyle, colour in zip(
                exps[1:],
                linestyles[1:],
                colours[1:]):
            rmsd: np.ndarray = np.load(
                f'data/{exp}/rmse_{varname}.npz')['RMSE']
            line: mlines.Line2D
            line, = ax.plot(t, (rmsd - rmsd_free)/rmsd_free, color=colour,
                            linestyle=linestyle,
                            label=config.exp_labels[exp], alpha=1)
            lines.append(line)

        # Set the major locator to be every day
        ax.xaxis.set_major_locator(locator)
        # Set the major formatter to display the date in 'Month-Day' format
        ax.xaxis.set_major_formatter(mdates.AutoDateFormatter(locator))
        ax.tick_params(axis='x', rotation=20)
        if i == 0:
            ax.set_ylabel('normalised RMSD difference')
        ax.set_xlabel('Time')
        vname = ''
        if varname == 'chlo':
            vname = 'Chl'
        if varname == 'nitrogen':
            vname = 'N'
        # plot the freerun RMSD as a reference
        line, = ax1.plot(
            t, rmsd_free,
            color='r', linestyle=':', label='Freerun', alpha=1)
        ax1.tick_params(axis='y', labelcolor='r')
        if i == 0:
            ax.set_title(f'a) phytoplankton {vname}')
            ax1.set_ylabel(r'Freerun RMSD (mg C m$^{-3}$)', color='r')
        if i == 1:
            ax.set_title(f'b) phytoplankton {vname}')
            ax1.set_ylabel(r'Freerun RMSD (mmol N m$^{-3}$)', color='r')
        lines.append(line)

    for i, varname in enumerate(['chlo-monthly', 'nitrogen-monthly']):
        metric = plot_crps.get_crps_series(varname)
        ax = fig.add_subplot(gs[2 + i])
        for exp, linestyle, colour in zip(exps, linestyles, colours):
            ax.plot(
                t, metric[f'{exp}_reli'],
                linestyle=linestyle, marker='.',
                color=colour)
            if varname == 'chlo-monthly':
                vname = 'Chl'
            if varname == 'nitrogen-monthly':
                vname = 'N'
            ax.set_title(f'Reliability of {vname}')

        ax.set_ylabel('Reliability')
        # Set the major locator to be every day
        ax.xaxis.set_major_locator(locator)
        # Set the major formatter to display the date in 'Month-Day' format
        ax.xaxis.set_major_formatter(mdates.AutoDateFormatter(locator))
        ax.tick_params(axis='x', rotation=20)
        ax.set_xlabel('Time')

    fig.legend(loc='outside lower center',
               handles=[mlines.Line2D(
                   [0],
                   [0],
                   linewidth=1,  linestyle=ls, color=colour)
                   for ls, colour in zip(linestyles, colours)],
               labels=[config.exp_labels[exp] for exp in exps],
               ncols=7, fontsize=12, markerscale=0.5)
    fig.savefig('RMSD_reliability_timeseries.pdf', dpi=300)
    plt.close(fig)


if __name__ == '__main__':
    expnames: list[str] = ['free', 'chlo', 'chlo-monthly',
                           'chlo-monthly-update', 'PC',
                           'PC-update', 'chlo-pc',]

    for expname in expnames:
        config.exp = expname
        obs_types: list[typing.Literal['chlo-monthly', 'pc']] = [
            'chlo-monthly', 'pc']
        for obs_t in obs_types:
            save_monthly(os.path.join(
                'data',
                config.exp, 'ensmean_{varname}-monthly_{year}{month}01.npz'),
                obs_t)

    plot_timeseries()
