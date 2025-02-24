"""Comaprison of ratio between phytoplankton constituents.
"""
import itertools
import os

import matplotlib.dates as mdates
import matplotlib.gridspec as mgs
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np

import config


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


def plot_timeseries() -> None:
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
    # linestyles: list[str] = ['-', ':', '-', ':', '-', ':', ':', ]
    # colours: list[str] = ['gray', '#1E88E5', '#FFC107', '#FFC107',
    #                       '#2CF8BA', '#2CF8BA', '#D81B1B']
    linestyles: list[str] = [':', '-', '-', '-', '-', '-', '-', ]
    colours: list[str] = ['r', '#1E88E5', '#FFC107', '#48B03E',
                          'k', '#2CF8BA', '#D81B1B']
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
    fig.set_size_inches(w*2, h*1.)
    fig.clf()
    gs: mgs.GridSpec = mgs.GridSpec(
        1, 2, figure=fig,
        wspace=0.45, hspace=0., left=0.07, right=0.915,
        bottom=0.38, top=0.93)

    locator: mdates.MonthLocator
    locator = mdates.MonthLocator(bymonth=range(2, 13, 3))

    # loop over the subplots
    for i, vname in enumerate(['OCN_PCO2', 'OXY']):
        ax: plt.Axes = fig.add_subplot(gs[i])
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

        if i == 0:
            ax1.set_ylabel(r'Freerun pCO$_2$', color='r')
            ax.set_title(r'a) pCO$_2$')
            ax.set_ylabel('Differences from freerun')
        if i == 1:
            ax1.set_ylabel(
                r'Freerun Oxygen ($mmol \cdot O_2 \cdot m^3$)', color='r')
            ax.set_title('b) Oxygen')
            ax.set_ylabel(
                r'Differences from freerun ($mmol \cdot O_2 \cdot m^3$)')
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
               ncols=7, fontsize=12, markerscale=0.5)
    fig.savefig('flux_gas.png', dpi=300)
    plt.close(fig)


if __name__ == '__main__':
    expnames: list[str] = ['free', 'chlo', 'chlo-monthly',
                           'chlo-monthly-update', 'PC',
                           'PC-update', 'chlo-pc',]
    variables = ['OCN_PCO2', 'OXY', ]
    # for expname in expnames:
    #     config.exp = expname
    #     config.output_path = config.output_path_format.format(exp=expname)
    #     for varname in variables:
    #         get_series(varname)

    plot_timeseries()
