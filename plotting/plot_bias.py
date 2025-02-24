"""Plot the observation - model bias
"""
import itertools
import os

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as mgs
import matplotlib.lines as mlines

import config
import diag
import file_info
import utils

plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']
plt.rcParams['font.size'] = 16


def get_bias(filename: str, obs_type: str,
             year: str, month: str, day: str) -> np.ma.MaskedArray:
    """Get the observation - model at given month day.

    Parameters
    ----------
    filename: str
        The climatology model filename in .npz format.
    obs_type: str
        observation type; it can be `pc`, `chlo`, or `chlo-monthly`.
    year: str
        The year.
    month: str
        The month.
    day: str
        The day.

    Returns
    -------
    np.ma.MaskedArray
        The observation - model.
    """
    f: np.lib.npyio.NpzFile = np.load(
        filename.format(
            year=year, month=month, day=day))

    model: np.ma.MaskedArray = np.ma.masked_array(
        f['model'], f['model_mask'])
    f_info = file_info.FileInfo(year, month, day, '')

    obs, _ = diag.read_observation(
        f_info.get_obsfilename(obs_type), obs_type)

    return obs - model


def save_monthly_bias(
        varname: str, obs_type: str, filename_model: str) -> None:
    """Save monthly bias for multiple years and months.

    Parameters
    ----------
    varname: str
        The variable name used for output filename as wel.
    obs_type: str
        observation type; it can be `pc`, `chlo`, or `chlo-monthly`.
    filename_model: str
        The climatology model filename in .npz format.
    """
    years: list[str] = ['2015', '2016']
    months: list[str] = [str(month).zfill(2)
                         for month in range(1, 13)]

    fname_output: str = os.path.join(
        'data', config.exp,
        'bias_{year}{month}{day}'f'_{varname}_pdaf.npz'
    )

    yearmonth = itertools.product(years, months)
    for year, month in yearmonth:
        if year == '2015' and month == '01':
            continue

        # Get the number of days in the month
        ndays = utils.get_month_days(int(year), int(month))
        days: range = range(ndays//2, ndays//2 + 1)
        if config.exp == 'chlo':
            days = range(1, ndays + 1)

        # Loop over days
        for day in days:
            day_str: str = str(day).zfill(2)

            b: np.ma.MaskedArray = get_bias(filename_model, obs_type,
                                            year, month, day_str)

            np.savez(
                fname_output.format(year=year, month=month, day=day_str
                                    ),
                bias=b.data,
                mask=b.mask
            )


def read_bias(varname: str) -> dict[str, np.ndarray]:
    """Read the bias for multiple years and months.

    Parameters
    ----------
    varname: str
        The variable name. It is either `nitrogen` or `chlorophyll`.

    Returns
    -------
    dict[str, np.ndarray]
        The bias for multiple experiments.
    """

    years: list[str] = ['2015', '2016']
    months: list[str] = [str(month).zfill(
        2) for month in range(1, 13)]

    exps: list[str]
    if varname == 'nitrogen':
        exps = ['PC', 'PC-update', 'chlo-pc',
                ]
    else:
        exps = ['chlo', 'chlo-monthly',
                'chlo-monthly-update', 'chlo-pc'
                ]

    fname_bias: str = os.path.join(
        'data',
        '{exp}/bias_{year}{month}{day}_{varname}_pdaf.npz'
    )

    exp_yearmonths = itertools.product(exps, years, months)

    bias: dict[str, np.ndarray] = {exp: np.array([]) for exp in exps}

    # get an array of o - b
    for exp, year, month in exp_yearmonths:
        if year == '2015' and month == '01':
            continue

        # Get the number of days in the month
        ndays: int = utils.get_month_days(int(year), int(month))
        days: range = range(ndays//2, ndays//2 + 1)
        if exp == 'chlo':
            days = range(1, ndays + 1)

        # Loop over days
        for day in days:
            f = np.load(fname_bias.format(
                exp=exp,
                year=year, month=month, day=str(day).zfill(2),
                varname=varname)
            )
            bias[exp] = np.append(
                bias[exp],
                np.ma.array(
                    f['bias'][0],
                    mask=f['mask'][0]).filled(np.nan).ravel()
            )
    return bias


def plot_axes(
        i: int, ax: plt.Axes, varname: str, exps: list[str],
        linestyles: list[str],
        colours: list[str]) -> plt.Axes:
    """Plot the histogram of o - b in PDAF output for single variable.

    Parameters
    ----------
    i: int
        The index of the variable.
    ax: plt.Axes
        The histogram axes object.
    varname: str
        The variable name.
    exps: list[str]
        The experiment names.
    linestyles: list[str]
        The linestyles for each experiment.
    colours: list[str]
        The colours for each experiment.

    Returns
    -------
    plt.Axes
        The histogram axes object.
    """
    bias = read_bias(varname)
    j = 0
    for exp, linestyle, colour in zip(exps,
                                      linestyles, colours):
        print(varname, exp, bias[exp].shape)
        if varname == 'nitrogen' and exp in ['chlo',
                                             'chlo-monthly',
                                             'chlo-monthly-update'
                                             ]:
            continue

        if varname == 'chlorophyll' and exp in ['PC', 'PC-update']:
            continue
        ax.hist(bias[exp], bins=100,
                density=True,
                histtype='step', color=colour,
                linestyle=linestyle,
                linewidth=3, label=exp)
        loc_text = 0.02 if i == 1 else 0.64
        ax.text(loc_text, 0.9 - 0.1*j,
                f'{config.exp_labels[exp]}:'
                f' {np.nanmean(bias[exp]):.2f}',
                transform=ax.transAxes, fontsize=12)
        j += 1
    return ax


def plot_bias_histogram(varnames: list[str]) -> None:
    """Plot the histogram of o - b in PDAF output.

    Here, the o - b is the mu parameter in the log-normal
    distribution.

    Parameters
    ----------
    varname: str
        The variable name.
    """
    assert len(varnames) == 2, 'Only two variables are plotted.'

    exps = ['chlo', 'chlo-monthly',
            'chlo-monthly-update', 'PC',
            'PC-update', 'chlo-pc',
            ]
    linestyles = [':', ':', ':', '-', ':', ':']
    colours = ['#1E88E5', '#FFC107', '#48B03E',
               'k', '#2CF8BA', '#D81B1B']

    # create the figure
    fig = plt.figure()
    fig.clf()
    w, h = fig.get_size_inches()
    fig.set_size_inches(w*2, h)

    gs = mgs.GridSpec(1, 2, figure=fig, wspace=0.17,
                      left=0.08, right=0.99, bottom=0.25, top=0.94)

    for i, varname in enumerate(varnames):
        ax = fig.add_subplot(gs[i])
        # plot the histogram for current varname
        ax = plot_axes(i, ax, varname, exps, linestyles, colours)
        unit = r'$mmol/m^3$' if varname == 'nitrogen' else r'$mg/m^3$'
        ax.set_xlabel(f'o - b ({unit})')
        ax.set_title(varname)
        ax.set_ylabel('Frequency')

    lines = [mlines.Line2D([], [], color=color, linestyle=linestyle)
             for linestyle, color in zip(linestyles, colours)
             ]

    fig.legend(
        loc='outside lower center', handles=lines,
        labels=[config.exp_labels[exp] for exp in exps],
        ncols=6, fontsize=14)

    fig.savefig('bias_hist.pdf', dpi=300)
    plt.close(fig)


if __name__ == '__main__':
    expnames: list[str] = ['chlo', 'chlo-monthly',
                           'chlo-monthly-update', 'PC',
                           'PC-update', 'chlo-pc',
                           ]
    for expname in expnames:
        config.exp = expname
        if expname in ['PC', 'PC-update', 'chlo-pc', ]:
            save_monthly_bias(
                'nitrogen', 'pc', os.path.join(
                    'data', expname,
                    'pdaf_ensmean_nitrogen_{year}{month}{day}.npz'))
        if expname in ['chlo-monthly', 'chlo-pc', 'chlo-monthly-update']:
            save_monthly_bias(
                'chlorophyll', 'chlo-monthly', os.path.join(
                    'data', expname,
                    'pdaf_ensmean_chlorophyll_{year}{month}{day}.npz'))
        if expname in ['chlo', ]:
            save_monthly_bias(
                'chlorophyll', 'chlo', os.path.join(
                    'data', expname,
                    'pdaf_ensmean_chlorophyll_{year}{month}{day}.npz'))

    plot_bias_histogram(['chlorophyll', 'nitrogen',])

    # plot_bias_timeseries(['chlorophyll', 'nitrogen',])
