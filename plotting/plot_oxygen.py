import matplotlib.pyplot as plt
import matplotlib.gridspec as mgs
import matplotlib.dates as mdates
import matplotlib.ticker as mticker
import numpy as np

plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']
plt.rcParams['font.size'] = 16

# Load the data
variables = ['oxygen_all_lvl_monthly', 'TCO2_monthly', 'CO2FLUX_monthly',
                'OCN_PCO2_monthly', 'O2FLUX_monthly',
                ]
variable_titles = ['Ocean oxygen', 'Total carbon', 'Air-sea CO2 flux',
                    'Ocean pCO2', 'Air-sea O2 flux', ]
ylabels = ['Ocean oxygen', 'Total carbon', 'Air-sea CO2 flux',
            'Ocean pCO2', 'Air-sea O2 flux',
             ]
varname:str ='oxygen_all_lvl_monthly'
varname_title:str = 'Total oxygen'
ylabel:str = 'Total oxygen'
for varname, varname_title, ylabel in zip(variables, variable_titles, ylabels):
    data: dict[str, list[float]] = {}
    EXPs : list[str] = ['PC', 'free', 'chlo', 'chlo-monthly', 'chlo-pc']

    for EXP in ['PC', 'free', 'chlo', 'chlo-monthly', 'chlo-pc']:
        data[EXP] = []
        for year in ['2015', ]:
            for month in range(2, 13):
                f:np.lib.npyio.NpzFile  = np.load(f'{EXP}/ensmean_{varname}_{year}{str(month).zfill(2)}01.npz')
                data[EXP].append( np.ma.sum(np.ma.masked_array(f['model'], f['model_mask'])) )

    t = [np.datetime64(f'{year}-{str(month).zfill(2)}-01') for year in ['2015', ] for month in range(2, 13)]
    fig: plt.Figure = plt.figure()
    fig.clf()
    gs = mgs.GridSpec(1, 1, figure=fig, wspace=0.01, left=0.14, right=0.99, bottom=0.17, top=0.94)
    ax:plt.Axes = fig.add_subplot(gs[0])
    lines = []
    for EXP, linestyle, color in zip(EXPs,
                                    ['--', '-', '-.', ':', '-'],
                                    ['k', 'r', 'k', 'k', 'k']):
            line, = ax.plot(t, data[EXP], color=color, linestyle=linestyle)
            lines.append(line)

    locator = mdates.MonthLocator(bymonth=range(2, 13, 2))
    ax.xaxis.set_major_locator(locator)

    # Set the major formatter to display the date in 'Month-Day' format
    ax.xaxis.set_major_formatter(mdates.AutoDateFormatter(locator))
    ax.tick_params(axis='x', rotation=15)
    # fig.legend(loc='outside lower center',
    #         handles=lines,
    #         labels=['Carbon', 'Freerun', 'daily Chlo', 'monthly Chlo', 'monthly Chlo & Carbon'],
    #         ncols=5, fontsize=10)
    ax.set_title(varname_title)
    ax.set_xlabel('time')
    ax.set_ylabel(ylabel)
    fig.savefig(f'timeseries_{varname}.png', dpi=300)
    plt.close(fig)

    fig = plt.figure(figsize=(8, 0.4))
    fig.clf()
    fig.legend(loc='center',
        handles=lines,
        labels=['Carbon', 'Freerun', 'daily Chlo', 'monthly Chlo', 'monthly Chlo & Carbon'],
        ncols=5, fontsize=10)
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
    fig.savefig('legend.png', dpi=300)