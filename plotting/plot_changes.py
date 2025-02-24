"""plot yearly average"""
import itertools
import multiprocessing as mp
import os

import cartopy.crs as ccrs
import cartopy.feature as cfeature  # type: ignore
import cmocean  # type: ignore # pylint: disable=unused-import
import matplotlib.colors as mcolors
import matplotlib.gridspec as mgs
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


def save_obs_data(vname: str) -> None:
    """Sum up data seasonally from monthly observation data.

    Parameters
    ----------
    vname : str
        The variable name.
    """
    config.exp = ''
    years: list[str]
    years = ['2015', '2016']
    months: list[str]
    months = [str(month).zfill(2) for month in range(1, 13)]

    data: np.ndarray = np.zeros((config.ny, config.nx))
    valid_counts: np.ndarray = np.zeros((config.ny, config.nx),
                                        dtype=int)
    yearmonths = itertools.product(years, months)
    for year, month in yearmonths:
        if year == '2015' and month == '01':
            continue
        f_info = file_info.FileInfo(year, month, '01', '')
        obs_data, _ = diag.read_physical_observation(
            f_info.get_obsfilename(vname), vname)
        obs_mask = np.logical_not(obs_data.mask)
        data += np.where(obs_mask, obs_data, 0)
        valid_counts += obs_mask.astype(int)

    # Calculate the temporal average
    valid_mask = valid_counts > 0
    data = np.where(
        valid_mask, data / valid_counts,
        np.nan)

    # Save the seasonal data to a file or process it as needed
    if vname == 'pc':
        vname = 'nitrogen-monthly'
    np.savez(
        os.path.join('data', 'obs', f'clim_{vname}.npz'),
        clim=data)


def sum_data(
        clim, varname: str, n_months: int,
        f_info: file_info.FileInfo) -> None:
    """sum the data"""

    fname = f_info.get_nemo_filename()
    data = diag.get_model_output_ensemble_mean(fname, f_info.stats,
                                               1, np.arange(1),
                                               config.varinfo[varname])
    clim[:] = clim[:] + data.ravel()/n_months


def save_data_per_var(varname: str) -> None:
    """save yearly observation and model output data"""
    # get model output file format
    file_format: str = str(config.varinfo[varname]['file_format'])
    # initialise the model climatology
    clim_model = mp.Array('d', np.ma.zeros(config.ny*config.nx))
    # clim_obs: np.ma.MaskedArray = np.ma.zeros((config.ny, config.nx))
    # # get observation type
    # if varname == 'nitrogen-monthly':
    #     obs_type = 'pc'
    # elif varname == 'chlo-monthly':
    #     obs_type = 'chlo-monthly'

    years: list[str] = ['2015', '2016']
    months: list[str] = [str(month).zfill(2)
                         for month in range(1, 13)]
    n_months: int = len(months)*len(years) - 1

    processes: list[mp.Process] = []

    yearmonth = itertools.product(years, months)
    for year, month in yearmonth:
        if year == '2015' and month == '01':
            continue
        f_info: file_info.FileInfo = file_info.FileInfo(year, month,
                                                        '01', file_format)
        p = mp.Process(target=sum_data, args=(
            clim_model, varname, n_months, f_info))

        processes.append(p)

    processes_batch = itertools.batched(processes, config.n_process)
    for batch in processes_batch:
        for p in batch:
            p.start()
        for p in batch:
            p.join()
    # save the data
    np.savez(f'data/{config.exp}/clim_{varname}.npz',
             clim=np.asarray(clim_model))


def save_data() -> None:
    """save climatology data for all variables and experiments.
    """
    os.makedirs('data', exist_ok=True)
    expnames: list[str] = ['free', 'chlo', 'chlo-monthly',
                           'chlo-monthly-update', 'PC',
                           'PC-update', 'chlo-pc',
                           ]

    for expname in expnames:
        os.makedirs(f'data/{expname}', exist_ok=True)
        config.exp = expname
        config.output_path = config.output_path_format.format(exp=expname)
        for vname in ['nitrogen-monthly', 'chlo-monthly']:
            save_data_per_var(vname)


def plot_axes(fig, ax, expname: str, vname: str) -> None:
    """plot the axes"""
    lons, lats = utils.get_coord()
    ax.set_global()
    if expname in ['free', 'obs']:
        data = np.load(os.path.join('data', expname,
                       f'clim_{vname}.npz'))['clim']
        if expname == 'free':
            data = data.reshape(config.ny, config.nx)
        vmax = 0.7 if vname == 'chlo-monthly' else 0.6
        norm = mcolors.Normalize(vmin=0., vmax=vmax)
        cmap = 'cmo.algae'
    else:
        clim = np.load(os.path.join('data', expname,
                       f'clim_{vname}.npz'))['clim']
        clim_ref = np.load(os.path.join(
            'data', 'free', f'clim_{vname}.npz'))['clim']
        clim = clim.reshape(config.ny, config.nx)
        clim_ref = clim_ref.reshape(config.ny, config.nx)
        data = clim - clim_ref
        vmax = 0.05 if vname == 'chlo-monthly' else 0.1
        if expname == 'chlo-monthly' and vname == 'nitrogen-monthly':
            vmax = 0.01
        if vname == 'chlo-monthly':
            if expname == 'PC-update' or expname == 'chlo':
                vmax = 0.2
        norm = mcolors.TwoSlopeNorm(vcenter=0., vmin=-vmax, vmax=vmax)
        cmap = 'cmo.balance'
    print(np.nanpercentile(np.abs(data), 95))
    pc = ax.pcolormesh(
        lons, lats, data, transform=ccrs.PlateCarree(),
        cmap=cmap, norm=norm)
    ax.set_title(config.exp_labels[expname])
    ax.coastlines(color='k', linewidth=.8)
    ax.add_feature(cfeature.LAND, zorder=3)
    # locator = mticker.FixedLocator(
    #     [-0.3, -0.01, -0.001, 0, 0.001, 0.01, 0.3])
    # formatter = mticker.FixedFormatter(
    #     [str(i).zfill(2)
    #      if i != 0 else '0'
    #      for i in [-0.3, -0.01, -0.001, 0, 0.001, 0.01, 0.3]])
    fig.colorbar(pc, ax=ax, orientation='horizontal',
                 pad=0.05, shrink=0.5)  # , ticks=locator, format=formatter)


def plot() -> None:
    """plot the differences between DA experiment and free run"""
    os.makedirs('figs', exist_ok=True)
    expnames: list[str] = ['obs', 'chlo', 'chlo-monthly',
                           'chlo-monthly-update', 'free', 'PC',
                           'PC-update', 'chlo-pc',
                           ]
    vnames: list[str] = ['nitrogen-monthly', 'chlo-monthly']

    for vname in vnames:
        fig: plt.Figure = plt.figure()
        w, h = fig.get_size_inches()
        fig.set_size_inches(w*2.5, h*1.2)
        gs: mgs.GridSpec = mgs.GridSpec(2, 4,
                                        figure=fig,
                                        wspace=0.01,
                                        hspace=0.17,
                                        left=0., right=1.,
                                        bottom=0.01, top=0.93)
        for i, expname in enumerate(expnames):
            ax = fig.add_subplot(
                gs[i], projection=ccrs.Robinson())
            plot_axes(fig, ax, expname, vname)
        fig.savefig(f'figs/clim_diff_{vname}.png', dpi=300)


if __name__ == '__main__':
    # save_data()
    save_obs_data('pc')
    save_obs_data('chlo-monthly')
    plot()
