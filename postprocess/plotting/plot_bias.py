import multiprocessing as mp

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as mgs
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cmocean

import config
import utils
import diag


def save_monthly_bias(year: str, cycle_month: str, month: str, level: int, varname: str) -> None:
    """Save monthly bias for a given year, cycling month, and month.

    Args:
        year (str): The year of the PDAF output.
        cycle_month (str): The cycling month of the PDAF output.
        month (str): The month of the PDAF output.

    Returns:
        None
    """
    f = np.load(f'clim_{varname}/clim_{year}{month}.npz')
    monthly_model_data = np.ma.masked_array(f['model'], f['model_mask'])  # type: np.ndarray
    monthly_obs_data = np.ma.masked_array(f['obs'], f['obs_mask'])  # type: np.ndarray
    b: np.ndarray = diag.calculate_innovation(monthly_model_data, monthly_obs_data)  # type: np.ndarray
    np.savez(f'bias/bias_{year}{month}.npz', bias=b.data, mask=b.mask)  # Save monthly bias to a numpy archive


def save_monthly_bias_multiprocessed() -> None:
    """Save monthly bias for multiple years and months using multiprocessing.

    Args:
        None

    Returns:
        None
    """
    varname = 'chlorophyll'
    level = 0
    years: list[str] = ['2015', '2016']  # type: List[str]
    cycle_months: list[str] = [str(month).zfill(2) for month in range(1, 13, 2)]  # type: List[str]

    processes: list[mp.Process] = []  # type: List[mp.Process]
    for year in years:
        for cycle_month in cycle_months:
            if year == '2015' and cycle_month == '01':
                months: list[str] = [str(int(cycle_month) + 1).zfill(2), ]  # type: List[str]
            else:
                months: list[str] = [cycle_month, str(int(cycle_month) + 1).zfill(2), ]  # type: List[str]

            for month in months:
                process: mp.Process = mp.Process(target=save_monthly_bias, args=(year, cycle_month, month, level, varname))  # type: mp.Process
                processes.append(process)
                process.start()

    for process in processes:
        process.join()


def plot_monthly_bias(year:str, month: str) -> None:
    """Plot the monthly bias for multiple months in a single year.

    Args:
        None

    Returns:
        None
    """
    mask = utils.get_land_mask()
    lons, lats = utils.get_coord()

    fig = plt.figure()
    fig.clf()
    gs = mgs.GridSpec(1, 1, figure=fig, wspace=0.01, left=0., right=1., bottom=0.0, top=1.)
    ax = fig.add_subplot(gs[0], projection=ccrs.Robinson())
    f = np.load(f'bias/bias_{year}{month}.npz', allow_pickle=True)
    bias = f['bias']
    mask = f['mask']
    ax.set_global()
    pc = ax.pcolormesh(lons, lats, bias, cmap=cmocean.cm.balance,
                        norm=mcolors.CenteredNorm(vcenter=0.),
                        transform=ccrs.PlateCarree())
    ax.coastlines(color='k', linewidth=.8)
    ax.set_title(f'bias in {year}-{month}')
    fig.colorbar(pc, ax=ax, orientation='horizontal', shrink=0.6, pad=0.02)
    fig.savefig()
    plt.close(fig)


def plot_monthly_bias_multiprocessed() -> None:
    """Plot monthly bias for multiple years and months using multiprocessing.

    Args:
        None

    Returns:
        None
    """
    years: list[str] = ['2015', '2016']  # type: List[str]
    months: list[str] = [str(month).zfill(2) for month in range(1, 13)]  # type: List[str]

    processes: list[mp.Process] = []  # type: List[mp.Process]
    for year in years:
        for month in months:
            process: mp.Process = mp.Process(target=plot_monthly_bias, args=(year, month))  # type: mp.Process
            processes.append(process)
            process.start()

    for process in processes:
        process.join()


if __name__ == '__main__':
    plot_monthly_bias_multiprocessed()