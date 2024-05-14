import multiprocessing as mp
from typing import Union
import datetime

import numpy as np
import matplotlib.pyplot as plt

import utils
import diag


def save_spread(year: str, month: str, day: Union[str, int], cycle_month: str) -> None:
    """Save standard deviation of chlorophyll for a given year, month, day, and cycling month to a numpy archive.

    Args:
        year: str
            Year of the PDAF output.
        month: str
            Month of the PDAF output.
        day: Union[str, int]
            Day of the PDAF output.
        cycling_month: str
            Cycling month of the PDAF output.
    """
    spread_data = {
        'carbon': diag.get_pdaf_output_std(year, cycle_month, month, day, 'carbon')
    }
    file_name = f'std/std_{year}{month}{str(day).zfill(2)}.npz'
    np.savez(file_name, **spread_data)


def save_spread_multiprocessed() -> None:
    """Save standard deviation of chlorophyll for multiple years and months using multiprocessing.
    """
    years = ['2015', '2016']
    cycle_months = [str(month).zfill(2) for month in range(1, 13, 2)]

    for year in years:
        for cycle_month in cycle_months:
            if year == '2015' and cycle_month == '01':
                months = [str(int(cycle_month) + 1).zfill(2), ]
            else:
                months = [cycle_month, str(int(cycle_month) + 1).zfill(2), ]

            processes = []
            for month in months:
                days = range(1, utils.get_month_days(int(year), int(month)) + 1)
                for day in days:
                    process = mp.Process(target=save_spread, args=(year, month, day, cycle_month))
                    processes.append(process)
                    process.start()

            for process in processes:
                process.join()


def plotSpreadSpatialMean(years: list[str], months: list[str]) -> None:
    """Plot the spatial mean of chlorophyll for multiple months in a single year.

    Args:
        years (List[str]): The years of the PDAF output.
        months (List[str]): The months of the PDAF output.

    Returns:
        None
    """
    std_f: list[np.ndarray] = []  # Standard forecast spread
    std_a: list[np.ndarray] = []  # Standard analysis spread
    t: list[datetime.date] = []  # Timestamps
    mask = utils.get_land_mask()

    for year in years:
        for month in months:
            if year == '2015' and month == '01':
                rng: range = range(31, utils.get_month_days(int(year), int(month)))
            else:
                rng: range = range(utils.get_month_days(int(year), int(month)))
            for it in rng:
                t.append(datetime.date(int(year), int(month), it + 1))
                x = np.load(f'std/std_{year}{month}{str(it + 1).zfill(2)}.npz',
                            allow_pickle=True)['chlorophyll'][0].astype(np.float64)
                invalid = np.logical_or(mask, x > 1e16)
                std = np.nanmean(x[~invalid])
                print(std)
                std_f.append(std)
                x = np.load(f'std/std_{year}{month}{str(it + 1).zfill(2)}.npz',
                            allow_pickle=True)['chlorophyll'][1].astype(np.float64)
                invalid = np.logical_or(mask, x > 1e16)
                std = np.nanmean(x[~invalid])
                std_a.append(std)

    fig: plt.Figure = plt.figure(1)
    fig.clf()
    ax: plt.Axes = fig.add_subplot(111)
    ax.plot(t, std_f, label='forecast spread')
    ax.plot(t, std_a, label='analysis spread')
    ax.legend()
    fig.savefig('spreadSeries.png', dpi=300)
    plt.close(fig)


if __name__ == '__main__':
    save_spread_multiprocessed()
