import multiprocessing as mp

import numpy as np

import diag
import config
import utils

def calculate_acc(year: str, cycle_month: str, month: str, day: int,
                  level: int, varname: str) -> np.ndarray:
    """Calculate anomaly correlation coefficient.

    Args:
        year (str): Year of the PDAF output.
        cycle_month (str): Cycling month of the PDAF output.
        month (str): Month of the PDAF output.
        day (int): Day of the PDAF output.
        level (int): Vertical level of the PDAF output.
        variable_info (dict): Dictionary containing the variable information.

    Returns:
        numpy.ndarray: Anomaly correlation coefficient.
    """
    f = np.load(f'clim_{varname}/clim_{year}{month}.npz')
    clim = np.ma.masked_array(f['model'], f['model_mask'])  # type: np.ndarray

    # Get the variable information
    variable_info = config.varinfo[varname]
    x_anom = diag.get_ensemble_mean_from_model(year, cycle_month, month, day, level, variable_info) - clim
    y, _ = diag.read_physical_observation(year, month, str(day).zfill(2), 'pc')
    y_anom = y - clim
    # Calculate the anomaly correlation coefficient
    acc = np.corrcoef(x_anom.ravel(), y_anom.ravel())
    return acc

def save_acc_multiprocessed() -> None:
    """Save monthly bias for multiple years and months using multiprocessing.

    Args:
        None

    Returns:
        None
    """
    varname = 'PC'
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
                # Calculate the number of days in the month
                days_in_month = utils.get_month_days(int(year), int(month))
                # Loop over the days in the month and calculate the monthly climatology
                rng = range(1, days_in_month + 1)
                for day in rng:
                    process: mp.Process = mp.Process(target=calculate_acc, args=(year, cycle_month, month, day, level, varname))  # type: mp.Process
                    processes.append(process)
                    process.start()

    for process in processes:
        process.join()


if __name__ == '__main__':
    save_acc_multiprocessed()