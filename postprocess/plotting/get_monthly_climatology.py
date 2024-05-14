import os
import multiprocessing as mp

import numpy as np

import config
import utils
import diag

def get_monthly_climatology(year: str, cycle_month: str, month: str, level: int, varname: str) -> np.ndarray:
    """Calculate (monthly) climatology of variable in a file.

    Args:
        year: PDAF output year.
        cycle_month: PDAF output cycling month.
        month: PDAF output month.
        level: Vertical level.
        varname: Name of the variable.

    Returns:
        numpy.ndarray: Array containing the monthly climatology of the model data.
    """
    # Calculate the number of days in the month
    days_in_month = utils.get_month_days(int(year), int(month))

    # Initialize arrays to store monthly model and observation data
    monthly_model_data = np.zeros((config.ny, config.nx), dtype=float)
    monthly_obs_data = np.zeros((config.ny, config.nx), dtype=float)

    # Determine the observation type based on the variable name
    observation_type = 'chlo' if varname == 'chlorophyll' else 'pc'

    # Get the variable information
    variable_info = config.varinfo[varname]

    # Generate the filename parameters
    # numdays is only used to calculate the filename
    numdays = utils.get_month_days(int(year), int(cycle_month) + 1)
    if year == '2015' and cycle_month == '01':
        filerange = f'{year}{str(int(cycle_month)+1).zfill(2)}01_{year}{str(int(cycle_month)+1).zfill(2)}{numdays}'
    else:
        filerange = f'{year}{cycle_month}01_{year}{str(int(cycle_month)+1).zfill(2)}{numdays}'
    variable_info['filename'] = variable_info['file_format'].format(filerange=filerange, year=year, month=month)

    # Loop over the days in the month and calculate the monthly climatology
    rng = range(1, days_in_month + 1)
    for day in rng:
        day_str = str(day).zfill(2)
        # Calculate the ensemble mean for the model data
        monthly_model_data += diag.get_ensemble_mean_from_model(year, cycle_month, month, day, level, variable_info) / days_in_month
        # Read the physical observation
        obs, _ = diag.read_physical_observation(year, month, day_str, observation_type)
        # Calculate the monthly climatology for the observation data
        monthly_obs_data +=  obs / days_in_month

    return monthly_model_data, monthly_obs_data


def save_monthly_climatology(
    year: str,
    cycle_month: str,
    month: str,
    level: int,
    varname: str
) -> None:
    """Save monthly climatology for a given year, cycling month, month, and variable.

    Args:
        year (str): The year of the PDAF output.
        cycle_month (str): The cycling month of the PDAF output.
        month (str): The month of the PDAF output.
        level (int): The vertical level of the PDAF output.
        varname (str): The name of the variable.

    Returns:
        None
    """
    # Create the output directory if it doesn't exist
    os.makedirs(f'clim_{varname}', exist_ok=True)

    # Calculate the monthly climatology of the variable
    monthly_model_data, monthly_obs_data = get_monthly_climatology(year, cycle_month, month, level, varname)

    # Save the monthly climatology to a numpy archive
    np.savez(
        f'clim_{varname}/clim_{year}{month}.npz',
        model=monthly_model_data.data,  # type: np.ndarray
        model_mask=monthly_model_data.mask,  # type: np.ndarray
        obs=monthly_obs_data,  # type: np.ndarray
        obs_mask=monthly_obs_data.mask  # type: np.ndarray
    )


def save_monthly_climatology_multiprocessed() -> None:
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
                process: mp.Process = mp.Process(target=save_monthly_climatology, args=(year, cycle_month, month, level, varname))  # type: mp.Process
                processes.append(process)
                process.start()

    for process in processes:
        process.join()