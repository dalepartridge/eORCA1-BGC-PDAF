import os
from typing import Union

import config
import utils
import operators

import numpy as np
import netCDF4


def read_observation(
    year: str,
    month: str,
    day: str,
    observation_type: str
) -> tuple[np.ndarray, np.ndarray]:
    """Read observation files.

    Args:
        year (str): Year of the observation.
        month (str): Month of the observation.
        day (str): Day of the observation.
        observation_type (str): Type of the observation. Can be either 'pc' or 'chlo'.

    Returns:
        Tuple[np.ndarray, np.ndarray]: Array containing the observation data and its uncertainty.
    """
    if not all([year, month, day, observation_type]):
        raise ValueError("Year, month, day, and observation_type cannot be empty")

    file_path: str = os.path.join(config.OBS_PATH, f"{observation_type}_{year}{month}{day}.nc")

    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Observation file '{file_path}' does not exist.")

    # this could be changed to adapt to more observation types
    varname = 'chlorophyll' if observation_type == 'chlo' else 'carbon'
    with netCDF4.Dataset(file_path, 'r') as dataset:
        if varname not in dataset.variables:
            raise ValueError(f"Variable 'chlorophyll' not found in file '{file_path}'")

        observation_array: np.ndarray = dataset[varname][:]
        observation_unc_array: np.ndarray = dataset[f'{varname}_unc'][:]

    return observation_array, observation_unc_array


def read_physical_observation(
    year: str,
    month: str,
    day: str,
    observation_type: str
) -> tuple[np.ndarray, np.ndarray]:
    """Read observations in the physical space.

    Args:
        year (str): Year of the observation.
        month (str): Month of the observation.
        day (str): Day of the observation.
        observation_type (str): Type of the observation. Can be either 'pc' or 'chlo'.

    Returns:
        Tuple[np.ndarray, np.ndarray]: Array containing the observation data and its uncertainty.
    """
    obs, obs_unc = read_observation(year, month, day, observation_type)
    obs, obs_unc = operators.unlog_observation(obs, obs_unc)
    return obs, obs_unc


def get_model_output_3d(
    year: str,
    cycle_month: str,
    day: int,
    level: int,
    ensemble_member_index: int,
    variable_info: dict
) -> np.ndarray:
    """Get model variable at given year/month/day at certain vertical level and ensemble member index.

    Args:
        year: Year of the model output.
        cycle_month: Cycling of the model output.
        day: Day of the model output.
        level: Vertical level of the model output.
        ensemble_member_index: Index of the ensemble member.
        variable_info: Dictionary containing the filename, operation, and variable names.

    Returns:
        numpy.ndarray: Array containing the model data.
    """
    if level < 0 or ensemble_member_index < 0:
        raise ValueError("Level and ensemble member index must be positive integers.")

    if {'filename', 'op', 'varname'} - variable_info.keys():
        raise ValueError(f"Variable info must contain 'filename', 'op', and 'varname' keys.")

    file_path = os.path.join(config.OUTPUT_PATH, year, cycle_month, f'ensemble_{str(ensemble_member_index)}', variable_info['filename'])
    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"Model file '{file_path}' does not exist.")


    with netCDF4.Dataset(file_path, 'r') as dataset:
        variable_names = [dataset[var_name] for var_name in variable_info['varname']]
        if any(var is None for var in variable_names):
            raise ValueError(f"Could not find variables '{variable_info['varname']}' in file '{file_path}'.")
        val = variable_info['op'](np.zeros((config.ny, config.nx)), level, day - 1, *variable_names)

    return val


def read_model_ensemble_mean_output_3d(
    year: str,
    cycle_month: str,
    day: int,
    level: int,
    variable_info: dict
) -> np.ndarray:
    """read ensemble mean model variable at given year/month/day at certain vertical level.

    Args:
        year: Year of the model output.
        cycle_month: Cycling of the model output.
        day: Day of the model output.
        level: Vertical level of the model output.
        variable_info: Dictionary containing the filename, operation, and variable names.

    Returns:
        numpy.ndarray: Array containing the model data.
    """
    if level < 0:
        raise ValueError("Level must be positive integers.")

    if {'filename', 'op', 'varname'} - variable_info.keys():
        raise ValueError(f"Variable info must contain 'filename', 'op', and 'varname' keys.")

    file_path = os.path.join(config.OUTPUT_PATH, year, cycle_month, variable_info['filename'])
    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"Model file '{file_path}' does not exist.")


    with netCDF4.Dataset(file_path, 'r') as dataset:
        variable_names = [dataset[var_name] for var_name in variable_info['varname']]
        if any(var is None for var in variable_names):
            raise ValueError(f"Could not find variables '{variable_info['varname']}' in file '{file_path}'.")
        val = variable_info['op'](np.zeros((config.ny, config.nx)), level, day - 1, *variable_names)

    return val


def read_model_std_output_3d(year: str, cycle_month: str, day: int, level: int, variable_info: dict) -> np.ndarray:
    """Read model variable at given year/month/day at certain vertical level.

    Args:
        year: Year of the model output.
        cycle_month: Month of the model output.
        day: Day of the model output.
        level: Vertical level of the model output.
        variable_info: Dictionary containing the filename, operation, and variable names.

    Returns:
        numpy.ndarray: Array containing the model data with shape (config.ny, config.nx).
    """
    it = int(day) - 1
    filename = os.path.join(year, cycle_month, variable_info['filename'])
    with netCDF4.Dataset(filename, 'r') as dataset:
        var_names = [dataset[var_name] for var_name in variable_info['varname']]
        if any(var is None for var in var_names):
            raise ValueError(f"Could not find variables '{variable_info['varname']}' in file '{filename}'.")
        array = variable_info['op'](np.zeros((config.ny, config.nx)), level, it, *var_names)
    return array


def read_pdaf_output_3d(year: str, cycle_month: str, month: str, day: str, ensemble_index: int, varname: str) -> np.ndarray:
    """Retrieve model variable at given year/month/day for a specific ensemble member.

    Args:
        year: Year of the PDAF output.
        cycle_month: Cycling month of the PDAF output.
        month: Month of the PDAF output.
        day: Day of the PDAF output.
        ensemble_index: Index of the ensemble member.

    Returns:
        numpy.ndarray: Array containing the model data.
    """
    filename = os.path.join(config.OUTPUT_PATH, year, cycle_month, 'analysis',
                            f'state_{varname}_{year}{month}{day}_{str(ensemble_index).zfill(3)}.nc')
    with netCDF4.Dataset(filename, 'r') as dataset:
        data = dataset[varname][:]
    return data


def get_std_pdaf_output_3d(year: str, cycle_month: str, month: str, day: str, varname: str) -> np.ndarray:
    """Calculate standard deviation of model variable at given year/month/day for all ensemble members.

    Args:
        year: Year of the PDAF output.
        cycle_month: Cycling month of the PDAF output.
        month: Month of the PDAF output.
        day: Day of the PDAF output.

    Returns:
        numpy.ndarray: Array containing the standard deviation of the model data for all ensemble members.
    """
    std_data = np.zeros((config.ne, 2, config.ny, config.nx), dtype=float)
    for i in range(config.ne):
        std_data[i] = read_pdaf_output_3d(year, cycle_month, month, day, i+1, varname)
    return std_data.std(axis=0, ddof=1)


def read_pdaf_std_output_3d(year: str, cycle_month: str, month: str, day: str, varname: str) -> np.ndarray:
    """Read the standard deviation of the model variable at given year, month, and day.

    Args:
        year: Year of the PDAF output.
        cycle_month: Cycling month of the PDAF output.
        month: Month of the PDAF output.
        day: Day of the PDAF output.

    Returns:
        numpy.ndarray: Array containing the standard deviation of the model data for all ensemble members.
    """
    filename = os.path.join(config.OUTPUT_PATH, year, cycle_month, 'analysis',
                            f'state_{varname}_{year}{month}{day}.nc')
    try:
        with netCDF4.Dataset(filename, 'r') as dataset:
            std_data = dataset[f'{varname}_std'][:]
    except OSError:
        std_data = np.empty((2, config.ny, config.nx))
        std_data.fill(np.nan)
    return std_data


def read_pdaf_ensemble_mean_output_3d(year: str, cycle_month: str, month: str, day: str, varname: str) -> np.ndarray:
    """Read the standard deviation of the model variable at given year, month, and day.

    Args:
        year: Year of the PDAF output.
        cycle_month: Cycling month of the PDAF output.
        month: Month of the PDAF output.
        day: Day of the PDAF output.

    Returns:
        numpy.ndarray: Array containing the standard deviation of the model data for all ensemble members.
    """
    filename = os.path.join(config.OUTPUT_PATH, year, cycle_month, 'analysis',
                            f'state_{varname}_{year}{month}{day}.nc')
    try:
        with netCDF4.Dataset(filename, 'r') as dataset:
            data = dataset[varname][:]
    except OSError:
        data = np.empty((2, config.ny, config.nx))
        data.fill(np.nan)
    return data


def get_pdaf_output_std(year: str, cycle_month: str, month: str, day: Union[str, int], varname: str) -> np.ndarray:
    """Retrieve the standard deviation of the model variable at given year, month, and day.

    Args:
        year: Year of the PDAF output.
        cycle_month: Cycling month of the PDAF output.
        month: Month of the PDAF output.
        day: Day of the PDAF output.

    Returns:
        numpy.ndarray: Array containing the standard deviation of the model data for all ensemble members.
    """
    day_str = str(day).zfill(2)
    if month in config.__monthsEns and year == '2016':
        return get_std_pdaf_output_3d(year, cycle_month, month, day_str, varname)
    else:
        return read_pdaf_std_output_3d(year, cycle_month, month, day_str, varname)


def calculate_innovation(observation: np.ndarray, model: np.ndarray) -> np.ndarray:
    """Calculate the difference between observation and model data.

    Parameters
    ----------
    observation : numpy.ndarray
        Array containing the observation data.
    model : numpy.ndarray
        Array containing the model data.

    Returns
    -------
    numpy.ndarray
        Array containing the difference between observation and model data.
    """
    mask = (observation > 1e14) | (model > 1e14)
    return np.ma.masked_where(mask, observation) - np.ma.masked_where(mask, model)


def get_ensemble_mean_from_model(year: str, cycle_month: str, month: str, day: int, level: int, variable_info: dict) -> np.ndarray:
    if month in config.__monthsEns and year == '2016':
        model_data = np.zeros((config.ny, config.nx), dtype=float)
        for ens_index in range(1, config.ne + 1):
            model_data += get_model_output_3d(year, cycle_month, day, level, ens_index, variable_info)/config.ne
    else:
        model_data = read_model_ensemble_mean_output_3d(year, cycle_month, day, level, variable_info)
    return model_data

