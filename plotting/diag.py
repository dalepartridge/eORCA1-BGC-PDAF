"""Data reading and diagnostic functions.
"""
import os
import typing

import netCDF4  # type: ignore
import numpy as np

import config
import operators


# functions dealing with observation files
def read_observation(filename: str, obs_type: str
                     ) -> tuple[np.ma.MaskedArray,
                                np.ma.MaskedArray]:
    """Read observation files.

    Here the observation files stores the parameters of
    base-10 log-normal distribution (mean and std. of a Gaussian).

    Parameters
    ----------
    filename : str
        observation filename.
    obs_type : str
        Type of the observation.
                It can be `pc`, `chlo`, or `chlo-monthly`.

    Returns
    -------
    Tuple[np.ma.MaskedArray, np.ma.MaskedArray]
        Array containing the observation data and its uncertainty.
    """

    if not os.path.exists(filename):
        raise FileNotFoundError(
            f"Observation file '{filename}' does not exist."
        )

    varnames: dict[str, str] = {'pc': 'carbon',
                                'chlo': 'chlorophyll',
                                'chlo-monthly': 'chlorophyll'}

    # this could be changed to adapt to more observation types
    with netCDF4.Dataset(filename, 'r') as dataset:  # pylint: disable=no-member
        observation_array: np.ma.MaskedArray = dataset[
            varnames[obs_type]
        ][:]
        observation_unc_array: np.ma.MaskedArray = dataset[
            f'{varnames[obs_type]}_unc'
        ][:]

    return observation_array, observation_unc_array


def read_physical_observation(filename: str, obs_type: str
                              ) -> tuple[np.ma.MaskedArray,
                                         np.ma.MaskedArray]:
    """Read observation files.

    Here the parameters of the base-10 log-normal distribution are
    transformed back to the mean of the random variable.

    Parameters
    ----------
    filename : str
        observation filename.
    obs_type : str
        Type of the observation.
        It can be `pc`, `chlo`, or `chlo-monthly`.

    Returns
    -------
    Tuple[np.ma.MaskedArray, np.ma.MaskedArray]
        Array containing the observation data and its uncertainty.
    """
    obs, obs_unc = read_observation(filename, obs_type)
    obs, obs_unc = operators.get_lognormal_mean_std(obs, obs_unc)
    return obs, obs_unc


def read_physical_observation_mode(
        filename: str, obs_type: str) -> tuple[
        np.ma.MaskedArray, np.ma.MaskedArray]:
    """Read observation files.

    Here the parameters of the base-10 log-normal distribution are
    transformed to the mode of the random variable.

    Parameters
    ----------
    filename : str
        observation filename
    obs_type : str
        Type of the observation.
        It can be `pc`, `chlo`, or `chlo-monthly`.

    Returns
    -------
    Tuple[np.ma.MaskedArray, np.ma.MaskedArray]
        Array containing the observation data and its uncertainty.
    """
    obs, obs_unc = read_observation(filename, obs_type)
    obs_mode = operators.get_lognormal_mode(obs, obs_unc)
    _, obs_unc = operators.get_lognormal_mean_std(obs, obs_unc)
    return obs_mode, obs_unc

# functions dealing with model output


def read_model_output_3d(filename: str, day: int, level: np.ndarray,
                         variable_info: dict) -> np.ma.MaskedArray:
    """Reading variable from NEMO output.

    The variable will be transformed by
    given operators provided in `variable_info`.

    Parameters
    ----------
    filename: str
        Filename of the model output.
    day : int
        Day in the month of the model output.
    level : np.ndarray
        Vertical level of the model output.
    variable_info : dict
        Dictionary containing the operator and variable names.

    Returns
    -------
    np.ma.MaskedArray
        Array containing the model data.
    """
    if {'op', 'varname'} - variable_info.keys():
        raise ValueError(
            "Variable info must contain 'op', and 'varname' keys."
            "Currently, it has {variable_info.keys()}"
        )

    with netCDF4.Dataset(filename, 'r') as dataset:  # pylint: disable=no-member
        var_handles = [dataset[var_name]
                       for var_name in variable_info['varname']]
        val = variable_info['op'](np.ma.zeros(
            (len(level), config.ny, config.nx)),
            level, day - 1, *var_handles
        )
    return val


def read_model_output_3d_stats(filename: str, day: int,
                               level: np.ndarray,
                               variable_info: dict
                               ) -> np.ma.MaskedArray:
    """Reading ensemble mean or standard deviation of
    the model variable (processed .nc files)
    at given day and vertical levels.

    Parameters
    ----------
    filename: str
        Filename of the model output.
    day : int
        Day in the month of the model output.
    level : np.ndarray
        Vertical level of the model output.
    variable_info : dict
        Dictionary containing the operators and variable names.

    Returns
    -------
    np.ma.MaskedArray
        Array containing the model data.
    """

    if {'op', 'varname'} - variable_info.keys():
        raise ValueError(
            "Variable info must contain 'op', and 'varname' keys.")

    with netCDF4.Dataset(filename, 'r') as dataset:  # pylint: disable=no-member
        var_handles = [dataset[var_name]
                       for var_name in variable_info['varname']]
        val: np.ma.MaskedArray = variable_info['op'](np.ma.zeros(
            (len(level), config.ny, config.nx)),
            level, day - 1, *var_handles
        )
    return val


def calculate_model_std_output_3d(filename: str,
                                  day: int, level: np.ndarray,
                                  variable_info: dict
                                  ) -> np.ma.MaskedArray:
    """Calculate standard deviation of model variable at
    given year/month/day from all ensemble members.

    Parameters
    ----------
    filename: str
        Filename of the model output.
    day : int
        Day in the month of the model output.
    level : np.ndarray
        Vertical level of the model output.
    variable_info : dict
        Dictionary containing the operators and variable names.

    Returns
    -------
    np.ma.MaskedArray
        Array containing the standard deviation of given variable.
    """
    std_data: np.ma.MaskedArray = np.ma.zeros(
        (config.ne, len(level), config.ny, config.nx), dtype=float)
    for i in range(config.ne):
        std_data[i] = read_model_output_3d(
            filename.format(i + 1), day, level, variable_info)
    return np.ma.std(std_data, axis=0, ddof=1)


def get_model_output_std(
        filename: str, read_stats: bool, day: int, level: np.ndarray,
        variable_info: dict) -> np.ma.MaskedArray:
    """Retrieve the standard deviation of
    the model variable at given year, month, and day.

    This function either reads directly from a file using
    :func:`read_model_output_3d_stats`,
    or calculate the standard deviation using
    :func:`calculate_model_std_output_3d`.

    Parameters
    ----------
    filename: str
        Filename of the model output.
    read_stats: bool
        If True, read the standard deviation from the PDAF output.
        Otherwise, calculate the standard deviation.
    day : int
        Day in the month of the model output.
    level : np.ndarray
        Vertical level of the model output.
    variable_info : dict
        Dictionary containing the operators and variable names.

    Returns
    -------
    np.ma.MaskedArray
        Array containing the model data.
    """

    if not read_stats:
        return calculate_model_std_output_3d(
            filename, day, level, variable_info)
    else:
        var_info: dict[str, typing.Union[str,
                                         typing.Union[
                                             list[str],
                                             typing.Callable
                                         ]
                                         ]
                       ] = dict()
        var_info['varname'] = variable_info['varname']
        var_info['filename'] = variable_info['filename']
        var_info['op'] = variable_info['op']
        if variable_info['op'] == operators.transform_sum:  # pylint: disable=comparison-with-callable
            # getting the std of summed random variable
            var_info['op'] = operators.identity

            if 'PHD' in variable_info['varname']:
                var_info['varname'] = ['nitrogen_std',]
            else:
                var_info['varname'] = ['chlorophyll_std',]

            return read_model_output_3d_stats(
                filename, day, level, var_info) * np.sqrt(
                config.ne / (config.ne - 1))
        elif variable_info['op'] == operators.identity:  # pylint: disable=comparison-with-callable
            # getting the std of the random variable
            return read_model_output_3d_stats(
                filename, day, level, var_info) * np.sqrt(
                config.ne / (config.ne - 1))
        else:
            raise ValueError(
                'unsupported operator in variable_info["op"]'
            )


def get_model_output_ensemble_mean(
        filename: str, read_stats: bool, day: int, level: np.ndarray,
        variable_info: dict) -> np.ma.MaskedArray:
    """Retrieve the ensemble mean of the model variable at
    given year, month, and day.

    This function either reads directly from a file using
    :func:`read_model_output_3d_stats`,
    or calculate the ensemble mean directly.

    Parameters
    ----------
    filename: str
        Filename of the model output.
    read_stats: bool
        If True, read the standard deviation from the PDAF output.
        Otherwise, calculate the standard deviation.
    day : int
        Day in the month of the model output.
    level : np.ndarray
        Vertical level of the model output.
    variable_info : dict
        Dictionary containing the filename, operation, and variable names.

    Returns
    -------
    np.ma.MaskedArray
        Array containing the model data.
    """
    model_data: np.ma.MaskedArray
    if not read_stats:
        model_data = np.ma.zeros(
            (len(level), config.ny, config.nx), dtype=float)
        for i in range(1, config.ne + 1):
            model_data = model_data + read_model_output_3d(
                filename.format(i), day, level, variable_info)/config.ne
    else:
        var_info: dict[str, typing.Union[str,
                                         typing.Union[list[str], typing.Callable]]] = dict()
        var_info['varname'] = variable_info['varname']
        var_info['op'] = variable_info['op']
        if variable_info['op'] == operators.transform_sum:  # pylint: disable=comparison-with-callable
            var_info['varname'] = [
                'nitrogen',] if 'PHD' in variable_info['varname'] else ['chlorophyll',]
            var_info['op'] = operators.identity
            model_data = read_model_output_3d_stats(
                filename, day, level, var_info)
        elif variable_info['op'] == operators.identity:  # pylint: disable=comparison-with-callable
            model_data = read_model_output_3d_stats(
                filename, day, level, var_info)
        else:
            raise ValueError('unsupported operator in variable_info["op"]')

    return model_data

# functions dealing with PDAF output


def read_pdaf_output_3d(filename: str, varname: str) -> np.ma.MaskedArray:
    """Read PDAF output file.

    Parameters
    ----------
    filename: str
        Filename of the PDAF output.
    varname: str
        Name of the variable to read.

    Returns
    -------
    np.ma.MaskedArray
        Array containing the PDAF output data.
    """
    with netCDF4.Dataset(filename, 'r') as dataset:  # pylint: disable=no-member
        data = dataset[varname][:]
    return data


def read_pdaf_std_output_3d(filename: str, varname: str) -> np.ma.MaskedArray:
    """Read the standard deviation of the PDAF output

    Parameters
    ----------
    filename: str
        Filename of the PDAF output.
    varname: str
        Name of the variable to read.

    Returns
    -------
    np.ma.MaskedArray
        Array containing the PDAF output of standard deviation
    """
    with netCDF4.Dataset(filename, 'r') as dataset:  # pylint: disable=no-member
        std_data = dataset[f'{
            varname}_std'][:]*np.sqrt(config.ne/(config.ne-1))
    return std_data


def calculate_std_pdaf_output_3d(
        filename: str, varname: str) -> np.ma.MaskedArray:
    """Calculate standard deviation of PDAF output.

    Parameters
    ----------
    filename: str
        Filename of the PDAF output.
    varname: str
        Name of the variable to read.

    Returns
    -------
    np.ma.MaskedArray
        Array containing the PDAF output of standard deviation.
    """
    std_data = np.ma.zeros((config.ne, 2, config.ny, config.nx), dtype=float)
    for i in range(config.ne):
        std_data[i] = read_pdaf_output_3d(
            filename.format(str(i+1).zfill(3)), varname)
    return std_data.std(axis=0, ddof=1)


def calculate_physical_std_pdaf_output_3d(
        filename: str, varname: str) -> np.ma.MaskedArray:
    """Calculate standard deviation of PDAF output in physical space.

    Parameters
    ----------
    filename: str
        Filename of the PDAF output.
    varname: str
        Name of the variable to read.

    Returns
    -------
    np.ma.MaskedArray
        Array containing the PDAF output of standard deviation.
    """
    std_data = np.ma.zeros((config.ne, 2, config.ny, config.nx), dtype=float)
    for i in range(config.ne):
        std_data[i] = read_pdaf_output_3d(
            filename.format(str(i+1).zfill(3)), varname)
        std_data[i] = 10**std_data[i]
    return std_data.std(axis=0, ddof=1)


def get_pdaf_output_std(
        filename: str, varname: str, read_stats: bool) -> np.ma.MaskedArray:
    """Retrieve the standard deviation of the PDAF output

    Parameters
    ----------
    filename: str
        Filename of the PDAF output.
    varname: str
        Name of the variable to read.
    read_stats: bool
        If True, read the standard deviation from the PDAF output directly.
        Otherwise, calculate the standard deviation.

    Returns
    -------
    np.ma.MaskedArray
        Array containing the PDAF output of standard deviation.
    """
    if not read_stats:
        return calculate_std_pdaf_output_3d(filename, varname)
    else:
        return read_pdaf_std_output_3d(filename, varname)


def get_pdaf_output_physical_std(
        filename: str, varname: str, read_stats: bool) -> np.ma.MaskedArray:
    """Retrieve the standard deviation of PDAF output in physical space.

    Parameters
    ----------
    filename: str
        Filename of the PDAF output.
    varname: str
        Name of the variable to read.
    read_stats: bool
        If True, read the standard deviation from the PDAF output directly.
        Otherwise, calculate the standard deviation.

    Returns
    -------
    np.ma.MaskedArray
        Array containing the PDAF output of standard deviation.
    """
    if not read_stats:
        return calculate_physical_std_pdaf_output_3d(filename, varname)
    else:
        model_data = read_pdaf_output_3d(filename, varname)
        std = read_pdaf_std_output_3d(filename, varname)
        _, std = operators.get_lognormal_mean_std(model_data, std)
        return std


def get_ensemble_mean_from_pdaf(
        filename: str, varname: str, read_stats: bool) -> np.ma.MaskedArray:
    """Retrieve the ensemble mean of PDAF output.

    Parameters
    ----------
    filename: str
        Filename of the PDAF output.
    varname: str
        Name of the variable to read.
    read_stats: bool
        If True, read the standard deviation from the PDAF output directly.
        Otherwise, calculate the standard deviation.

    Returns
    -------
    np.ma.MaskedArray
        Array containing the ensemble mean of PDAF output.
    """
    if not read_stats:
        model_data = np.ma.zeros((2, config.ny, config.nx), dtype=float)
        for i in range(1, config.ne + 1):
            model_data = model_data + read_pdaf_output_3d(
                filename.format(str(i).zfill(3)), varname)/config.ne
    else:
        model_data = read_pdaf_output_3d(filename, varname)
    return model_data


def get_physical_ensemble_mean_from_pdaf(
        filename: str, varname: str, read_stats: bool) -> np.ma.MaskedArray:
    """Retrieve the ensemble mean of PDAF output in physical space.

    Parameters
    ----------
    filename: str
        Filename of the PDAF output.
    varname: str
        Name of the variable to read.
    read_stats: bool
        If True, read the standard deviation from the PDAF output directly.
        Otherwise, calculate the standard deviation.

    Returns
    -------
    np.ma.MaskedArray
        Array containing the ensemble mean of PDAF output in physical space.
    """
    if not read_stats:
        model_data = np.ma.zeros((config.ny, config.nx), dtype=float)
        for i in range(1, config.ne + 1):
            model_data = model_data + 10**read_pdaf_output_3d(
                filename.format(str(i).zfill(3)), varname)/config.ne
    else:
        model_data = read_pdaf_output_3d(filename, varname)
        std = read_pdaf_std_output_3d(filename, varname)
        model_data, _ = operators.get_lognormal_mean_std(model_data, std)
    return model_data


def get_physical_ensemble_mode_from_pdaf(
        filename: str, varname: str, read_stats: bool) -> np.ma.MaskedArray:
    """Retrieve the ensemble mode of PDAF output in physical space.

    Parameters
    ----------
    filename: str
        Filename of the PDAF output.
    varname: str
        Name of the variable to read.
    read_stats: bool
        If True, read the standard deviation from the PDAF output directly.
        Otherwise, calculate the standard deviation.

    Returns
    -------
    np.ma.MaskedArray
        Array containing the ensemble mode of PDAF output in physical space.
    """
    model_data = get_ensemble_mean_from_pdaf(filename, varname, read_stats)
    model_unc = get_pdaf_output_std(filename, varname, read_stats)
    model_mode = operators.get_lognormal_mode(model_data, model_unc)
    return model_mode


def get_ensemble_anomaly_from_pdaf(
        filename: str, varname: str, read_stats: bool) -> np.ma.MaskedArray:
    """Calculate the anomaly of the pdaf ensemble at assimilation step

    Parameters
    ----------
    filename: str
        Filename of the PDAF output.
    varname: str
        Name of the variable to read.
    read_stats: bool
        If True, read the standard deviation from the PDAF output directly.
        Otherwise, calculate the standard deviation.

    Returns
    -------
    np.ma.MaskedArray
        Array containing the anomaly of the model data for all ensemble members.
    """
    mean_val: np.ma.MaskedArray = get_ensemble_mean_from_pdaf(
        filename, varname, read_stats)
    anomaly_mat: np.ma.MaskedArray = np.ma.zeros(
        (config.ne, 2, config.ny, config.nx), dtype=float)
    for i in range(config.ne):
        anomaly_mat[i] = read_pdaf_output_3d(
            filename.format(str(i + 1).zfill(3)),
            varname) - mean_val
    return anomaly_mat

# diagnostics functions


def calculate_innovation(observation: np.ma.MaskedArray,
                         model: np.ma.MaskedArray) -> np.ma.MaskedArray:
    """Calculate the difference between observation and model data.

    Parameters
    ----------
    observation : np.ma.MaskedArray
        Array containing the observation data.
    model : np.ma.MaskedArray
        Array containing the model data.

    Returns
    -------
    np.ma.MaskedArray
        Array containing the difference between observation and model data.
    """
    return observation - model


def calc_r_pp(chd: np.ma.MaskedArray, chn: np.ma.MaskedArray,
              phd: np.ma.MaskedArray, phn: np.ma.MaskedArray,
              prd: np.ma.MaskedArray, prn: np.ma.MaskedArray,
              med_xpar: np.ma.MaskedArray) -> tuple[np.ma.MaskedArray, np.ma.
                                                    MaskedArray]:
    """Calculate an estimate of the ratio of the
    chlorophyll primary production to the phytoplankton nitrogen.

    Parameters
    ----------
    CHD : np.ndarray
        Chlorophyll of the diatoms.
    CHN : np.ndarray
        Chlorophyll of the non-diatoms.
    PHD : np.ndarray
        Nitrogen of the diatoms.
    PHN : np.ndarray
        Nitrogen of the non-diatoms.
    PRD : np.ndarray
        Primary production of the diatoms.
    PRN : np.ndarray
        Primary production of the non-diatoms.
    MED_XPAR : np.ndarray
        Median of the photosynthetically active radiation.

    Returns
    -------
    Rn : np.ndarray
        Ratio of the chlorophyll primary production
        to the phytoplankton nitrogen of the non-diatoms.
    Rd : np.ndarray
        Ratio of the chlorophyll primary production
        to the phytoplankton nitrogen of the diatoms.
    """
    xthetam_max: float = 0.05
    xxi: float = 0.01257
    xaln: float = 15.0
    xthetam: np.ma.MaskedArray = chn*xxi/phn
    alpha_pn: np.ma.MaskedArray = xaln*xthetam
    rn: np.ma.MaskedArray = xthetam_max*prn/xthetam/alpha_pn/med_xpar
    xald: float = 11.25
    xthetad_max: float = 0.05
    xthetad = chd*xxi/phd
    alpha_pd: np.ma.MaskedArray = xald*xthetad
    rd: np.ma.MaskedArray = xthetad_max*prd/xthetad/alpha_pd/med_xpar
    return rn, rd
