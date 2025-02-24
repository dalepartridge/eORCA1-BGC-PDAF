"""transformations used in diagnostics"""
import numpy as np


def transform_sum_log10(val: np.ma.MaskedArray, lev: np.ndarray,
                        it: int, *variables
                        ) -> np.ma.MaskedArray:
    """Sum over variables and
    apply log10 function to the summed variable

    Parameters
    ----------
    val : np.ndarray, shape (M, N), dtype float
        resultant variable
    lev : int
        level of the variable (vertical coordinate), int
    it : int
        time step of the variable, int
    *variables: Tuple[np.ndarray, ...]
        variable names, all of shape (T, M), dtype float

    Returns
    -------
    np.ndarray, shape (M, N), dtype float
        resultant variable, log10(sum(variables))
    """
    assert len(lev) == 1, 'transformSumLog10 operator can only have one level'
    val[:] = 0.
    for var in variables:
        # Check if the variable is not None
        if var is None:
            raise ValueError("Variable cannot be None")
        # Check if the indices are within bounds
        if it >= var.shape[0] or lev[0] >= var.shape[1]:
            raise IndexError("Index out of bounds")
        val = val + var[it, lev[0]]

    # Check if the value is non-negative
    val[val <= 0.0] = -14.0
    # Apply log10 to positive values only
    val[val > 0.0] = np.log10(val[val > 0.0])

    return val


def transform_sum(val: np.ma.MaskedArray, lev: np.ndarray,
                  it: int, *variables) -> np.ma.MaskedArray:
    """Summation over variables

    Parameters
    ----------
    val : np.ndarray, shape (M, N), dtype float
        resultant variable
    lev : int
        level of the variable (vertical coordinate), int
    it : int
        time step of the variable, int
    *variables: Tuple[np.ndarray, ...]
        variable names, all of shape (T, M), dtype float

    Returns
    -------
    np.ndarray, shape (M, N), dtype float
        resultant variable, log10(sum(variables))
    """
    assert len(lev) == 1, 'transformSum operator can only have one level'
    val[:] = 0
    for var in variables:
        # Check if the variable is not None
        if var is None:
            raise ValueError("Variable cannot be None")
        # Check if the indices are within bounds
        if it >= var.shape[0] or lev[0] >= var.shape[1]:
            raise IndexError("Index out of bounds")
        # to retain the masked array property, do not use +=
        val = val + var[it, lev[0]]
    return val


def identity(val: np.ma.MaskedArray, lev: np.ndarray,
             it: int, *variables) -> np.ma.MaskedArray:
    """Return the val by given `lev` and `it` for 3d variables,
    or `it` for 2d variables.

    Parameters
    ----------
    val : np.ma.MaskedArray, shape (M, N), dtype float
        resultant variable
    lev : int
        level of the variable (vertical coordinate), int
    it : int
        time step of the variable, int
    *variables: Tuple[netCDF4.Variable, ...]
        netCDF4 variable handlers

    Returns
    -------
    np.ma.MaskedArray, shape (M, N), dtype float
        variable[it, lev] for 3d variables
    """
    assert len(variables) == 1, 'identity operator can only have one varname'
    var = variables[0]
    if var is None:
        raise ValueError("Variable cannot be None")
    if len(var.shape) == 3:
        val[:] = var[it]
    else:
        val[:] = var[it, lev]
    return val


def get_lognormal_mean_std(
        mu: np.ma.MaskedArray, sigma: np.ma.MaskedArray) -> tuple[
        np.ma.MaskedArray, np.ma.MaskedArray]:
    """Return the expectation and standard deviation of
    a base 10 lognormal distribution from its parameters.

    Parameters
    ----------
    mu :np.ma.MaskedArray, shape (M,), dtype float
        mean values of the corresponding normal distribution
    sigma : np.ma.MaskedArray, shape (M,), dtype float
        standard deviation of the corresponding normal distribution

    Returns
    -------
    mean : np.ma.MaskedArray, shape (M,), dtype float
        expactation of the base 10 lognormal distribution
    std : np.ma.MaskedArray, shape (M,), dtype float
        standard deviation of the base 10 lognormal distribution
    """
    mean: np.ma.MaskedArray = 10**(mu + 0.5*np.log(10)*sigma*sigma)
    std: np.ma.MaskedArray = np.sqrt(np.exp((np.log(10)*sigma)**2) - 1)
    std = mean*std
    return mean, std


def get_lognormal_params(mean: np.ma.MaskedArray,
                         std: np.ma.MaskedArray
                         ) -> tuple[np.ma.
                                    MaskedArray, np.ma.MaskedArray]:
    """Compute parameters of a based 10 lognormal distribution from
    its mean and standard deviation.

    Parameters
    ----------
    mean : np.ma.MaskedArray, shape (M,), dtype float
        Mean values of lognormal distribution.
    std : np.ma.MaskedArray, shape (M,), dtype float
        Standard deviations of lognormal distribution.

    Returns
    -------
    mu : np.ma.MaskedArray, shape (M,), dtype float
        Mean values of normal distribution.
    sigma: np.ma.MaskedArray, shape (M,), dtype float
        Standard deviations of normal distribution.
    """
    sigma: np.ma.MaskedArray = std*std/mean/mean
    sigma = (np.log(sigma) + 1)/np.log(10)/np.log(10)

    mu: np.ma.MaskedArray = np.log10(mean) - 0.5*np.log(10)*sigma
    sigma = np.sqrt(sigma)

    return mu, sigma


def get_lognormal_mode(mu: np.ma.MaskedArray,
                       sigma: np.ma.MaskedArray
                       ) -> np.ma.MaskedArray:
    """Compute the mode of a base 10 lognormal distribution
    from its parameters.

    Parameters
    ----------
    mean : np.ma.MaskedArray, shape (M,), dtype float
        Mean values of the normal distribution.
    std : np.ma.MaskedArray, shape (M,), dtype float
        Standard deviations of the normal distribution.

    Returns
    -------
    np.ma.MaskedArray, shape (M,), dtype float
        Mode values.
    """
    return 10**(mu - np.log(10)*sigma*sigma)
