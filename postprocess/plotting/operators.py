import numpy as np


def transformSumLog10(
    val: np.ndarray,  # type: np.ndarray
    lev: int,  # type: int
    it: int,  # type: int
    *variables: np.ndarray  # type: Tuple[np.ndarray, ...]
) -> np.ndarray:
    """Summation over variables and log10 the summed variable

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
    val = 0.0
    for var in variables:
        # Check if the variable is not None
        if var is None:
            raise ValueError("Variable cannot be None")
        # Check if the indices are within bounds
        if it >= var.shape[0] or lev >= var.shape[1]:
            raise IndexError("Index out of bounds")
        val += var[it, lev]

    # Check if the value is non-negative
    val[val <= 0.0] = -14.0
    # Apply log10 to positive values only
    val[val > 0.0] = np.log10(val[val > 0.0])

    return val


def transformSum(
    val: np.ndarray,  # type: np.ndarray
    lev: int,  # type: int
    it: int,  # type: int
    *variables: np.ndarray  # type: Tuple[np.ndarray, ...]
) -> np.ndarray:
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
    val = 0.0
    for var in variables:
        # Check if the variable is not None
        if var is None:
            raise ValueError("Variable cannot be None")
        # Check if the indices are within bounds
        if it >= var.shape[0] or lev >= var.shape[1]:
            raise IndexError("Index out of bounds")
        val += var[it, lev]

    return val



def identity(val: np.ndarray, lev: int, it: int, *variables: np.ndarray) -> np.ndarray:
    """Summation over variables and log10 the summed variable

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
        variable[it, lev]
    """
    assert len(variables) == 1, 'identity operator can only have one varname'
    var = variables[0]
    if var is None:
        raise ValueError("Variable cannot be None")
    if it >= var.shape[0] or lev >= var.shape[1]:
        raise IndexError("Index out of bounds")
    return var[it, lev]

def unlog_observation(obs: np.ndarray, obs_unc: np.ndarray) -> np.ndarray:
    """Unlog observation with given uncertainty.

    Parameters
    ----------
    obs : np.ndarray, shape (M,), dtype float
        Observation values.
    obs_unc : np.ndarray, shape (M,), dtype float
        Observation uncertainties.

    Returns
    -------
    np.ndarray, shape (M,), dtype float
        Unlogged observation values.
    """
    obs = 10**(obs + 0.5*np.log(10)*obs_unc)
    obs_unc = obs*np.sqrt(np.exp((np.log(10)*obs_unc)**2) - 1)
    return obs, obs_unc
