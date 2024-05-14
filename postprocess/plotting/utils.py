import calendar

import netCDF4
import numpy as np

import config


def get_month_days(year: int, month: int) -> int:
    """Get the number of days in a specific month of a specific year.

    Parameters
    ----------
    year : int
        The year.
    month : int
        The month.

    Returns
    -------
    int
        The number of days in the month.
    """

    # Get the number of days in the month using the calendar module
    return calendar.monthrange(year, month)[1]


def get_land_mask() -> np.ndarray:
    """Get mask of land grid cells.

    Returns
    -------
    numpy.ndarray
        Boolean array of mask.
    """
    with netCDF4.Dataset(config.DOMAIN_PATH, 'r') as dataset:
        top_level_array = dataset['top_level'][0]
    return (top_level_array < 1e-3).data


def get_coord() -> tuple[np.ndarray, np.ndarray]:
    """Get coordinates of the domain.

    Returns:
        tuple(numpy.ndarray, numpy.ndarray): Longitude and latitude arrays.
    """
    with netCDF4.Dataset(config.DOMAIN_PATH, 'r') as f:
        lons: np.ndarray = wrap_data(f['nav_lon'][:])
        lats: np.ndarray = f['nav_lat'][:]
    return lons, lats


def wrap_data(lons: np.ndarray) -> np.ndarray:
    """Wrap longitude data.

    Args:
        lons (numpy.ndarray): Array of longitude values.

    Returns:
        numpy.ndarray: Array of wrapped longitude values.
    """
    fixed_lons = lons.copy()
    for i, start in enumerate(np.argmax(np.abs(np.diff(lons)) > 180, axis=1)):
        fixed_lons[i, start+1:] += 360
    return fixed_lons
