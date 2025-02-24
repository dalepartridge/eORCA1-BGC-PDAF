"""This module calculates the 5-year climatology of the chlorophyll
and nitrogen variables for the eORCA1-BGC-PDAF model output.

The output can be used to calculate the anomalies of the chlorophyll
and nitrogen variables for the ACC calculation.
"""

import numpy as np
import netCDF4
import os
import config
import utils
import typing


def get_5year_climatology() -> None:
    """Calculate the 5-year climatology for a given variable and month.

    Args:
        years (List[str]): The years of the PDAF output.
        months (List[str]): The months of the PDAF output.

    Returns:
        np.ma.MaskedArray: The 5-year climatology for the given variable and month.
    """
    clim_CHD: np.ma.MaskedArray = np.ma.zeros(
        (config.ny, config.nx), dtype=float)
    clim_CHN: np.ma.MaskedArray = np.ma.zeros(
        (config.ny, config.nx), dtype=float)
    clim_PHD: np.ma.MaskedArray = np.ma.zeros(
        (config.ny, config.nx), dtype=float)
    clim_PHN: np.ma.MaskedArray = np.ma.zeros(
        (config.ny, config.nx), dtype=float)

    years: typing.Iterable = range(2000, 2005)
    months: typing.Iterable = range(1, 13)
    ndays: int = 0
    history_path: str = '/work/n01/n01/dapa/NCEO/eORCA1-BGC-PDAF/OUTPUTS/'
    for year in years:
        for month in months:
            ndays_in_month = utils.get_month_days(year, month)
            filename = 'eORCA1_1d_{year}{month}01_{year}{month}{day}_ptrc1_T.nc'.format(
                year=year, month=str(month).zfill(2), day=ndays_in_month)
            f: netCDF4.Dataset = netCDF4.Dataset(os.path.join(
                history_path, year, month.zfill(2), filename), 'r')
            CHD: np.ma.MaskedArray = f.variables['CHD'][:, 0, :, :].sum(
                0)  # Diatom chlorophyll
            CHN: np.ma.MaskedArray = f.variables['CHN'][:, 0, :, :].sum(
                0)  # Non-diatom chlorophyll
            PHD: np.ma.MaskedArray = f.variables['PHD'][:, 0, :, :].sum(
                0)  # Diatom nitrogen
            PHN: np.ma.MaskedArray = f.variables['PHN'][:, 0, :, :].sum(
                0)  # Non-diatom nitrogen
            f.close()
            clim_CHD = clim_CHD + CHD
            clim_CHN = clim_CHN + CHN
            clim_PHD = clim_PHD + PHD
            clim_PHN = clim_PHN + PHN
            ndays = ndays + ndays_in_month
    clim_CHD = clim_CHD / ndays
    clim_CHN = clim_CHN / ndays
    clim_PHD = clim_PHD / ndays
    clim_PHN = clim_PHN / ndays

    np.savez('5year_clim.npz', clim_CHD=clim_CHD.data, clim_CHN=clim_CHN.data, clim_PHD=clim_PHD.data, clim_PHN=clim_PHN.data,
             CHD_mask=clim_CHD.mask, CHN_mask=clim_CHN.mask, PHD_mask=clim_PHD.mask, PHN_mask=clim_PHN.mask)


if __name__ == '__main__':
    get_5year_climatology()
