"""Derive filenames based on date and file format.

This module is specifically coded for PDAF-NEMO-MEDUSA experiments.
"""
import os

import config
import utils


class FileInfo:
    """A class storing all information of a file
    """

    def __init__(self, year: str, month: str,
                 day: str, nemo_fileformat: str) -> None:
        self.year: str = year
        self.month: str = month
        self.cycle_month: str = month if int(month) % 2 == 1 else \
            str(int(month) - 1).zfill(2)
        self.day: str = day
        self.stats: bool
        if '__1m__' not in nemo_fileformat:
            self.stats = False
        else:
            self.stats = self.month not in config.monthsEns
        self.stats_pdaf: bool = self.month not in config.monthsEns and \
            config.exp == 'chlo'
        self.nemo_fileformat = nemo_fileformat

    def get_nemo_filename(self) -> str:
        """Get the filename of the NEMO output

        Returns
        -------
        str
            The filename of the NEMO output.
        """
        # get the time period of the experiment cycle
        numdays = utils.get_month_days(
            int(self.year), int(self.cycle_month) + 1)
        filerange: str
        if self.year == '2015' and self.cycle_month == '01':
            filerange = f'{self.year}' \
                f'{str(int(self.cycle_month) + 1).zfill(2)}01_'\
                f'{self.year}' \
                f'{str(int(self.cycle_month) + 1).zfill(2)}{numdays}'
        else:
            filerange = f'{self.year}{self.cycle_month}01_' \
                f'{self.year}' \
                f'{str(int(self.cycle_month) + 1).zfill(2)}{numdays}'
        # get the nemo filename
        fname = self.nemo_fileformat.format(filerange=filerange,
                                            year=self.year,
                                            month=self.month)
        # combine the path and filename
        if self.stats:
            return os.path.join(
                config.output_path, self.year, self.cycle_month, fname)

        return os.path.join(
            config.output_path, self.year, self.cycle_month,
            'ensemble_{0}', fname)

    def get_obsfilename(self, obs_type: str) -> str:
        """Get the filename of the observations.

        Parameters
        ----------
        obs_type : str
            Type of the observations. It can be 'chlo', 'pc', 'chlo-monthly'.

        Returns
        -------
        str
            The filename of the observations.
        """
        if obs_type == 'chlo':
            return os.path.join(
                config.OBS_PATH, f"chlo_{self.year}{self.month}{self.day}.nc")
        elif obs_type == 'pc':
            return os.path.join(
                config.OBS_PATH, f"pc_{self.year}{self.month}.nc")
        elif obs_type == 'chlo-monthly':
            return os.path.join(
                config.OBS_PATH, f"chlo_{self.year}{self.month}.nc")
        else:
            raise ValueError(f"Unknown observation type: {obs_type}")

    def get_pdaf_filename(self, varname: str) -> str:
        """Get the filename of the PDAF output

        Parameters
        ----------
        varname : str
            the variable name.

        Returns
        -------
        str
            The filename of the PDAF output.
        """
        if self.stats_pdaf:
            return os.path.join(
                config.output_path, self.year,
                self.cycle_month, 'analysis',
                f'state_{varname}_{self.year}{self.month}{self.day}.nc')
        else:
            return os.path.join(config.output_path,
                                self.year, self.cycle_month, 'analysis',
                                f'state_{varname}_{self.year}{
                                    self.month}{self.day}_'
                                '{0}.nc')
