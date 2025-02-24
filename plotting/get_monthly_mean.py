"""Calculate monthly climatology of model and observation data.
"""
import multiprocessing as mp
import os
import typing

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.gridspec as mgs
import matplotlib.colors as mcolors
import cartopy  # type: ignore
import cartopy.crs as ccrs  # type: ignore
import cmocean  # type: ignore

import config
import diag
import file_info
import utils


def get_monthly_model(year: str, month: str, varname: str,
                      input_filename: str) -> None:
    """Calculate monthly climatology of variable from
    ensemble mean and save to file.

    Parameters
    ----------
    year : str
        year of model data
    month : str
        month of model data
    varname : str
        name of the variable in output filename.
    filename : str
        Name of the npz file from `get_ensmean` module.
    """
    monthly_model_data: np.ma.MaskedArray
    # Calculate the number of days in the month
    days_in_month = utils.get_month_days(int(year), int(month))
    # Initialize arrays to store monthly model and observation data
    monthly_model_data = np.ma.zeros((config.ny, config.nx),
                                     dtype=float)

    # Calculate the monthly mean
    rng = range(1, days_in_month + 1)
    for day in rng:
        # Calculate the ensemble mean for the model data
        f = np.load(
            input_filename.format(
                year=year, month=month, day=str(day).zfill(2)))
        monthly_model_data = monthly_model_data + \
            np.ma.masked_array(f['model'], f['model_mask']) / \
            days_in_month

    # Create the output directory if it doesn't exist
    os.makedirs(os.path.join(
        'data', config.exp), exist_ok=True)

    # Save the monthly climatology to a numpy archive
    np.savez(
        os.path.join('data',
                     config.exp,
                     f'monthly_{year}{month}_model_{varname}.npz'
                     ),
        model=monthly_model_data.data,
        model_mask=monthly_model_data.mask
    )


def get_monthly_ens_pdaf(year: str, month: str,
                         varname: str, i_ens: int
                         ) -> None:
    """Calculate (monthly) climatology of PDAF output
    for a given ensemble member `i_ens` and save to file.

    Parameters
    ----------
    year : str
        PDAF output year.
    month: str
        PDAF output month.
    varname: str
        Name of the variable.
    i_ens: int
        Ensemble index.
    """
    model_data: np.ma.MaskedArray
    # Calculate the number of days in the month
    days_in_month = utils.get_month_days(int(year), int(month))
    # Initialize arrays to store monthly model and observation data
    monthly_model_data = np.ma.zeros((config.ny, config.nx),
                                     dtype=float)

    f_info: file_info.FileInfo
    f_info = file_info.FileInfo(year, month, '00', '')
    if f_info.stats_pdaf:
        return

    rng = range(1, days_in_month + 1)
    if config.exp in ['PC', 'chlo-monthly', 'chlo-pc']:
        rng = range(days_in_month//2, days_in_month//2 + 1)
    for day in rng:
        f_info.day = str(day).zfill(2)
        model_data = diag.read_pdaf_output_3d(
            f_info.get_pdaf_filename(varname).format(i_ens),
            varname)
        monthly_model_data = monthly_model_data + \
            model_data / days_in_month

    # Create the output directory if it doesn't exist
    os.makedirs(os.path.join('data', config.exp), exist_ok=True)

    # Save the monthly climatology to a numpy archive
    np.savez(
        os.path.join(
            'data',
            config.exp,
            f'pdaf_monthly_{year}{month}_ens_{i_ens}_{varname}.npz'),
        model=monthly_model_data.data,
        model_mask=monthly_model_data.mask
    )


def get_monthly_obs(
        year: str, month: str, obs_type: typing.Literal
        ['chlo-monthly', 'pc']) -> None:
    """Obtain monthly composite of observations in log10 space.

    Parameters
    ----------
    year : str
        observation year.
    month : str
        observation month.
    obs_type : Literal['chlo-monthly', 'pc']
        Name of the observation type. It can be chlo-monthly or pc.
    """
    # Initialize arrays to store monthly model and observation data
    monthly_obs_data: np.ma.MaskedArray = np.ma.zeros(
        (config.ny, config.nx), dtype=float)
    f_info: file_info.FileInfo
    f_info = file_info.FileInfo(year, month, '01', '')
    monthly_obs_data, _ = diag.read_observation(
        f_info.get_obsfilename(obs_type), obs_type)
    # Create the output directory if it doesn't exist
    os.makedirs('monthly_obs', exist_ok=True)
    # Save the monthly climatology to a numpy archive
    np.savez(
        os.path.join(
            'monthly_obs',
            f'monthly_{year}{month}_obs_{obs_type}.npz'
        ),
        obs=monthly_obs_data.data,
        obs_mask=monthly_obs_data.mask
    )


def get_monthly_obs_physical(
        year: str, month: str, obs_type: typing.Literal
        ['chlo-monthly', 'pc']) -> None:
    """Obtain monthly composite of observations in physical space.

    Parameters
    ----------
    year : str
        observation year.
    month : str
        observation month.
    obs_type : Literal['chlo-monthly', 'pc']
        Name of the observation type. It can be chlo-monthly or pc.
    """
    # Initialize arrays to store monthly model and observation data
    monthly_obs_data: np.ma.MaskedArray = np.ma.zeros(
        (config.ny, config.nx), dtype=float)
    f_info: file_info.FileInfo
    f_info = file_info.FileInfo(year, month, '01', '')
    monthly_obs_data, _ = diag.read_physical_observation(
        f_info.get_obsfilename(obs_type), obs_type)
    # Create the output directory if it doesn't exist
    os.makedirs('monthly_obs', exist_ok=True)
    # Save the monthly climatology to a numpy archive
    np.savez(
        os.path.join(
            'monthly_obs',
            f'monthly_{year}{month}_obs_{obs_type}_phys.npz'
        ),
        obs=monthly_obs_data.data,
        obs_mask=monthly_obs_data.mask
    )


def save_monthly_multiprocessed(func: typing.Callable,
                                kwargs: dict) -> None:
    """Save monthly mean over multiple years and months
    using multiprocessing.

    Parameters
    ----------
    func : typing.Callable
        Function to be called.
    kwargs: dict
        Dictionary of keyword arguments to be passed to
          the function, `func`.
        This depends on the function being called.
        The input should refer to individual functions.
    """

    years: list[str] = ['2015', '2016']
    months: list[str] = [str(month).zfill(2)
                         for month in range(1, 13)]

    for year in years:
        processes: list[mp.Process] = []
        process: mp.Process
        for month in months:
            if year == '2015' and month == '01':
                continue
            kwargs['year'] = year
            kwargs['month'] = month
            process = mp.Process(target=func, kwargs=kwargs)
            processes.append(process)
            process.start()

        for process in processes:
            process.join()

        if config.exp in ['PC', 'chlo-pc',]:
            break


def plot_exp_from_free(
        year: str, month: str, varnames: list[str]) -> None:
    """Plot the deviation of model from free run
    """

    lons: np.ma.MaskedArray
    lats: np.ma.MaskedArray
    lons, lats = utils.get_coord()

    exps: list[str] = ['chlo', 'chlo-monthly', 'chlo-monthly-update',
                       'PC', ]
    exp_labels: list[str] = ['Daily Chl', 'Monthly Chl',
                             'Monthly Chl+', 'Monthly C', ]

    for exp, explabel in zip(exps, exp_labels):
        os.makedirs(f'exp_model_deviation/{exp}', exist_ok=True)
        fig: matplotlib.figure.Figure = plt.figure()
        fig.clf()
        w, h = fig.get_size_inches()
        fig.set_size_inches(w*2, h)
        gs = mgs.GridSpec(1, 2, figure=fig, wspace=0.01,
                          left=0., right=1., bottom=0.0, top=1.)
        for i, varname in enumerate(varnames):
            ax: cartopy.mpl.geoaxes.GeoAxes = fig.add_subplot(
                gs[i], projection=ccrs.Robinson())

            # read the monthly model exp ensemble mean
            f: np.lib.npyio.NpzFile = np.load(os.path.join(
                'data',
                exp, f'monthly_{year}{month}_model_{varname}.npz')
            )
            model: np.ma.MaskedArray = np.ma.masked_array(
                f['model'], f['model_mask']
            )
            # read the monthly model free run ensemble mean
            f = np.load(
                os.path.join(
                    'data',
                    'free',
                    f'monthly_{year}{month}_model_{varname}.npz'
                )
            )
            model_free: np.ma.MaskedArray = np.ma.masked_array(
                f['model'], f['model_mask']
            )

            ax.set_global()
            pc = ax.pcolormesh(
                lons, lats, np.squeeze(model - model_free),
                # pylint: disable=no-member
                cmap=cmocean.cm.balance, norm=mcolors.CenteredNorm(
                    vcenter=0., halfrange=0.1),
                transform=ccrs.PlateCarree())
            ax.coastlines(color='k', linewidth=.8)
            ax.set_title(f'{year}-{month} {varname} ({explabel} - Free)')

            fig.colorbar(pc, ax=ax, orientation='horizontal',
                         shrink=0.6, pad=0.02)
        fig.savefig(
            f'exp_model_deviation/{exp}/bias_{year}{month}.png',
            dpi=300)
        plt.close(fig)


if __name__ == '__main__':
    # get ensemble mean of model output climatology
    vnames = ['chlorophyll', 'nitrogen']
    # for expname in ['chlo', 'PC', 'chlo-monthly', 'free',
    #                 'chlo-monthly-update']:
    for expname in ['chlo-pc']:
        config.exp = expname
        config.output_path = config.output_path_format.format(
            exp=expname)
        for vname in vnames:
            monthly_model_kwargs = {'varname': vname,
                                    'input_filename':
                                    os.path.join(
                                        'data', config.exp,
                                        f'ensmean_{vname}_'
                                        '{year}{month}{day}.npz'
                                    )
                                    }
            save_monthly_multiprocessed(
                get_monthly_model, monthly_model_kwargs)

    # for year in ['2015', ]:
    #     for month in [str(month).zfill(2) for month in range(2, 13)]:
    #         plot_exp_from_free(year, month, ['chlorophyll', 'nitrogen'])
