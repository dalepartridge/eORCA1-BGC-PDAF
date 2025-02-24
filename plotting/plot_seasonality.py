"""Comaprison of ratio between phytoplankton constituents.
"""
import itertools
import os

import cartopy
import cartopy.mpl.geoaxes  # type: ignore
import cartopy.crs as ccrs
import cartopy.feature as cfeature  # type: ignore
import cmocean  # type: ignore # pylint: disable=unused-import
import matplotlib.colors as mcolors
import matplotlib.gridspec as mgs
import matplotlib.pyplot as plt
import numpy as np

import config
import utils
import file_info
import diag


plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman'] + \
    plt.rcParams['font.serif']
plt.rcParams['font.size'] = 12


def save_seasonal_model_data(vname: str) -> None:
    """Sum up data seasonally from monthly data.

    Parameters
    ----------
    vname : str
        The variable name.
    """
    years: list[str]
    years = ['2015', '2016']
    seasons = {
        'winter': ['12', '01', '02'],
        'spring': ['03', '04', '05'],
        'summer': ['06', '07', '08'],
        'autumn': ['09', '10', '11']
    }

    data: dict[str, np.ndarray] = {season: np.zeros(
        (config.ny, config.nx)) for season in seasons}
    n_months: dict[str, int]
    n_months = {season: len(years) * 3 for season in seasons}
    n_months['winter'] = n_months['winter'] - 1

    for year in years:
        for season, months in seasons.items():
            for month in months:
                if year == '2015' and month == '01':
                    continue
                f = np.load(
                    os.path.join(
                        'data', config.exp,
                        'ensmean_'f'{vname}_{year}{month}01.npz'
                    )
                )
                data[season] = data[season] + np.ma.masked_array(
                    f['model'], f['model_mask'])/n_months[season]
    np.savez(os.path.join('data', config.exp, f'season_{vname}_map.npz'),
             **data)  # type: ignore


def save_seasonal_obs_data(vname: str) -> None:
    """Sum up data seasonally from monthly observation data.

    Parameters
    ----------
    vname : str
        The variable name.
    """
    years: list[str]
    years = ['2015', '2016']
    seasons = {
        'winter': ['12', '01', '02'],
        'spring': ['03', '04', '05'],
        'summer': ['06', '07', '08'],
        'autumn': ['09', '10', '11']
    }

    data: dict[str, np.ndarray] = {season: np.zeros(
        (config.ny, config.nx)) for season in seasons}
    valid_counts: dict[str, np.ndarray] = {season: np.zeros(
        (config.ny, config.nx), dtype=int) for season in seasons}

    for year in years:
        for season, months in seasons.items():
            for month in months:
                if year == '2015' and month == '01':
                    continue
                f_info = file_info.FileInfo(year, month, '01', '')
                obs_data, _ = diag.read_physical_observation(
                    f_info.get_obsfilename(vname), vname)
                obs_mask = np.logical_not(obs_data.mask)
                data[season] += np.where(obs_mask, obs_data, 0)
                valid_counts[season] += obs_mask.astype(int)

    # Calculate the temporal average
    for season in seasons:
        valid_mask = valid_counts[season] > 0
        data[season] = np.where(
            valid_mask, data[season] / valid_counts[season],
            np.nan)
    if vname == 'pc':
        vname = 'nitrogen-monthly'
    # Save the seasonal data to a file or process it as needed
    np.savez(
        os.path.join('data', 'obs', f'season_{vname}_map.npz'),
        **data)  # type: ignore


def plot_axes(
        ax: cartopy.mpl.geoaxes.GeoAxes, vname: str, exp: str, season: str):
    """plot specific axis for the seasonal map.
    """
    norm_free: dict[str, mcolors.Normalize] = {}
    norm_free['chlo-monthly'] = mcolors.TwoSlopeNorm(
        vcenter=0., vmin=-0.35, vmax=0.35)
    norm_free['nitrogen-monthly'] = mcolors.TwoSlopeNorm(
        vcenter=0., vmin=-0.5, vmax=0.5)
    norm_free['OCN_PCO2'] = mcolors.TwoSlopeNorm(
        vcenter=0., vmin=-20., vmax=20)
    norm_free['OXY'] = mcolors.TwoSlopeNorm(vcenter=0., vmin=-20., vmax=20)
    norm_free['PHD'] = mcolors.TwoSlopeNorm(vcenter=0., vmin=-0.3, vmax=0.3)
    norm_free['PHN'] = mcolors.TwoSlopeNorm(vcenter=0., vmin=-0.1, vmax=0.1)
    norm_free['CHD'] = mcolors.TwoSlopeNorm(vcenter=0., vmin=-0.1, vmax=0.1)
    norm_free['CHN'] = mcolors.TwoSlopeNorm(vcenter=0., vmin=-.1, vmax=0.1)
    norm_free['PRD'] = mcolors.TwoSlopeNorm(vcenter=0., vmin=-0.8, vmax=0.8)
    norm_free['PRN'] = mcolors.TwoSlopeNorm(vcenter=0., vmin=-1.5, vmax=1.5)

    norm_diff: dict[str, mcolors.Normalize] = {}
    norm_diff['chlo-monthly'] = mcolors.TwoSlopeNorm(
        vcenter=0., vmin=-0.06, vmax=0.06)
    norm_diff['nitrogen-monthly'] = mcolors.TwoSlopeNorm(
        vcenter=0., vmin=-0.08, vmax=0.08)
    if exp == 'chlo-monthly-update':
        norm_diff['nitrogen-monthly'] = mcolors.TwoSlopeNorm(
            vcenter=0., vmin=-0.03, vmax=0.03)
    if exp == 'PC-update':
        norm_diff['chlo-monthly'] = mcolors.TwoSlopeNorm(
            vcenter=0., vmin=-0.3, vmax=0.3)
    norm_diff['OCN_PCO2'] = mcolors.TwoSlopeNorm(
        vcenter=0., vmin=-1.2, vmax=1.2)
    norm_diff['OXY'] = mcolors.TwoSlopeNorm(vcenter=0., vmin=-0.2, vmax=0.2)
    norm_diff['PHD'] = mcolors.TwoSlopeNorm(vcenter=0., vmin=-0.02, vmax=0.02)
    norm_diff['PHN'] = mcolors.TwoSlopeNorm(vcenter=0., vmin=-0.01, vmax=0.01)
    norm_diff['CHD'] = mcolors.TwoSlopeNorm(vcenter=0., vmin=-0.03, vmax=0.03)
    norm_diff['CHN'] = mcolors.TwoSlopeNorm(vcenter=0., vmin=-.03, vmax=0.03)
    norm_diff['PRD'] = mcolors.TwoSlopeNorm(vcenter=0., vmin=-0.2, vmax=0.2)
    norm_diff['PRN'] = mcolors.TwoSlopeNorm(vcenter=0., vmin=-0.5, vmax=0.5)

    seasons: dict[str, str] = {
        'DJF': 'winter',
        'MAM': 'spring',
        'JJA': 'summer',
        'SON': 'autumn'
    }

    lons, lats = utils.get_coord()

    f = np.load(os.path.join('data', exp, f'season_{vname}_map.npz'))
    if vname in ['chlo-monthly', 'nitrogen-monthly']:
        f_clim = np.load(os.path.join(
            'data', exp, f'clim_{vname}.npz'))
        clim = f_clim['clim'].reshape(config.ny, config.nx)

        f_clim_free = np.load(os.path.join(
            'data', 'free', f'clim_{vname}.npz'))
        clim_free = f_clim_free['clim'].reshape(config.ny, config.nx)
    else:
        clim = np.load(os.path.join(
            'data', exp, f'flux_{vname}_map.npz'))['data']

        clim_free = np.load(os.path.join(
            'data', 'free', f'flux_{vname}_map.npz'))['data']

    if exp in ['obs', 'free']:
        d = np.squeeze(f[seasons[season]]) - clim
        norm = norm_free[vname]
    else:
        f_free = np.load(os.path.join('data', 'free',
                                      f'season_{vname}_map.npz'))
        d = np.squeeze(f[seasons[season]] - clim
                       ) - np.squeeze(f_free[seasons[season]] - clim_free)
        norm = norm_diff[vname]
    print(np.nanmin(d), np.nanmax(d), np.nanpercentile(d, [5, 95]))
    pc = ax.pcolormesh(
        lons, lats, np.squeeze(d), transform=ccrs.PlateCarree(),
        cmap='cmo.balance', norm=norm)

    ax.coastlines(color='k', linewidth=.8)
    ax.add_feature(cfeature.LAND, zorder=3)

    return pc, norm


def plot_diff_map(vname: str) -> None:
    """plot the differences between DA experiment and free run"""
    os.makedirs('figs', exist_ok=True)

    exps: list[str] = ['obs', 'free', 'chlo-monthly-update',
                       'PC-update', 'chlo-pc',
                       ]

    fig: plt.Figure = plt.figure()
    w, h = fig.get_size_inches()
    fig.set_size_inches(w*1.9, h*1.2)

    gs: mgs.GridSpec = mgs.GridSpec(4, len(exps), figure=fig,
                                    wspace=0.0,
                                    hspace=0.04,
                                    left=0.03, right=1.,
                                    bottom=0.1, top=0.94)

    ax: cartopy.mpl.geoaxes.GeoAxes
    season_exp = itertools.product(
        ['DJF', 'MAM', 'JJA', 'SON'], exps)
    for i, (season, exp) in enumerate(season_exp):
        ax = fig.add_subplot(gs[i], projection=ccrs.Robinson())
        ax.set_global()
        pc, norm = plot_axes(ax, vname, exp, season)
        if season == 'DJF':
            ax.set_title(config.exp_labels[exp])
        if exp == 'obs':
            ax.text(0.005, 0.82 - 0.22*(i//5), season, rotation=90,
                    transform=fig.transFigure)
        if season == 'SON':
            ax = fig.add_axes((0.065 + 0.197*(i - 15), 0.07, 0.12, 0.02))
            fig.colorbar(pc, cax=ax, orientation='horizontal', pad=0.05,
                         shrink=0.5, norm=norm)

    fig.savefig(f'figs/season_{vname}_diff_map.png', dpi=300)


if __name__ == '__main__':
    expnames: list[str] = ['free', 'chlo', 'chlo-monthly',
                           'chlo-monthly-update', 'PC',
                           'PC-update', 'chlo-pc',]
    variables = ['chlo-monthly', 'nitrogen-monthly',]
    # variables += ['OCN_PCO2', 'OXY', ]
    # variables += ['PHD', 'PHN', 'CHD', 'CHN']
    # variables += ['PRD', 'PRN']
    # for expname in expnames:
    #     config.exp = expname
    #     config.output_path = config.output_path_format.format(exp=expname)
    #     for varname in variables[7:]:
    #         save_seasonal_model_data(varname)
    #         if varname == 'chlo-monthly':
    #             save_seasonal_obs_data(varname)
    #         if varname == 'nitrogen-monthly':
    #             save_seasonal_obs_data('pc')

    for varname in variables[:1]:
        print(varname)
        plot_diff_map(varname)
