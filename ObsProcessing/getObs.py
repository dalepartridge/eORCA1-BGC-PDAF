"""This file is used to process BGC observation data

It links each model grid with a set of observations
then the median of observations are used as superobing. 
"""
import xarray as xr
import numpy as np
import mod_gridsearch
import mod_obs
import mpi4py.MPI as MPI


def writeFineGridModelIndexofObs(filename, outputFilename, lonName, latName):
    """This function links the model grid with observations.
    """
    fo = xr.open_dataset(filename)
    fm = xr.open_dataset('state_20150116.nc')
    lons_m = fm['nav_lon']
    lats_m = fm['nav_lat']
    ny, nx = lons_m.shape

    # double the model grid resolution
    ny, nx = 2*ny - 1, 2*nx - 1
    mlon = np.zeros((ny, nx))
    mlat = np.zeros((ny, nx))
    mmask = np.ones_like(mlon)
    mlon[::2, ::2] = lons_m[:]
    mlon[1::2, ::2] = 0.5*(lons_m[1:, :] + lons_m[:-1, :])
    mlon[::2, 1::2] = 0.5*(lons_m[:, 1:] + lons_m[:, :-1])
    mlon[1::2, 1::2] = 0.5*(mlon[1::2, 2::2] + mlon[1::2, :-2:2])

    mlat[::2, ::2] = lats_m[:]
    mlat[1::2, ::2] = 0.5*(lats_m[1:, :] + lats_m[:-1, :])
    mlat[::2, 1::2] = 0.5*(lats_m[:, 1:] + lats_m[:, :-1])
    mlat[1::2, 1::2] = 0.5*(mlat[1::2, 2::2] + mlat[1::2, :-2:2])

    lons_o = fo[lonName]
    lats_o = fo[latName]

    nlons = len(lons_o)
    nlats = len(lats_o)
    mobsi = np.zeros((nlats, nlons), dtype=int)
    mobsj = np.zeros((nlats, nlons), dtype=int)

    for i, lat in enumerate(lats_o.to_numpy()):
        olon = lons_o.to_numpy()
        print (lat)
        olat = lat*np.ones_like(olon)
        mobsi[i], mobsj[i], _ = mod_gridsearch.mod_gridsearch.obs_grd_bruteforce( nx, ny,
                                                                                  1, nx, 1, ny, 1, 1,
                                                                                  mlon.T, mlat.T, mmask.T,
                                                                                  olon, olat)
    ds = xr.Dataset(
        data_vars=dict(
            mobsi=(['ilat', 'ilon'], mobsi),
            mobsj=(["ilat", 'ilon'], mobsj),
        ),
        attrs=dict(description="obs to model indices"),
    )
    ds.to_netcdf(outputFilename)


def getObs(nx, ny, obstype, obsfilename, indexFilename):
    # get land mask
    f = xr.open_dataset('/work/n01/n01/ymchen/eORCA1-BGC-PDAF/INPUTS/PHYSICS/DOM/eORCA_R1_zps_domcfg.nc')
    mask = f['top_level'][0].to_numpy().ravel() < 1e-3
    f.close()

    # get processor information
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nPE = comm.Get_size()

    # decompose domain
    nTotalGp = np.sum(~mask)
    ngridpoints = nTotalGp//nPE

    #ngridpoints = nx*ny//nPE
    istart = rank*ngridpoints
    iend = istart + ngridpoints
    if rank == nPE - 1: iend = nTotalGp
    ngridpoints = iend - istart

    # define model grid index
    model_obs = np.zeros(ngridpoints, dtype='d')
    model_obs_unc = np.zeros(ngridpoints, dtype='d')
    obs_index = np.arange(len(mask))[~mask][istart:iend]
    ind_j = obs_index//nx
    ind_i = obs_index - ind_j*nx
    print (rank, len(obs_index), ngridpoints)

    # get observations
    mod_obs.mod_obs.getobs(obstype, obsfilename, indexFilename,
                           ind_i, ind_j,
                           model_obs=model_obs, model_obs_unc=model_obs_unc,
                           ngridpoints=ngridpoints)
    model_obs[model_obs>1e10] = np.nan

    if obstype == 'pc':
        varname = 'carbon'
    elif obstype == 'chlo':
        varname = 'chlorophyll'
    else:
        print ('unrecognised obstype')


    # gather data into one processor
    sendcounts = np.array(comm.gather(ngridpoints, 0))
    obs_full = None
    obsunc_full = None
    if rank == 0:
        print("sendcounts: {}, total: {}, nTotalGp: {}".format(sendcounts, sum(sendcounts), nTotalGp))
        obs_full=np.empty(nTotalGp, dtype='d')
        obsunc_full=np.empty(nTotalGp, dtype='d')

    comm.Gatherv(sendbuf=model_obs, recvbuf=(obs_full, sendcounts), root=0)
    comm.Gatherv(sendbuf=model_obs_unc, recvbuf=(obsunc_full, sendcounts), root=0)

    # scatter ocean data into the global dataset
    if rank == 0:
        data = np.zeros(len(mask))
        print (varname, mask.shape, data.shape, obs_full.shape)
        data[~mask]  = obs_full
        data[mask]  = 1e20
        obs_full = data.copy().reshape(ny, nx)

        data[~mask]  = obsunc_full
        data[mask]  = 1e20
        obsunc_full = data.copy().reshape(ny, nx)

    return obs_full, obsunc_full


def writeObs(obstype, obs, obsunc, outputFilename):
    if obstype == 'pc':
        varname = 'carbon'
    elif obstype == 'chlo':
        varname = 'chlorophyll'
    else:
        print ('unrecognised obstype')

    data_vars = dict()
    data_vars[varname] = (['y', 'x'], obs)
    data_vars[varname+'_unc'] = (['y', 'x'], obsunc)

    ds = xr.Dataset(
        data_vars,
        attrs=dict(description=f'processed observation on {year}{month}'),
    )

    encoding = dict()
    encoding[varname] = {'_FillValue': 1e20,}
    encoding[varname+'_unc'] = {'_FillValue': 1e20,}
    ds.to_netcdf(outputFilename, encoding=encoding)


if __name__ == '__main__':
    year = '2015'
    days = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    for month in [9, ]: # range(5, 9):
        # for day in range(1, days[month - 1] + 1):
        # for day in range(1, 15):
        for day in range(15, 31):
            MPI.COMM_WORLD.barrier()
            obstype = 'chlo'
            # obsfilename = f'/work/n01/n01/ymchen/obs/BICEP_NCEO_PC_ESA-OC-L3S-MERGED-1M_MONTHLY_9km_mapped_{year}{str(month).zfill(2)}-fv5.nc'
            obsfilename = f'/work/n01/n01/ymchen/obs/ESACCI-OC-L3S-CHLOR_A-MERGED-1D_DAILY_4km_GEO_PML_OCx-{year}{str(month).zfill(2)}{str(day).zfill(2)}-fv5.0.nc'
            obs, obsunc = getObs(362, 332, obstype, obsfilename, f'obs_model_indice-{obstype}.nc')
            if  MPI.COMM_WORLD.Get_rank() == 0:
                writeObs(obstype, obs, obsunc, f'log-{obstype}_{year}{str(month).zfill(2)}{str(day).zfill(2)}.nc')
