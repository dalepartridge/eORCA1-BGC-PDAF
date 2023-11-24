from mpi4py import MPI
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import mod_obs


_Ne = 30
_range = dict()
_range['01'] = ('0101', '0331')
_range['04'] = ('0401', '0630')
_range['07'] = ('0701', '0930')
_range['10'] = ('1001', '1231')


class RankHistogram:
    def __init__(self, obstype, expname, suffix, year, month, cycle_month, varnames=['PHD', 'PHN']):
        # get parallel information
        self.comm = MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.nPE = self.comm.Get_size()
        self.getMask()

        # get model fields
        filename = f'eORCA1_1m_{year}{_range[cycle_month][0]}_{year}' \
                    f'{_range[cycle_month][1]}_ptrc_T_{year}{month}-{year}{month}.nc'
        fens = [xr.open_dataset(f'{expname}/{year}{suffix}/{cycle_month}/ensemble_{i}/{filename}') for i in range(1, _Ne+1)]
        fo = dict()
        f = xr.open_dataset(f'processed_obs/{obstype}_{year}{month}.nc')
        for varname in varnames:
            fo[varname] = f[varname].to_numpy().ravel()
            fo[varname][fo[varname] > 1e10] = np.nan
        f.close()

        # remove masked grid points
        variable = dict()
        for varname in varnames:
            mask = self.mask | np.isnan(fo[varname])
            fo[varname] = fo[varname][~mask]
            variable[varname] = np.array([f[varname][0, 0].to_numpy().ravel()[~mask] for f in fens], dtype=float)

        # rank histogram
        self.getRH(variable, fo, varnames)

    def getRH(self, variable, fo, varnames):
        self.rh = dict()
        for varname in varnames:
            self.rh[varname] = None
            nTotalGp = len(fo[varname])
            ngridpoints = nTotalGp//self.nPE
            istart = self.rank*ngridpoints
            iend = istart + ngridpoints
            if self.rank == self.nPE - 1: iend = nTotalGp
            ngridpoints = iend - istart


            fo[varname] = fo[varname][istart:iend]
            variable[varname] = variable[varname][:, istart:iend] + \
                np.clip(np.random.normal(0, 0.3*fo[varname], size=ngridpoints), 0, None)
            variable[varname] = np.sort(np.vstack((variable[varname],
                                                   fo[varname]
                                                   )
                                                  ),
                                        axis=0, kind='heapsort'
                                        )
            length = variable[varname].shape[1]
            d = np.abs(variable[varname] - fo[varname]) < 1e-10
            d = d.ravel(order='F')
            rh = np.zeros(length, dtype=np.int32)
            mod_obs.mod_obs.getrank(_Ne, d, rh, length, length*(_Ne+1))
            if self.rank == 0:
                self.rh[varname] = np.empty(nTotalGp, dtype='i')
            # gather ranks from each processes
            sendcounts = np.array(self.comm.gather(length, 0))
            self.comm.Gatherv(sendbuf=rh, recvbuf=(self.rh[varname], sendcounts), root=0)

    def getMask(self):
        f = xr.open_dataset('../INPUTS/PHYSICS/DOM/eORCA_R1_zps_domcfg.nc')
        self.mask = f['top_level'][0].to_numpy().ravel() < 1e-3

    def plot(self, month):
        fig = plt.figure(1)
        fig.clf()
        ax = fig.add_subplot(121)
        ax.set_title('CHD')
        ax.hist(self.rh['CHD'], bins=_Ne+1)
        ax = fig.add_subplot(122)
        ax.hist(self.rh['CHN'], bins=_Ne+1)
        ax.set_title('CHN')
        plt.savefig(f'hist_{month}.png', dpi=300)


if __name__ == '__main__':
    obstype = 'chlo'
    expname = 'exp-logmean-chlo'
    varnames = ['CHD', 'CHN']
    for ic in ['01', '04', '07', '10']:
        for month in range(int(ic), int(ic)+3):
            rh = RankHistogram(obstype, expname, '', '2015', str(month).zfill(2), ic, varnames=varnames)
            rank = MPI.COMM_WORLD.Get_rank()
            if rank == 0:
                rh.plot(str(month).zfill(2))
