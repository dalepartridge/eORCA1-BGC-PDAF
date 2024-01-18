import numpy as np
import netCDF4
import matplotlib.pyplot as plt
import matplotlib.gridspec


_colors = ['#0173b2', '#de8f05', '#029e73', '#d55e00', '#cc78bc', '#ca9161', '#fbafe4', '#949494']
def plotObsDistribution(obstype):
    fo = np.load(f'processed_obs/{obstype}_201506.npz')
    chlo = dict()
    for varname in fo:
        chlo[varname] = fo[varname][:].ravel()
        chlo[varname][chlo[varname] > 1e10] = np.nan
        chlo[varname][chlo[varname] < 1e-10] = np.nan
    fo.close()

    dataName = 'Chlorophyll' if obstype == 'chlo' else 'Carbon'

    gs = matplotlib.gridspec.GridSpec(1, 2)#, wspace=0.1, left=0.1, right=1., bottom=0.2, top=1.)
    fig = plt.figure(1)
    fig.clf()
    ax = fig.add_subplot(gs[0])
    for varname, color in zip(chlo, _colors):
        ax.hist(chlo[varname], bins=100, label=varname, alpha=0.5, color=color)
    ax.hist(sum(chlo.values()), bins=100, label='total', alpha=0.5, color='r')
    ax.legend()
    ax.set_title(dataName)

    ax = fig.add_subplot(gs[1])
    for varname, color in zip(chlo, _colors):
        ax.hist(np.log(chlo[varname]), bins=100, label=varname, alpha=0.5, color=color)
    ax.hist(np.log(sum(chlo.values())), bins=100, label='total', alpha=0.5, color='r')
    ax.set_title(f'log({dataName})')
    ax.legend()
    fig.savefig(f'{obstype}_dist.png', dpi=300)


def plotObsDistributionNC(obstype):
    chlo = dict()
    if obstype == 'chlo':
        fo = netCDF4.Dataset('/work/n01/n01/ymchen/obs/ESACCI-OC-L3S-CHLOR_A-MERGED-1M_MONTHLY_4km_GEO_PML_OCx-201506-fv5.0.nc', 'r')
        chlo_t = 10**(np.log10(fo['chlor_a'][:].ravel() + fo['chlor_a_log10_bias'][:].ravel()))
        nano = 1.057*(1-np.exp(-0.851*chlo_t))
        micro = chlo_t - nano

        chlo['CHD'] = micro*2/3.
        chlo['CHN'] = chlo_t - chlo['CHD']
    else:
        fo = netCDF4.Dataset('/work/n01/n01/ymchen/obs/BICEP_NCEO_PC_ESA-OC-L3S-MERGED-1M_MONTHLY_9km_mapped_201506-fv5.nc', 'r')
        micro = fo['C_microphyto'][:].ravel()
        chlo_t =  fo['C_microphyto'][:].ravel() + fo['C_nanophyto'][:].ravel() + fo['C_picophyto'][:].ravel()
        chlo['PHD'] = micro*2/3.
        chlo['PHN'] = chlo_t - chlo['PHD']
    fo.close()


    chlo = dict()
    for varname in chlo:
        chlo[varname][chlo[varname] > 1e10] = np.nan
        chlo[varname][chlo[varname] < 1e-10] = np.nan

    dataName = 'Chlorophyll' if obstype == 'chlo' else 'Carbon'

    gs = matplotlib.gridspec.GridSpec(1, 2)#, wspace=0.1, left=0.1, right=1., bottom=0.2, top=1.)
    fig = plt.figure(1)
    fig.clf()
    ax = fig.add_subplot(gs[0])
    for varname, color in zip(chlo, _colors):
        ax.hist(chlo[varname], bins=100, label=varname, alpha=0.5, color=color)
    ax.hist(chlo_t, bins=100, label='total', alpha=0.5, color='r')
    ax.legend()
    ax.set_title(dataName)

    ax = fig.add_subplot(gs[1])
    for varname, color in zip(chlo, _colors):
        ax.hist(np.log(chlo[varname]), bins=100, label=varname, alpha=0.5, color=color)
    ax.hist(np.log(chlo_t), bins=100, label='total', alpha=0.5, color='r')
    ax.set_title(f'log({dataName})')
    ax.legend()
    fig.savefig(f'{obstype}_dist_nc.png', dpi=300)

if __name__ == '__main__':
    plotObsDistribution('chlo')
    plotObsDistributionNC('chlo')
    plotObsDistribution('pc')
    plotObsDistributionNC('pc')