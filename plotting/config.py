"""Configurations of the experiment diagnostics.
"""
import typing
import operators
# Define experiment name
exp: str
# Define file paths
DOMAIN_PATH: str = '/gws/nopw/j04/nceo_generic/users/yumengch/NEMO_PDAF/eORCA_R1_zps_domcfg.nc'
OBS_PATH: str = '/gws/nopw/j04/nceo_generic/users/yumengch/NEMO_PDAF/obs'
output_path_format: str = '/gws/nopw/j04/nceo_generic/users/yumengch/NEMO_PDAF/{exp}'
output_path: str
do_nanclim_obs: bool = True
n_process: int = 16

exp_labels: dict[str, str] = {
    'chlo': 'Daily Chl', 'chlo-monthly': 'Monthly Chl',
    'chlo-monthly-update': 'Monthly Chl+',
    'PC': 'Monthly C', 'chlo-pc': 'Monthly C & Chl',
    'PC-update': 'Monthly C+',
    'obs': 'Obs',
    'free': 'Freerun'}

# Define grid dimensions
ne: int = 30
nx: int = 362
ny: int = 332
monthsEns: list[str] = ['04', '08', '12']

# Define variable information
varinfo: dict[str, dict[str,
                        typing.Union[str, typing.Union
                                     [list[str],
                                      typing.Callable
                                      ]
                                     ]
                        ]
              ] = {}
varinfo['chlorophyll'] = {
    'varname': ['CHN', 'CHD'],
    'op': operators.transform_sum,
    'file_format': 'eORCA1_1d_{filerange}_ptrc1_T_{year}{month}-{year}{month}.nc'
}
varinfo['nitrogen'] = {
    'varname': ['PHN', 'PHD'],
    'op': operators.transform_sum,
    'file_format': 'eORCA1_1d_{filerange}_ptrc1_T_{year}{month}-{year}{month}.nc'
}
varinfo['chlo-monthly'] = {
    'varname': ['CHN', 'CHD'],
    'op': operators.transform_sum,
    'file_format': 'eORCA1_1m_{filerange}_ptrc_T_{year}{month}-{year}{month}.nc'
}
varinfo['nitrogen-monthly'] = {
    'varname': ['PHN', 'PHD'],
    'op': operators.transform_sum,
    'file_format': 'eORCA1_1m_{filerange}_ptrc_T_{year}{month}-{year}{month}.nc'
}
varinfo['SST'] = {
    'varname': ['toce_con',],
    'op': operators.identity,
    'file_format': 'eORCA1_1d_{filerange}_grid_T_{year}{month}-{year}{month}.nc'
}
varinfo['SSS'] = {
    'varname': ['soce_abs',],
    'op': operators.identity,
    'file_format': 'eORCA1_1d_{filerange}_grid_T_{year}{month}-{year}{month}.nc'
}
varinfo['CHN'] = {
    'varname': ['CHN', ],
    'op': operators.identity,
    'file_format': 'eORCA1_1m_{filerange}_ptrc_T_{year}{month}-{year}{month}.nc'
}
varinfo['CHD'] = {
    'varname': ['CHD', ],
    'op': operators.identity,
    'file_format': 'eORCA1_1m_{filerange}_ptrc_T_{year}{month}-{year}{month}.nc'
}
varinfo['PHN'] = {
    'varname': ['PHN', ],
    'op': operators.identity,
    'file_format': 'eORCA1_1m_{filerange}_ptrc_T_{year}{month}-{year}{month}.nc'
}
varinfo['PHD'] = {
    'varname': ['PHD',],
    'op': operators.identity,
    'file_format': 'eORCA1_1m_{filerange}_ptrc_T_{year}{month}-{year}{month}.nc'
}
varinfo['PDS'] = {
    'varname': ['PDS',],
    'op': operators.identity,
    'file_format': 'eORCA1_1m_{filerange}_ptrc_T_{year}{month}-{year}{month}.nc'
}
varinfo['chlorophyll-log10'] = {
    'varname': ['CHN', 'CHD'],
    'op': operators.transform_sum_log10,
    'file_format': 'eORCA1_1d_{filerange}_ptrc1_T_{year}{month}-{year}{month}.nc'
}
varinfo['nitrogen-log10'] = {
    'varname': ['PHN', 'PHD'],
    'op': operators.transform_sum_log10,
    'file_format': 'eORCA1_1d_{filerange}_ptrc1_T_{year}{month}-{year}{month}.nc'
}
varinfo['temperature'] = {
    'varname': ['toce_con',],
    'op': operators.identity,
    'file_format': 'eORCA1_1d_{filerange}_grid_T_{year}{month}-{year}{month}.nc'
}
varinfo['OXY'] = {
    'varname': ['OXY',],
    'op': operators.identity,
    'file_format': 'eORCA1_1m_{filerange}_ptrc_T_{year}{month}-{year}{month}.nc'
}
varinfo['DIN'] = {
    'varname': ['DIN',],
    'op': operators.identity,
    'file_format': 'eORCA1_1m_{filerange}_ptrc_T_{year}{month}-{year}{month}.nc'
}
varinfo['DIC'] = {
    'varname': ['DiC',],
    'op': operators.identity,
    'file_format': 'eORCA1_1m_{filerange}_ptrc_T_{year}{month}-{year}{month}.nc'
}
varinfo['SIL'] = {
    'varname': ['SIL',],
    'op': operators.identity,
    'file_format': 'eORCA1_1m_{filerange}_ptrc_T_{year}{month}-{year}{month}.nc'
}
varinfo['FER'] = {
    'varname': ['FER',],
    'op': operators.identity,
    'file_format': 'eORCA1_1m_{filerange}_ptrc_T_{year}{month}-{year}{month}.nc'
}
varinfo['ALK'] = {
    'varname': ['ALK',],
    'op': operators.identity,
    'file_format': 'eORCA1_1m_{filerange}_ptrc_T_{year}{month}-{year}{month}.nc'
}
varinfo['DET'] = {
    'varname': ['DET',],
    'op': operators.identity,
    'file_format': 'eORCA1_1m_{filerange}_ptrc_T_{year}{month}-{year}{month}.nc'
}
varinfo['ZME'] = {
    'varname': ['ZME',],
    'op': operators.identity,
    'file_format': 'eORCA1_1m_{filerange}_ptrc_T_{year}{month}-{year}{month}.nc'
}
varinfo['ZMI'] = {
    'varname': ['ZMI',],
    'op': operators.identity,
    'file_format': 'eORCA1_1m_{filerange}_ptrc_T_{year}{month}-{year}{month}.nc'
}
varinfo['TCO2'] = {
    'varname': ['TCO2',],
    'op': operators.identity,
    'file_format': 'eORCA1_1m_{filerange}_diad_T_{year}{month}-{year}{month}.nc'
}
varinfo['CO2FLUX'] = {
    'varname': ['CO2FLUX',],
    'op': operators.identity,
    'file_format': 'eORCA1_1m_{filerange}_diad_T_{year}{month}-{year}{month}.nc'
}
varinfo['ATM_PCO2'] = {
    'varname': ['ATM_PCO2',],
    'op': operators.identity,
    'file_format': 'eORCA1_1m_{filerange}_diad_T_{year}{month}-{year}{month}.nc'
}
varinfo['OCN_PCO2'] = {
    'varname': ['OCN_PCO2',],
    'op': operators.identity,
    'file_format': 'eORCA1_1m_{filerange}_diad_T_{year}{month}-{year}{month}.nc'
}
varinfo['O2FLUX'] = {
    'varname': ['O2FLUX',],
    'op': operators.identity,
    'file_format': 'eORCA1_1m_{filerange}_diad_T_{year}{month}-{year}{month}.nc'
}
varinfo['PRN'] = {
    'varname': ['PRN',],
    'op': operators.identity,
    'file_format': 'eORCA1_1m_{filerange}_diad_T_{year}{month}-{year}{month}.nc'
}
varinfo['PRD'] = {
    'varname': ['PRD',],
    'op': operators.identity,
    'file_format': 'eORCA1_1m_{filerange}_diad_T_{year}{month}-{year}{month}.nc'
}
varinfo['N_PROD'] = {
    'varname': ['N_PROD',],
    'op': operators.identity,
    'file_format': 'eORCA1_1m_{filerange}_diad_T_{year}{month}-{year}{month}.nc'
}
varinfo['N_CONS'] = {
    'varname': ['N_CONS',],
    'op': operators.identity,
    'file_format': 'eORCA1_1m_{filerange}_diad_T_{year}{month}-{year}{month}.nc'
}
varinfo['C_CONS'] = {
    'varname': ['C_CONS',],
    'op': operators.identity,
    'file_format': 'eORCA1_1m_{filerange}_diad_T_{year}{month}-{year}{month}.nc'
}
varinfo['O2_PROD'] = {
    'varname': ['O2_PROD',],
    'op': operators.identity,
    'file_format': 'eORCA1_1m_{filerange}_diad_T_{year}{month}-{year}{month}.nc'
}
varinfo['O2_CONS'] = {
    'varname': ['O2_CONS',],
    'op': operators.identity,
    'file_format': 'eORCA1_1m_{filerange}_diad_T_{year}{month}-{year}{month}.nc'
}
varinfo['MPN'] = {
    'varname': ['MPN',],
    'op': operators.identity,
    'file_format': 'eORCA1_1m_{filerange}_diad_T_{year}{month}-{year}{month}.nc'
}
varinfo['MPD'] = {
    'varname': ['MPD',],
    'op': operators.identity,
    'file_format': 'eORCA1_1m_{filerange}_diad_T_{year}{month}-{year}{month}.nc'
}
varinfo['GMIPn'] = {
    'varname': ['GMIPn',],
    'op': operators.identity,
    'file_format': 'eORCA1_1m_{filerange}_diad_T_{year}{month}-{year}{month}.nc'
}
varinfo['GMEPN'] = {
    'varname': ['GMEPN',],
    'op': operators.identity,
    'file_format': 'eORCA1_1m_{filerange}_diad_T_{year}{month}-{year}{month}.nc'
}
varinfo['GMEPD'] = {
    'varname': ['GMEPD',],
    'op': operators.identity,
    'file_format': 'eORCA1_1m_{filerange}_diad_T_{year}{month}-{year}{month}.nc'
}
varinfo['GMEZMI'] = {
    'varname': ['GMEZMI',],
    'op': operators.identity,
    'file_format': 'eORCA1_1m_{filerange}_diad_T_{year}{month}-{year}{month}.nc'
}
varinfo['PD_JLIM'] = {
    'varname': ['PD_JLIM',],
    'op': operators.identity,
    'file_format': 'eORCA1_1m_{filerange}_diad_T_{year}{month}-{year}{month}.nc'
}
varinfo['PN_JLIM'] = {
    'varname': ['PN_JLIM',],
    'op': operators.identity,
    'file_format': 'eORCA1_1m_{filerange}_diad_T_{year}{month}-{year}{month}.nc'
}
varinfo['PD_NLIM'] = {
    'varname': ['PD_NLIM',],
    'op': operators.identity,
    'file_format': 'eORCA1_1m_{filerange}_diad_T_{year}{month}-{year}{month}.nc'
}
varinfo['PD_FELIM'] = {
    'varname': ['PD_FELIM',],
    'op': operators.identity,
    'file_format': 'eORCA1_1m_{filerange}_diad_T_{year}{month}-{year}{month}.nc'
}
varinfo['PD_SILIM'] = {
    'varname': ['PD_SILIM',],
    'op': operators.identity,
    'file_format': 'eORCA1_1m_{filerange}_diad_T_{year}{month}-{year}{month}.nc'
}
varinfo['PDSILIM2'] = {
    'varname': ['PDSILIM2',],
    'op': operators.identity,
    'file_format': 'eORCA1_1m_{filerange}_diad_T_{year}{month}-{year}{month}.nc'
}
varinfo['PDSILIM2'] = {
    'varname': ['PDSILIM2',],
    'op': operators.identity,
    'file_format': 'eORCA1_1m_{filerange}_diad_T_{year}{month}-{year}{month}.nc'
}
varinfo['PN_NLIM'] = {
    'varname': ['PN_NLIM',],
    'op': operators.identity,
    'file_format': 'eORCA1_1m_{filerange}_diad_T_{year}{month}-{year}{month}.nc'
}
varinfo['PN_FELIM'] = {
    'varname': ['PN_FELIM',],
    'op': operators.identity,
    'file_format': 'eORCA1_1m_{filerange}_diad_T_{year}{month}-{year}{month}.nc'
}
varinfo['PN_LLOSS'] = {
    'varname': ['PN_LLOSS',],
    'op': operators.identity,
    'file_format': 'eORCA1_1m_{filerange}_diad_T_{year}{month}-{year}{month}.nc'
}
varinfo['PD_LLOSS'] = {
    'varname': ['PD_LLOSS',],
    'op': operators.identity,
    'file_format': 'eORCA1_1m_{filerange}_diad_T_{year}{month}-{year}{month}.nc'
}
varinfo['MED_XPAR'] = {
    'varname': ['MED_XPAR',],
    'op': operators.identity,
    'file_format': 'eORCA1_1m_{filerange}_diad_T_{year}{month}-{year}{month}.nc'
}
