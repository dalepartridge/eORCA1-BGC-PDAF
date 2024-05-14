# Define experiment name
EXP = 'PC'
# Define file paths
DOMAIN_PATH = '/work/n01/n01/ymchen/eORCA1-BGC-PDAF/INPUTS/PHYSICS/DOM/eORCA_R1_zps_domcfg.nc'
OBS_PATH = '/work/n01/n01/ymchen/obs'
OUTPUT_PATH = f'/work/n01/n01/ymchen/eORCA1-BGC-PDAF/OUTPUTS/{EXP}'

# Define grid dimensions
ne, nx, ny = 30, 362, 332
__monthsEns = ['04', '08', '12']

# Define variable information
import operators
varinfo = {
    'chlorophyll': {
        'varname': ['CHN', 'CHD'],
        'op': operators.transformSum,
        'file_format': 'eORCA1_1d_{filerange}_ptrc1_T_{year}{month}-{year}{month}.nc'
    },
    'nitrogen': {
        'varname': ['PHN', 'PHD'],
        'op': operators.transformSum,
        'file_format': 'eORCA1_1d_{filerange}_ptrc1_T_{year}{month}-{year}{month}.nc'
    },
    'SST': {
        'varname': ['toce_con'],
        'op': operators.identity,
        'file_format': 'eORCA1_1d_{filerange}_grid_T_{year}{month}-{year}{month}.nc'
    },
    'SSS': {
        'varname': ['soce_abs'],
        'op': operators.identity,
        'file_format': 'eORCA1_1d_{filerange}_grid_T_{year}{month}-{year}{month}.nc'
    }
}
