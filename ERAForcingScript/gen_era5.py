import os, sys, glob
import numpy  as np
import datetime
import xarray as xr

class era5(object):
    """
    Generate ERA5 atmospheric forcing for regional NEMO config
    Loosly based on code by Nico.

    Defaults parameters are for AMM15
    """

    def __init__(self, hourly):
        self.year_init = 2016                      ## First year to process
        self.year_end  = 2016                     ## Last one [included]
        self.east      = 360              ## East Border
        self.west      = 0              ## West Border
        self.north     = 90             ## North Border
        self.south     = -90             ## South Border
        # ROOT PATH OF ERA5 DATA
        # self.path_ERA5 = '/work/n01/n01/dapa/NCEO/eORCA1-ERA5/download_2.grib' 
        self.path_ERA5 = '/work/n01/n01/ymchen/Forcing2016/download.grib' 
        # NEMO FORCING
        self.clean        = False            ## Clean extraction (longest bit)
        self.chunks={'time':50}

        if hourly == 3:
            self.sph_ON       = True             ## Compute specific humidity or not
            self.var_path = { 
                   "10m_u_component_of_wind" : "u10",
                   "10m_v_component_of_wind" : "v10",
                   "2m_temperature"          : "t2m", 
                   "mean_sea_level_pressure" : "msl", }

            if self.sph_ON :
               self.var_path[ "surface_pressure"  ] = 'sp'
               self.var_path[ "2m_dewpoint_temperature" ] = 'd2m'
        else:
            self.var_path = {
                   "mean_snowfall_rate"      : "msr" ,
                   "mean_surface_downward_long_wave_radiation_flux"  : "msdwlwrf",
                   "mean_surface_downward_short_wave_radiation_flux" : "msdwswrf",
                   "mean_total_precipitation_rate" : "mtpr" }

    def timeit(func):
        """ decorator for timing a function """ 

        def inner():
            t0 = datetime.datetime.now()
            func()
            t1 = datetime.datetime.now()
            print ('time elapsed = ', t1-t0)
            
        return inner
    
    def add_global_attrs(self, ds):
        """ set global attributes for netcdf """
    
        fmt = "%Y-%m-%d %H:%M:%S"
        ds.attrs['Created'] = datetime.datetime.now().strftime(fmt)
        ds.attrs['Description'] = 'ERA5 Atmospheric conditions for eORCA1 NEMO'
    
        return ds
    
    def interpolate_by_year(self, ds, nameVar):
        """
        Loop over each extracted year interpolating to the half
        time-step, saving each year.
        """
    
        # output name
        fout = self.path_FORCING + '/ERA5_' + nameVar + '_y2016.nc'

        if self.clean : os.system( "rm {0}".format( fout ) )
        if not os.path.exists( fout ) :

            # interpolate to half time-step
            Time = ds.time.values
            dt = (Time[1] - Time[0]) / 2
            half_time = (ds.time + dt).values
            ds = ds.interp(time=half_time)

            # format indexes and coords
            ds = self.format_nc(ds, nameVar)
            ds.time.encoding['units'] = 'hours since 2000-01-01 00:00:00'

            # save with encoding
            ds.to_netcdf(fout,unlimited_dims='time')

    def format_nc(self, da, nameVar):

        # mesh lat and lon
        mlon, mlat = np.meshgrid(da.longitude, da.latitude)
        lon_attrs={'long_name':'lon','units':'degrees_east'}
        lat_attrs={'long_name':'lat', 'units':'degrees_north'}
        mlon = xr.DataArray(mlon, dims=['Y','X'], attrs=lon_attrs)
        mlat = xr.DataArray(mlat, dims=['Y','X'], attrs=lat_attrs)
      
        # assign X/Y as indexes
        da = da.drop(['longitude','latitude'])
        da = da.rename({'longitude':'X','latitude':'Y'})
        da = da.assign_coords({'lon':mlon,'lat':mlat})
        
        # file information
        self.add_global_attrs(da)
 
        return da

    def process_all(self, ds, ens):
        # self.path_FORCING = '/work/n01/n01/dapa/NCEO/eORCA1-ERA5/ERA5_ens_'+str(ens)
        self.path_FORCING = '/work/n01/n01/ymchen/Forcing2016/ERA5_ens_'+str(ens)
        os.system("mkdir {0}".format(self.path_FORCING ) )
        
        ## Loop over each variable
        for dirVar, nameVar in self.var_path.items() :
        
            print ("================== {0} - {1} ==================".format(
                    dirVar, nameVar ))
            self.interpolate_by_year(ds[nameVar],nameVar)
        
        ##---------- PROCESS SPECIFIC HUMIDITY ----------------------     
        ## Compute Specific Humidity according to ECMWF documentation
        
        if self.sph_ON : 
        
            for iY in range(self.year_init, self.year_end+1) :
        
                # read
                d2m_path = self.path_FORCING + '/ERA5_d2m_y'\
                           + str(iY) + '.nc'
                sp_path  = self.path_FORCING + '/ERA5_sp_y'\
                           + str(iY) + '.nc'
                d2m = xr.open_dataarray(d2m_path, chunks=self.chunks)
                sp  = xr.open_dataarray(sp_path,  chunks=self.chunks) 
        
                # calculate sph
                esat = 611.21 * np.exp( 17.502 * (d2m-273.16) / (d2m-32.19) )
                dyrvap = 287.0597 / 461.5250
                sph = dyrvap * esat / ( sp - (1-dyrvap) * esat)
                sph.attrs = {'units':'1', 'standard_name':'specific humidity'}
         
                # save
                fout = self.path_FORCING + '/ERA5_SPH_y' + str(iY) + '.nc'
                sph.to_netcdf(fout)

if __name__ == '__main__':
    import cfgrib
    era = era5(3)
    ds = cfgrib.open_datasets(era.path_ERA5, chunks=era.chunks)[0]
    for ens in range(10):
        print("############# ENSEMBLE NUMBER {0} #########################".format(ens))
        era.process_all(ds.isel(number=ens),ens)
        for year in range(era.year_init, era.year_end+1):
            os.system(f'ncrename -v __xarray_dataarray_variable__,SPH ERA5_ens_${ens}/ERA5_SPH_y{year}.nc')
            os.system(f'ncks --mk_rec_dmn time ERA5_ens_${ens}/ERA5_SPH_y{year}.nc -o ERA5_ens_${ens}/ERA5_SPH_y{year}-unlimited.nc')
            os.system(f'mv ERA5_ens_${ens}/ERA5_SPH_y{year}-unlimited.nc ERA5_ens_${ens}/ERA5_SPH_y{year}.nc')

    era = era5(12)
    ds = cfgrib.open_datasets(era.path_ERA5, chunks=era.chunks)[1].isel(step=1)
    ds['time'] = ds.time + np.timedelta64(6,'h')
    for ens in range(10):
        print("############# ENSEMBLE NUMBER {0} #########################".format(ens))
        era.process_all(ds.isel(number=ens),ens)

