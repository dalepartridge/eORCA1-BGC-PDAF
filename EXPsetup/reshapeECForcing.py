import numpy as np
from netCDF4 import Dataset, Variable


def reshapeNC(input_filename, output_filename, varname):
    # Open the original netCDF file in read mode
    print (f'open {input_filename}...')
    original_file = Dataset(input_filename, 'r')

    t0 = original_file['time'][0]
    t1 = original_file['time'][-1]
    t = np.arange(t0, t1+4*3, 3)
    nt = len(t)
    print (nt)

    # Create a new netCDF file in write mode
    print (f'open {output_filename}...')
    new_file = Dataset(output_filename, 'w')

    # Copy global attributes from the original file to the new file
    for attr_name in original_file.ncattrs():
        new_file.setncattr(attr_name, original_file.getncattr(attr_name))

    # Define the new dimensions (time1, y, x)
    time = new_file.createDimension('time', nt)
    y = new_file.createDimension('Y', original_file.dimensions['Y'].size)
    x = new_file.createDimension('X', original_file.dimensions['X'].size)

    # Iterate through all variables in the original file
    for var_name, var in original_file.variables.items():
        if var_name == 'step': continue
        if var_name == varname:
            print (f'writing {var_name}...')
            new_var = new_file.createVariable(var_name, var.dtype, ('time', 'Y', 'X'), fill_value=var.getncattr("_FillValue"))
            new_var[:] = var[:].reshape(nt, original_file.dimensions['Y'].size, original_file.dimensions['X'].size)
            print (np.max(np.abs(new_var[:].ravel() - var[:].ravel())))
            
            # Copy variable attributes
            for attr_name in var.ncattrs():
                if attr_name == '_FillValue':continue
                new_var.setncattr(attr_name, var.getncattr(attr_name))
        else:
            # Copy all other variables to the new file
            if hasattr(var, '_FillValue'):
                new_var = new_file.createVariable(var_name, var.dtype, var.dimensions, fill_value=var.getncattr('_FillValue'))
            else:
                new_var = new_file.createVariable(var_name, var.dtype, var.dimensions)
            new_var[:] = t[:] if var_name == 'time' else var[:]
            
            # Copy variable attributes
            for attr_name in var.ncattrs():
                if attr_name == '_FillValue':continue
                new_var.setncattr(attr_name, var.getncattr(attr_name))

    # Close both files
    original_file.close()
    new_file.close()


def DoReshape(i):
    reshapeNC(f'INPUT/ERA5_ens_{i}/ERA5_msr_y2015.nc',
              f'INPUT/ERA5_ens_{i}/ERA5_msr_y2015_reshape.nc',
              'msr')
    reshapeNC(f'INPUT/ERA5_ens_{i}/ERA5_msdwlwrf_y2015.nc',
              f'INPUT/ERA5_ens_{i}/ERA5_msdwlwrf_y2015_reshape.nc',
              'msdwlwrf')
    reshapeNC(f'INPUT/ERA5_ens_{i}/ERA5_msdwswrf_y2015.nc',
              f'INPUT/ERA5_ens_{i}/ERA5_msdwswrf_y2015_reshape.nc',
              'msdwswrf')
    reshapeNC(f'INPUT/ERA5_ens_{i}/ERA5_mtpr_y2015.nc',
              f'INPUT/ERA5_ens_{i}/ERA5_mtpr_y2015_reshape.nc',
              'mtpr')

if __name__ == '__main__':
    import multiprocessing as mp
    process = []
    for i in range(10):
        p = mp.Process(target=DoReshape, args=(i,))
        p.start()
        process.append(p)

    for p in process:
        p.join()

