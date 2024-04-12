import cdsapi

c = cdsapi.Client()

c.retrieve(
    'reanalysis-era5-single-levels',
    {
        'product_type': 'ensemble_members',
        'format': 'grib',
        'variable': [
            '10m_u_component_of_wind', '10m_v_component_of_wind', '2m_dewpoint_temperature',
            '2m_temperature', 'mean_sea_level_pressure', 'mean_snowfall_rate',
            'mean_surface_downward_long_wave_radiation_flux', 'mean_surface_downward_short_wave_radiation_flux', 'mean_total_precipitation_rate',
            'surface_pressure',
        ],
        'day': [
            '01', '02', '03',
            '04', '05', '06',
            '07', '08', '09',
            '10', '11', '12',
            '13', '14', '15',
            '16', '17', '18',
            '19', '20', '21',
            '22', '23', '24',
            '25', '26', '27',
            '28', '29', '30',
            '31',
        ],
        'time': [
            '00:00', '03:00', '06:00',
            '09:00', '12:00', '15:00',
            '18:00', '21:00',
        ],
        'year': '2017',
        'month': [
            '01', # '02', '03',
          #   '04', '05', '06',
          #   '07', '08', '09',
          #   '10', '11', '12',
        ],
    },
    'download2017.grib')
