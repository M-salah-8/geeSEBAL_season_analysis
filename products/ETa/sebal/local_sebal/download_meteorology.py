import rasterio
import rasterio.warp
import datetime
import os

import cdsapi

def download_era5_re_hourly(tif_dr, date):
    previous_day = date - datetime.timedelta(days=1)
    year = date.year.__str__()
    months = [date.month.__str__().zfill(2)]
    days = [previous_day.day.__str__().zfill(2), date.day.__str__().zfill(2)]
    
    image_dr = os.path.split(tif_dr)[0]
    image_dr = os.path.split(image_dr)[0]
    src = rasterio.open(tif_dr)
    bbox = rasterio.warp.transform_bounds(src.crs, 'EPSG:4326', *src.bounds)
    bounds = int(bbox[3]+1), int(bbox[0]), int(bbox[1]), int(bbox[2]+1)

    c = cdsapi.Client()

    c.retrieve(
        'reanalysis-era5-single-levels',
        {
            'product_type': 'reanalysis',
            'format': 'grib',
            'variable': [
                '10m_u_component_of_wind', '10m_v_component_of_wind', '2m_dewpoint_temperature',
                '2m_temperature', 'surface_solar_radiation_downwards',
            ],
            'year': year,
            'month': months,
            'day': days,
            'time': [
                '00:00', '01:00', '02:00',
                '03:00', '04:00', '05:00',
                '06:00', '07:00', '08:00',
                '09:00', '10:00', '11:00',
                '12:00', '13:00', '14:00',
                '15:00', '16:00', '17:00',
                '18:00', '19:00', '20:00',
                '21:00', '22:00', '23:00',
            ],
            'area': bounds,
        },
        os.path.join(image_dr, 'era5-hourly.grib'))
