import rasterio
import rasterio.warp
import datetime
import os

import cdsapi

def download_era5_re_hourly(image_dr, date, gdf):
    era5_file = os.path.join(image_dr, 'era5-hourly.grib')
    if os.path.exists(era5_file):
        print("-----------era5 data already exists-----------")
    else:
        print("-----------downloading era5 data-----------")
        previous_day = date - datetime.timedelta(days=1)
        year = [date.year.__str__()]
        months = [date.month.__str__().zfill(2)]
        days = [previous_day.day.__str__().zfill(2), date.day.__str__().zfill(2)]
        gdf = gdf.to_crs("epsg:4326")
        bbox = gdf.total_bounds
        area = int(bbox[3]+1), int(bbox[0]), int(bbox[1]), int(bbox[2]+1)

        dataset = "reanalysis-era5-single-levels"
        request = {
            'product_type': ['reanalysis'],
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
            'data_format': 'grib',
            'download_format': 'unarchived',
            'variable': [
                    '10m_u_component_of_wind', '10m_v_component_of_wind', '2m_dewpoint_temperature',
                    '2m_temperature', 'surface_solar_radiation_downwards',
                ],
            'area': area
        }
        client = cdsapi.Client()
        client.retrieve(dataset, request, era5_file)
        print("-----------downloaded-----------")