import tempfile
import os
import json
import datetime

import ecmwf_utils as eu
## snappy_utils should be imported last, as it modifies the system path
import snappy_utils as su

with open(os.path.join(os.path.dirname(__file__), "ecmwf_data_preparation_parameters.json")) as f:
    parameters = json.load(f)

elevation_map = parameters["elevation_map"]
elevation_band = parameters["elevation_band"]
ecmwf_data_file = parameters["ecmwf_data_file"]
date_time_utc = datetime.datetime.strptime(parameters["date_time_utc"], "%Y-%m-%d %H:%M")
time_zone = parameters["time_zone"]
prepare_temperature = parameters["prepare_temperature"]
prepare_vapour_pressure = parameters["prepare_vapour_pressure"]
prepare_air_pressure = parameters["prepare_air_pressure"]
prepare_wind_speed = parameters["prepare_wind_speed"]
prepare_clear_sky_solar_radiation = parameters["prepare_clear_sky_solar_radiation"]
prepare_daily_solar_irradiance = parameters["prepare_daily_solar_irradiance"]
output_file = parameters["output_file"]

# Save elevation to GeoTIFF becasue it will need to be read by GDAL later
temp_file = tempfile.NamedTemporaryFile(suffix=".tif", delete=False)
temp_elev_path = temp_file.name
temp_file.close()
su.copy_bands_to_file(elevation_map, temp_elev_path, [elevation_band])

# Calculate required meteorological parameters
bands = []
if prepare_temperature:
    data = eu.get_ECMWF_data(ecmwf_data_file, 'air_temperature', date_time_utc, temp_elev_path,
                                time_zone)
    bands.append({'band_data': data, 'band_name': 'air_temperature', 'description':
                    'Air temperature at 100 m above surface(K)'})
if prepare_vapour_pressure:
    data = eu.get_ECMWF_data(ecmwf_data_file, 'vapour_pressure', date_time_utc, temp_elev_path,
                                time_zone)
    bands.append({'band_data': data, 'band_name': 'vapour_pressure', 'description':
                    'Surface vapour pressure (mb)'})
if prepare_air_pressure:
    data = eu.get_ECMWF_data(ecmwf_data_file, 'air_pressure', date_time_utc, temp_elev_path,
                                time_zone)
    bands.append({'band_data': data, 'band_name': 'air_pressure', 'description':
                    'Surface air pressure (mb)'})
if prepare_wind_speed:
    data = eu.get_ECMWF_data(ecmwf_data_file, 'wind_speed', date_time_utc, temp_elev_path,
                                time_zone)
    bands.append({'band_data': data, 'band_name': 'wind_speed', 'description':
                    'Wind speed at 100 m above surface (m/s)'})
if prepare_clear_sky_solar_radiation:
    data = eu.get_ECMWF_data(ecmwf_data_file, 'clear_sky_solar_radiation', date_time_utc,
                                temp_elev_path, time_zone)
    bands.append({'band_data': data, 'band_name': 'clear_sky_solar_radiation', 'description':
                    'Instantenous clear sky surface solar irradiance (W/m^2)'})
if prepare_daily_solar_irradiance:
    data = eu.get_ECMWF_data(ecmwf_data_file, 'average_daily_solar_irradiance', date_time_utc,
                                temp_elev_path, time_zone)
    bands.append({'band_data': data, 'band_name': 'average_daily_solar_irradiance',
                    'description': 'Average daily solar irradiance (W/m^2)'})

# Save the output file
geo_coding = su.read_snappy_product(elevation_map, elevation_band)[1]
su.write_snappy_product(output_file, bands, 'ecmwfData', geo_coding)