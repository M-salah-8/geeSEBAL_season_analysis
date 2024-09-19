import os
import json
import datetime

from ecmwf_utils import download_CDS_data

with open(os.path.join(os.path.dirname(__file__), "ecmwf_data_download_parameters.json")) as f:
    parameters = json.load(f)
area = parameters["area"]
date = datetime.datetime.strptime(parameters["date"], "%Y-%m-%d")
download_path = parameters["download_path"]
download_temperature = parameters["download_temperature"]
download_dewpoint = parameters["download_dewpoint"]
download_pressure = parameters["download_pressure"]
download_wind_speed = parameters["download_wind_speed"]
download_clear_sky_solar_radiation = parameters["download_clear_sky_solar_radiation"]
download_solar_radiation = parameters["download_solar_radiation"]
overwrite = parameters["overwrite"]



fields = []
if download_temperature:
    fields.extend(['2m_temperature', 'zero_degree_level', '2m_dewpoint_temperature', 'surface_pressure'])
if download_dewpoint and '2m_dewpoint_temperature' not in fields:
    fields.append('2m_dewpoint_temperature')
if download_pressure and 'surface_pressure' not in fields:
    fields.append('surface_pressure')
if download_wind_speed:
    fields.extend(['100m_v_component_of_wind', '100m_u_component_of_wind'])
if download_clear_sky_solar_radiation:
    fields.append("surface_solar_radiation_downward_clear_sky")
if download_solar_radiation:
    fields.append('surface_solar_radiation_downwards')

download_CDS_data(date, fields, download_path, overwrite, area)