# -*- coding: utf-8 -*-
"""
Created on Tue Jul 17 14:47:50 2018

@author: rmgu
"""

import sys, os
sys.path.append(os.getcwd())
import datetime
import numpy as np
from osgeo import gdal
import rasterio
import netCDF4

from pyTSEB import meteo_utils as met

import gdal_utils as gu
from funcs import geo

# Acceleration of gravity (m s-2)
GRAVITY = 9.80665
# Blending height of 100 m
Z_BH = 100.0


def download_CDS_data(date, variables, target, overwrite=False, area=None):
    previous_day = date - datetime.timedelta(days=1)
    next_day = date + datetime.timedelta(days=1)
    print(previous_day, date, next_day)
    year = [date.year.__str__()]
    months = [date.month.__str__().zfill(2)]
    days = [previous_day.day.__str__().zfill(2), date.day.__str__().zfill(2), next_day.day.__str__().zfill(2)]

    import cdsapi

    dataset = "reanalysis-era5-single-levels"
    request = {
        'product_type': ['reanalysis'],
        'year': year,
        'month': months,
        'day': days,
        'time': ['00:00', '01:00', '02:00', '03:00', '04:00', '05:00', '06:00', '07:00', '08:00', '09:00', '10:00', '11:00', '12:00', '13:00', '14:00', '15:00', '16:00', '17:00', '18:00', '19:00', '20:00', '21:00', '22:00', '23:00'],
        'data_format': 'netcdf',
        'download_format': 'unarchived',
        'variable': variables,
        'area': area
    }
    client = cdsapi.Client()
    era5_file = os.path.join(target, 'era5.nc')
    client.retrieve(dataset, request, era5_file).download()

    print("------------------Downloaded-------------------")


def get_ECMWF_data(ecmwf_data_file, field, timedate_UTC, elev, time_zone):

    ncfile = netCDF4.Dataset(ecmwf_data_file, 'r')
    era5_temp = os.path.join(os.path.dirname(ecmwf_data_file),"temp.tif")
    # Find the location of bracketing dates
    time = ncfile.variables['valid_time']
    dates = netCDF4.num2date(time[:], time.units, time.calendar)
    beforeI, afterI, frac = _bracketing_dates(dates, timedate_UTC)
    if beforeI is None:
        ncfile.close()
        return None
    time = None
    ncfile.close()

    if field == "air_temperature":
        print("air_temperature")
        t2m, gt, proj = _getECMWFTempInterpData(ecmwf_data_file, "t2m", beforeI, afterI, frac)
        # Get geopotential height at which Ta is calculated
        z, gt, proj = _getECMWFTempInterpData(ecmwf_data_file, "deg0l", beforeI, afterI, frac)
        z /= GRAVITY
        d2m, gt, proj = _getECMWFTempInterpData(ecmwf_data_file, "d2m", beforeI, afterI, frac)
        ea = calc_vapour_pressure(d2m)
        sp, gt, proj = _getECMWFTempInterpData(ecmwf_data_file, "sp", beforeI, afterI, frac)
        p = calc_pressure_mb(sp)
        # Calcultate temperature at 0m datum height
        T_datum = calc_air_temperature_blending_height(t2m, ea, p, 0, z_ta=z+2.0)
        # Resample dataset and calculate actual blendingh height temperature based on input
        # elevation data
        ea = _ECMWFRespampleData(ea, gt, proj, elev, era5_temp)
        p = _ECMWFRespampleData(p, gt, proj, elev, era5_temp)
        T_datum = _ECMWFRespampleData(T_datum, gt, proj, elev, era5_temp)
        elev_data = gu.raster_data(elev)
        data = calc_air_temperature_blending_height(T_datum, ea, p, elev_data+Z_BH, z_ta=0)

    elif field == "vapour_pressure":
        print("vp")
        d2m, gt, proj = _getECMWFTempInterpData(ecmwf_data_file, "d2m", beforeI, afterI, frac)
        data = calc_vapour_pressure(d2m)
        data = _ECMWFRespampleData(data, gt, proj, elev, era5_temp)

    elif field == "wind_speed":
        print("ws")
        u100, gt, proj = _getECMWFTempInterpData(ecmwf_data_file, "u100", beforeI, afterI, frac)
        v100, gt, proj = _getECMWFTempInterpData(ecmwf_data_file, "v100", beforeI, afterI, frac)
        # Combine the two components of wind speed and calculate speed at blending height
        ws100 = calc_wind_speed(u100, v100)
        data = _ECMWFRespampleData(ws100, gt, proj, elev, era5_temp)

    elif field == "air_pressure":
        print("ap")
        sp, gt, proj = _getECMWFTempInterpData(ecmwf_data_file, "sp", beforeI, afterI, frac)
        # Convert pressure from pascals to mb
        data = calc_pressure_mb(sp)
        data = _ECMWFRespampleData(data, gt, proj, elev, era5_temp)

    elif field == "clear_sky_solar_radiation":
        print("cssr")
        ssrdc, gt, proj = _getECMWFTempInterpData(ecmwf_data_file, "ssrdc", beforeI, afterI, frac)
        # Convert from Jules to Watts
        data = ssrdc / 3600.0
        data = _ECMWFRespampleData(data, gt, proj, elev, era5_temp)

    elif field == "average_daily_solar_irradiance":
        print("adsi")
        # Find midnight in local time and convert to UTC time
        date_local = (timedate_UTC + datetime.timedelta(hours=time_zone)).date()
        midnight_local = datetime.datetime.combine(date_local, datetime.time())
        midnight_UTC = midnight_local - datetime.timedelta(hours=time_zone)
        # Interpolate solar irradiance over 24 hour period starting at midnight local time
        data, gt, proj = _getECMWFIntegratedData(ecmwf_data_file, "ssrd", midnight_UTC,
                                                 time_window=24)
        data = _ECMWFRespampleData(data, gt, proj, elev, era5_temp)
    else:
        raise RuntimeError("Unknown field: %s!" % field)

    return data


def calc_air_temperature_blending_height(ta, ea, p, z_bh, z_ta=2.0):
    if type(ta) is np.ndarray:
        ta = ta.astype(np.float32)
    if type(ea) is np.ndarray:
        ea = ea.astype(np.float32)
    if type(p) is np.ndarray:
        p = p.astype(np.float32)
    if type(z_bh) is np.ndarray:
        z_bh = z_bh.astype(np.float32)
    if type(z_ta) is np.ndarray:
        z_ta = z_ta.astype(np.float32)
    lapse_rate = met.calc_lapse_rate_moist(ta, ea, p)
    ta_bh = ta - lapse_rate * (z_bh - z_ta)
    return ta_bh


def calc_vapour_pressure(td):
    # output in mb
    td = td - 273.15
    e = 6.11 * np.power(10, (7.5 * td)/(237.3 + td))
    return e


def calc_wind_speed(u, v):
    ws = (u**2 + v**2)**0.5
    ws = np.maximum(ws, 1.0)
    return ws


def calc_pressure_mb(sp):
    # Convert from pascals to mb
    sp_mb = sp/100.0
    return sp_mb


def calc_tcwv_cm(tcwv):
    # Conert from from kg/m**2 to g/cm**2
    return tcwv/10.0


def _getECMWFTempInterpData(ncfile, var_name, before_I, after_I, frac):
    ds = gdal.Open('NETCDF:"'+ncfile+'":'+var_name)
    if ds is None:
        raise RuntimeError("Variable %s does not exist in file %s." % (var_name, ncfile))

    # Get some metadata
    scale = ds.GetRasterBand(before_I+1).GetScale() if ds.GetRasterBand(before_I+1).GetScale() is not None else 1.0
    offset = ds.GetRasterBand(before_I+1).GetOffset() if ds.GetRasterBand(before_I+1).GetOffset() is not None else 0.0
    no_data_value = ds.GetRasterBand(before_I+1).GetNoDataValue()
    gt = ds.GetGeoTransform()
    sr = rasterio.crs.CRS.from_epsg(4326)
    proj = sr.to_wkt()

    # Read the right time layers
    try:
        data_before = ds.GetRasterBand(before_I+1).ReadAsArray()
        data_before = (data_before.astype(float) * scale) + offset
        data_before[data_before == no_data_value] = np.nan
        data_after = ds.GetRasterBand(after_I+1).ReadAsArray()
        data_after = (data_after.astype(float) * scale) + offset
        data_after[data_after == no_data_value] = np.nan
    except AttributeError:
        ds = None
        raise RuntimeError("ECMWF file does not contain data for the requested date.")

    # Perform temporal interpolation
    data = data_before*frac + data_after*(1.0-frac)

    return data, gt, proj


def _ECMWFRespampleData(data, gt, proj, template_file, era5_temp):
    # Subset and reproject to the template file extent and projection
    geo_in = geo()
    array, _ = geo_in.coregister(template_file, era5_temp, data, dtype=np.float32, resampling='cubic_spline')
    return array


def _getECMWFIntegratedData(ncfile, var_name, date_time, time_window=24,):

    # Open the netcdf time dataset
    fid = netCDF4.Dataset(ncfile, 'r')
    time = fid.variables['valid_time']
    dates = netCDF4.num2date(time[:], time.units, time.calendar)
    del fid

    # Get the time right before date_time, to use it as integrated baseline
    date_0, _, _ = _bracketing_dates(dates, date_time)
    # Get the time right before the temporal witndow set
    date_1, _, _ = _bracketing_dates(dates, date_time + datetime.timedelta(hours=time_window))

    ds = gdal.Open('NETCDF:"'+ncfile+'":'+var_name)
    if ds is None:
        raise RuntimeError("Variable %s does not exist in file %s." % (var_name, ncfile))
    # Get some metadata
    scale = ds.GetRasterBand(date_0+1).GetScale() if ds.GetRasterBand(date_0+1).GetScale() is not None else 1.0
    offset = ds.GetRasterBand(date_0+1).GetOffset() if ds.GetRasterBand(date_0+1).GetOffset() is not None else 0.0
    no_data_value = ds.GetRasterBand(date_0+1).GetNoDataValue()
    # Report geolocation of the top-left pixel of rectangle
    gt = ds.GetGeoTransform()
    sr = rasterio.crs.CRS.from_epsg(4326)
    proj = sr.to_wkt()

    # Forecasts of ERA5 the accumulations are since the previous post processing
    # (archiving)
    data_ref = 0

    # Initialize output variable
    cummulated_value = 0.

    try:
        for date_i in range(date_0+1, date_1+1):
            # Read the right time layers
            data = ds.GetRasterBand(date_i+1).ReadAsArray()
            data = (data.astype(float) * scale) + offset
            data[data == no_data_value] = 0
            # The time step value is the difference between  the actual timestep value and the
            # previous value
            cummulated_value += (data - data_ref)
    except AttributeError:
        ds = None
        raise RuntimeError("ECMWF file does not contain data for the requested date")

    # Convert to average W m^-2
    cummulated_value = cummulated_value / (time_window * 3600.)

    return cummulated_value, gt, proj


def _bracketing_dates(date_list, target_date):
    date_list = list(date_list)
    try:
        before = max([x for x in date_list if (target_date - x).total_seconds() >= 0])
        after = min([x for x in date_list if (target_date - x).total_seconds() <= 0])
    except ValueError:
        return None, None, np.nan
    if before == after:
        frac = 1
    else:
        frac = float((after - target_date).total_seconds())/float((after-before).total_seconds())
    return date_list.index(before), date_list.index(after), frac