#----------------------------------------------------------------------------------------#
#---------------------------------------//GEESEBAL//-------------------------------------#
#GEESEBAL - GOOGLE EARTH ENGINE APP FOR SURFACE ENERGY BALANCE ALGORITHM FOR LAND (SEBAL)
#CREATE BY: LEONARDO LAIPELT, RAFAEL KAYSER, ANDERSON RUHOFF AND AYAN FLEISCHMANN
#PROJECT - ET BRASIL https://etbrasil.org/
#LAB - HIDROLOGIA DE GRANDE ESCALA [HGE] website: https://www.ufrgs.br/hge/author/hge/
#UNIVERSITY - UNIVERSIDADE FEDERAL DO RIO GRANDE DO SUL - UFRGS
#RIO GRANDE DO SUL, BRAZIL

#DOI
#VERSION 0.1.1
#CONTACT US: leonardo.laipelt@ufrgs.br

#----------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------#

#PYTHON PACKAGES
#Call EE
# import ee
import datetime
from osgeo import gdal
import rasterio
import numpy as np
import os
import math
import pyproj
from rasterio.warp import reproject, Resampling

def save_data(image, output_dr, meta, array, name):
    with rasterio.open(os.path.join(output_dr, f'{name.lower()}.tif'), 'w', **meta) as dst:
        dst.write(array, 1)
    image[name.upper()] = os.path.join(output_dr, f'{name.lower()}.tif')

def get_meteorology(image, time_start, data_dr, cal_bands_dr, meta):
    col_meteorology = {}
    bands = ['ALFA']
    arrays = {}
    for band in bands:
        src = rasterio.open(image[band])
        array = src.read(1).astype(np.float32)
        array[array == src.nodata] = np.nan
        arrays[band] = array.copy()

    #LINEAR INTERPOLATION
    TIME_START_NUM = time_start
    PREVIOUS_TIME = time_start - datetime.timedelta(hours=1)
    NEXT_TIME = time_start + datetime.timedelta(hours=1)

    ds = gdal.Open(os.path.join(data_dr, "era5-hourly.grib"))
    grib_src = rasterio.open(os.path.join(data_dr, "era5-hourly.grib"))
    previous_bands = []
    next_bands = []
    for band in range(1, ds.RasterCount):
        seconds = ds.GetRasterBand(band).GetMetadata()['GRIB_VALID_TIME']
        time = datetime.datetime.fromtimestamp(int(seconds))
        if PREVIOUS_TIME <= time <= TIME_START_NUM:
            previous_bands.append(band)
        elif TIME_START_NUM <= time <= NEXT_TIME:
            next_bands.append(band)
    seconds = ds.GetRasterBand(previous_bands[0]).GetMetadata()['GRIB_VALID_TIME']
    IMAGE_PREVIOUS_TIME = datetime.datetime.fromtimestamp(int(seconds))
    seconds = ds.GetRasterBand(next_bands[0]).GetMetadata()['GRIB_VALID_TIME']
    IMAGE_NEXT_TIME = datetime.datetime.fromtimestamp(int(seconds))
    DELTA_TIME = (TIME_START_NUM - IMAGE_PREVIOUS_TIME) / (IMAGE_NEXT_TIME - IMAGE_PREVIOUS_TIME)

    #DAY OF THE YEAR
    doy = time_start.timetuple().tm_yday
    Pi= math.pi

    #INVERSE RELATIVE DISTANCE EARTH-SUN
    #ALLEN ET AL.(1998)
    d1 = 2 *Pi / 365
    d2 = d1 * doy
    d3 = np.cos(d2)
    dr = 1 + (0.033 * d3)

    #SOLAR DECLINATION [RADIANS]
    #ASCE REPORT (2005)
    e1 = 2 * Pi * doy
    e2 = e1 / 365
    e3 = e2 - 1.39
    e4 = np.sin(e3)
    solar_dec = 0.409 * e4
    
    #GET COORDINATES
    src = rasterio.open(image["ALFA"])
    height, width = src.height, src.width
    transform = src.transform
    cols, rows = np.meshgrid(np.arange(width), np.arange(height))
    xs, ys = transform * (cols, rows)
    proj = pyproj.Transformer.from_crs(src.crs, "EPSG:4326")  # Convert to WGS84
    lons, lats = proj.transform(xs, ys)
    lons[np.isnan(arrays["ALFA"])] = np.nan
    lats[np.isnan(arrays["ALFA"])] = np.nan
    save_data(image, cal_bands_dr, meta, lons, 'LONGITUDE')
    save_data(image, cal_bands_dr, meta, lats, 'LATITUDE')
    
    #SUNSET  HOUR ANGLE [RADIANS]
    #ASCE REPORT (2005)
    i_lat_rad = (lats * Pi) / 180
    i_sun_hour = np.arccos(- np.tan(i_lat_rad)* np.tan(solar_dec))    
    del(lons, lats)

    #SOLAR CONSTANT
    gsc = 4.92 #[MJ M-2 H-1]

    #EXTRATERRESTRIAL RADIATION 24H  [MJ M-2 D-1]
    #ASCE REPORT (2005)
    i_Ra_24h = (24/Pi)*gsc * dr * ( (i_sun_hour * np.sin(i_lat_rad)* np.sin(solar_dec)) +  (np.cos(i_lat_rad) * np.cos(solar_dec) * np.sin(i_sun_hour)))*11.5740
    del(i_sun_hour)

    i_Ra_24h= np.nanmean(i_Ra_24h)

    #INCOMING SHORT-WAVE RADIATION DAILY EAN [W M-2]
    sr_time = time_start - datetime.timedelta(hours=11)
    er_time = time_start + datetime.timedelta(hours=13)
    i_RS_sec = 0
    for band in range(1, ds.RasterCount):
        seconds = ds.GetRasterBand(band).GetMetadata()['GRIB_VALID_TIME']
        time = datetime.datetime.fromtimestamp(int(seconds))
        if sr_time <= time <= er_time:
            if ds.GetRasterBand(band).GetMetadata()['GRIB_COMMENT'] == 'Surface solar radiation downwards [W*s/m^2]':
                i_RS_sec = i_RS_sec + ds.GetRasterBand(band).ReadAsArray()
                # filtered_bands.append(band)
    i_Rs_24h = i_RS_sec / 86400
    i_Rs_24h, _ = reproject(i_Rs_24h, destination= np.zeros(src.shape), src_transform=grib_src.transform,
        src_crs=grib_src.crs, src_nodata=None, dst_transform=src.transform,
        dst_crs=src.crs, dst_nodata=None, dst_resolution=src.res, resampling=Resampling.bilinear)
    i_Rs_24h[np.isnan(arrays["ALFA"])] = np.nan
    save_data(col_meteorology, cal_bands_dr, meta, i_Rs_24h, 'SW_DOWN')

    # TASUMI
    i_albedo_ls = arrays['ALFA']

    #NET RADIATION 24H [W M-2]
    #BRUIN (1982)
    i_Rn_24h = ((1 - i_albedo_ls) * i_Rs_24h) - (110 * (i_Rs_24h / i_Ra_24h))
    save_data(col_meteorology, cal_bands_dr, meta, i_Rn_24h, 'RN24h_G')
    del(i_albedo_ls)
    
    next_image = {}
    for band in next_bands:
        comment = ds.GetRasterBand(band).GetMetadata()['GRIB_COMMENT']
        next_image[comment] = band

    previous_image = {}
    for band in previous_bands:
        comment = ds.GetRasterBand(band).GetMetadata()['GRIB_COMMENT']
        previous_image[comment] = band
    
    # AIR TEMPERATURE [K]
    tair_pre = ds.GetRasterBand(previous_image['2 metre temperature [C]']).ReadAsArray()
    tair_next = ds.GetRasterBand(next_image['2 metre temperature [C]']).ReadAsArray()
    tair_c = ((tair_next - tair_pre) * DELTA_TIME) + tair_pre
    
    tair_c, _ = reproject(tair_c, destination= np.zeros(src.shape), src_transform=grib_src.transform,
        src_crs=grib_src.crs, src_nodata=None, dst_transform=src.transform,
        dst_crs=src.crs, dst_nodata=None, dst_resolution=src.res, resampling=Resampling.bilinear)
    tair_c[np.isnan(arrays["ALFA"])] = np.nan
    del(tair_pre, tair_next)

    # WIND SPEED [M S-1]
    wind_u_pre = ds.GetRasterBand(previous_image['10 metre u wind component [m/s]']).ReadAsArray()
    wind_u_next = ds.GetRasterBand(next_image['10 metre u wind component [m/s]']).ReadAsArray()
    wind_u = ((wind_u_next - wind_u_pre) * DELTA_TIME) + wind_u_pre
    del(wind_u_pre, wind_u_next)

    wind_v_pre = ds.GetRasterBand(previous_image['10 metre v wind component [m/s]']).ReadAsArray()
    wind_v_next = ds.GetRasterBand(next_image['10 metre v wind component [m/s]']).ReadAsArray()
    wind_v = ((wind_v_next - wind_v_pre) * DELTA_TIME) + wind_v_pre
    del(wind_v_pre, wind_v_next)

    # TODO: CGM check if the select calls are needed
    wind_med = np.sqrt(wind_u ** 2 + wind_v ** 2)
    wind_med = wind_med * (4.87) / np.log(67.8 * 10 - 5.42)

    wind_med, _ = reproject(wind_med, destination= np.zeros(src.shape), src_transform=grib_src.transform,
        src_crs=grib_src.crs, src_nodata=None, dst_transform=src.transform,
        dst_crs=src.crs, dst_nodata=None, dst_resolution=src.res, resampling=Resampling.bilinear)
    wind_med[np.isnan(arrays["ALFA"])] = np.nan
    save_data(col_meteorology, cal_bands_dr, meta, wind_med, 'UX_G')
    del(wind_u, wind_v, wind_med)

    # PRESSURE [PA] CONVERTED TO KPA
    tdp_pre = ds.GetRasterBand(previous_image['2 metre dewpoint temperature [C]']).ReadAsArray()
    tdp_next = ds.GetRasterBand(next_image['2 metre dewpoint temperature [C]']).ReadAsArray()
    tdp = ((tdp_next - tdp_pre) * DELTA_TIME) + tdp_pre
    
    tdp, _ = reproject(tdp, destination= np.zeros(src.shape), src_transform=grib_src.transform,
        src_crs=grib_src.crs, src_nodata=None, dst_transform=src.transform,
        dst_crs=src.crs, dst_nodata=None, dst_resolution=src.res, resampling=Resampling.bilinear)
    tdp[np.isnan(arrays["ALFA"])] = np.nan
    del(tdp_next, tdp_pre)

    # ACTUAL VAPOR PRESSURE [KPA]
    ea = 0.6108 * (np.exp((17.27 * (tdp - 273.15)) / ((tdp - 273.15) + 237.3)))
    del(tdp)
    # SATURATED VAPOR PRESSURE [KPA]
    esat = 0.6108 * (np.exp((17.27 * (tair_c - 273.15)) / ((tair_c - 273.15) + 237.3)))

    # RELATIVE HUMIDITY (%)
    rh = ea / esat * 100
    save_data(col_meteorology, cal_bands_dr, meta, rh, 'RH_G')

    tair_c = tair_c - 273.15
    save_data(col_meteorology, cal_bands_dr, meta, tair_c, 'AIRT_G')
    del (arrays, tair_c, rh, esat, ea)

    return col_meteorology

if __name__ == "__main__":
    get_meteorology()
