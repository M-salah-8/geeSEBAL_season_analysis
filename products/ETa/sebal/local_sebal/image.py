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
import os
import datetime
import rasterio
import numpy as np
import shutil
import pickle

#FOLDERS
from .masks import (f_cloudMaskL8_SR,f_albedoL_8_9)
from .meteorology import get_meteorology
from .tools import (fexp_spec_ind, fexp_lst_export,fexp_radlong_up, LST_DEM_correction,
fexp_radshort_down, fexp_radlong_down, fexp_radbalance, fexp_soil_heat,fexp_sensible_heat_flux, co_dem)
from .endmembers import fexp_cold_pixel, fexp_hot_pixel
from .evapotranspiration import fexp_et
from .download_meteorology import download_era5_re_hourly

#IMAGE FUNCTION
class Image_local():

    #ENDMEMBERS DEFAULT
    #ALLEN ET AL. (2013)
    def __init__(self,
                 image_dr,
                 local_data_dr,
                 results_dr,
                 NDVI_cold=5,
                 Ts_cold=20,
                 NDVI_hot=10,
                 Ts_hot=20):
                #GET INFORMATIONS FROM IMAGE
        self.image_dr = image_dr
        self.local_data_dr = local_data_dr
        self.ls_data_dr = os.path.split(image_dr)[0]
        self.cal_bands_dr = os.path.join(self.image_dr, "calculated_bands")
        os.makedirs(self.cal_bands_dr, exist_ok= True)
        self.results_dr = results_dr
        meta_names = ['LANDSAT_PRODUCT_ID', 'SPACECRAFT_ID', 'SUN_ELEVATION', 'CLOUD_COVER', 'SCENE_CENTER_TIME', 'DATE_ACQUIRED']
        self.meta = {}
        for file in os.listdir(self.image_dr):
            if file.endswith('MTL.txt'):
                    with open(os.path.join(self.image_dr, file), 'r') as file:
                        for line in file:
                            if any(meta_name in line for meta_name in meta_names):
                                meta_name = line.split("=")[0].strip()
                                self.meta[meta_name] = line.split("=")[1].strip().replace('"', '')
                            if "LEVEL2_PROCESSING_RECORD" in line:
                                break
        #GET INFORMATIONS FROM IMAGE
        # self.image = ee.Image(image)
        self._index=self.meta['LANDSAT_PRODUCT_ID']
        self.cloud_cover=float(self.meta['CLOUD_COVER'])
        self.landsat_version=self.meta['SPACECRAFT_ID']
        self.sun_elevation=float(self.meta['SUN_ELEVATION'])
        utc_timestamp =f"{self.meta['DATE_ACQUIRED']} {self.meta['SCENE_CENTER_TIME'][:-2]}"
        self.time_start=datetime.datetime.strptime(utc_timestamp, "%Y-%m-%d %H:%M:%S.%f")
        self._year=self.time_start.year
        self._month=self.time_start.month
        self._day=self.time_start.day
        self._hour=self.time_start.hour
        self._minute = self.time_start.minute
        self.date_string=self.meta['DATE_ACQUIRED']

        #ENDMEMBERS
        self.p_top_NDVI=NDVI_cold
        self.p_coldest_Ts=Ts_cold
        self.p_lowest_NDVI=NDVI_hot
        self.p_hottest_Ts=Ts_hot

        #LANDSAT IMAGE
        if self.landsat_version == 'LANDSAT_8':
            bands = {'SR_B1': 'UB', 'SR_B2': 'B', 'SR_B3': 'GR', 'SR_B4': 'R',
                     'SR_B5': 'NIR', 'SR_B6': 'SWIR_1', 'SR_B7': 'SWIR_2',
                     'ST_B10': 'BRT', 'QA_PIXEL': 'pixel_qa'}
            self.image = {}
            for file in os.listdir(self.image_dr):
                if any(band in file for band in bands):
                    band_name = "_".join(file.split("_")[-2:]).split(".")[0]
                    self.image[bands[band_name]] = os.path.join(self.image_dr, file)
            self.ls_meta = rasterio.open(self.image['UB']).meta
            self.ls_meta.update(dtype= np.float32, nodata= np.nan)
            self.res = rasterio.open(self.image['UB']).res
            #CLOUD REMOVAL
            self.image=f_cloudMaskL8_SR(self.image, self.cal_bands_dr)

            # ALBEDO TASUMI ET AL. (2008) METHOD WITH KE ET AL. (2016) COEFFICIENTS
            self.image=f_albedoL_8_9(self.image, self.ls_meta, self.cal_bands_dr)

        elif self.landsat_version == 'LANDSAT_9':
            bands = {'SR_B1': 'UB', 'SR_B2': 'B', 'SR_B3': 'GR', 'SR_B4': 'R',
                     'SR_B5': 'NIR', 'SR_B6': 'SWIR_1', 'SR_B7': 'SWIR_2',
                     'ST_B10': 'BRT', 'QA_PIXEL': 'pixel_qa'}
            self.image = {}
            for file in os.listdir(self.image_dr):
                if any(band in file for band in bands):
                    band_name = "_".join(file.split("_")[-2:]).split(".")[0]
                    self.image[bands[band_name]] = os.path.join(self.image_dr, file)
            self.ls_meta = rasterio.open(self.image['UB']).meta
            self.ls_meta.update(dtype= np.float32, nodata= np.nan)
            self.res = rasterio.open(self.image['UB']).res
            #CLOUD REMOVAL
            self.image=f_cloudMaskL8_SR(self.image, self.cal_bands_dr)

            # ALBEDO TASUMI ET AL. (2008) METHOD WITH KE ET AL. (2016) COEFFICIENTS
            self.image=f_albedoL_8_9(self.image, self.ls_meta, self.cal_bands_dr)

        else:
            print('version error')

        # METEOROLOGY PARAMETERS
        download_era5_re_hourly(self.image['UB'], self.time_start)
        col_meteorology= get_meteorology(self.image, self.time_start, self.ls_data_dr, self.cal_bands_dr, self.ls_meta)

        #AIR TEMPERATURE [C]
        self.T_air = col_meteorology['AIRT_G']

        #WIND SPEED [M S-1]
        self.ux= col_meteorology['UX_G']

        #RELATIVE HUMIDITY [%]
        self.UR = col_meteorology['RH_G']

        #NET RADIATION 24H [W M-2]
        self.Rn24hobs = col_meteorology['RN24H_G']

        #SRTM DATA ELEVATION
        SRTM_ELEVATION = os.path.join(self.local_data_dr, 'strm_30m.tif')
        self.z_alt = co_dem(self.ls_meta.copy(), self.res, SRTM_ELEVATION, self.ls_data_dr)

        # SPECTRAL IMAGES (NDVI, EVI, SAVI, LAI, T_LST, e_0, e_NB, long, lat)
        self.image=fexp_spec_ind(self.image, self.ls_meta, self.cal_bands_dr, self.results_dr, self.date_string)

        #LAND SURFACE TEMPERATURE
        self.image=LST_DEM_correction(self.image, self.z_alt, self.T_air, self.UR,self.sun_elevation,self.time_start, self._hour,self._minute, self.ls_meta, self.cal_bands_dr)

        #COLD PIXEL
        self.d_cold_pixel=fexp_cold_pixel(self.image, self.p_top_NDVI, self.p_coldest_Ts)

        #COLD PIXEL NUMBER
        self.n_Ts_cold = self.d_cold_pixel['temp']

        # INSTANTANEOUS OUTGOING LONG-WAVE RADIATION [W M-2]
        self.image=fexp_radlong_up(self.image, self.cal_bands_dr, self.ls_meta)

        #INSTANTANEOUS INCOMING SHORT-WAVE RADIATION [W M-2]
        self.image=fexp_radshort_down(self.image,self.z_alt,self.T_air,self.UR, self.sun_elevation, self.time_start,self.ls_meta, self.cal_bands_dr)

        #INSTANTANEOUS INCOMING LONGWAVE RADIATION [W M-2]
        self.image=fexp_radlong_down(self.image,  self.n_Ts_cold, self.cal_bands_dr, self.ls_meta)

        #INSTANTANEOUS NET RADIATON BALANCE [W M-2]
        self.image=fexp_radbalance(self.image, self.cal_bands_dr, self.ls_meta)

        #SOIL HEAT FLUX (G) [W M-2]
        self.image=fexp_soil_heat(self.image, self.cal_bands_dr, self.ls_meta)

        #HOT PIXEL
        self.d_hot_pixel=fexp_hot_pixel(self.image, self.p_lowest_NDVI, self.p_hottest_Ts)

        #SENSIBLE HEAT FLUX (H) [W M-2]
        self.image=fexp_sensible_heat_flux(self.image, self.ux, self.UR,self.Rn24hobs,self.n_Ts_cold, self.d_hot_pixel, self.cal_bands_dr, self.ls_meta)

        #DAILY EVAPOTRANSPIRATION (ET_24H) [MM DAY-1]
        self.image=fexp_et(self.image, self.Rn24hobs, self.cal_bands_dr, self.ls_meta, self.results_dr, self.date_string)

        shutil.rmtree(self.cal_bands_dr)