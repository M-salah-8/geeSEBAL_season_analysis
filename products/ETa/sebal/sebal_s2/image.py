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
import xml.etree.ElementTree as ET
import zipfile
import pickle

#FOLDERS
from .masks import (f_cloudMask,f_albedo)
from .meteorology import get_meteorology
from .tools import (fexp_spec_ind, fexp_lst_export,fexp_radlong_up,
fexp_radshort_down, fexp_radlong_down, fexp_radbalance, fexp_soil_heat,fexp_sensible_heat_flux, co_dem)
from .endmembers import fexp_cold_pixel, fexp_hot_pixel
from .evapotranspiration import fexp_et
from .download_meteorology import download_era5_re_hourly
from .snap_process import run_sen_et

#IMAGE FUNCTION
class Image_s2():

    #ENDMEMBERS DEFAULT
    #ALLEN ET AL. (2013)
    def __init__(self,
                 s2_image_dir,
                 s3_image_dir,
                 gpt,
                 py39,
                 gdf,
                 local_data_dr,
                 results_dr,
                 NDVI_cold=5,
                 Ts_cold=20,
                 NDVI_hot=10,
                 Ts_hot=20):

        # metadata
        self.s2_image_dir = s2_image_dir
        self.s3_image_dir = s3_image_dir
        self.local_data_dr = local_data_dr
        self.image_dir = os.path.dirname(self.s2_image_dir)
        self.calculations_dr = os.path.join(self.image_dir, "calculated_bands")
        os.makedirs(self.calculations_dr, exist_ok= True)
        self.results_dr = results_dr        
        dt = os.path.basename(self.s3_image_dir).split("_")[-11]
        utc_timestamp = dt[:4] + "-" + dt[4:6] + "-" + dt[6:8] + " " + dt[9:11] + ":" + dt[11:13]
        self.date_string= dt[:4] + "-" + dt[4:6] + "-" + dt[6:8]
        self.time_start= datetime.datetime.strptime(utc_timestamp, "%Y-%m-%d %H:%M")
        self._year=self.time_start.year
        self._month=self.time_start.month
        self._day=self.time_start.day
        self._hour=self.time_start.hour
        self._minute = self.time_start.minute
        with zipfile.ZipFile(self.s2_image_dir, 'r') as zip:
            meta_file = [f for f in zip.namelist() if f.endswith('MTD_TL.xml')][0]
            with zip.open(meta_file, 'r') as f:
                tree = ET.parse(f)
                root = tree.getroot()
                zenith_angle = float(root.find(".//Mean_Sun_Angle/ZENITH_ANGLE").text)
                self.cloud_cover = float(root.find(".//CLOUDY_PIXEL_PERCENTAGE").text)
                self.sun_elevation = 90 - zenith_angle        

        #ENDMEMBERS
        self.p_top_NDVI=NDVI_cold
        self.p_coldest_Ts=Ts_cold
        self.p_lowest_NDVI=NDVI_hot
        self.p_hottest_Ts=Ts_hot

        with open(os.path.join('image.pkl'), 'rb') as f:
            self.image = pickle.load(f)
        with open(os.path.join('meteorology.pkl'), 'rb') as f:
            col_meteorology = pickle.load(f)

        strm_tif = os.path.join(self.local_data_dr, 'strm_30m.tif')
        reflectance, s2_mask, lai, lst, s3_mask, sharpened_LST, self.z_alt = run_sen_et(gpt, py39, self.s2_image_dir, self.s3_image_dir, self.calculations_dr, strm_tif, gdf)

        self.image = {}
        self.image['reflectance'] = reflectance
        self.image['LAI'] = lai
        self.image['LST'] = lst
        self.image['T_LST_DEM'] = sharpened_LST
        bands = {'B1': 'UB', 'B2': 'B', 'B3': 'GR', 'B4': 'R',
                'B5': "RG1", 'B6': "RG2", 'B7': "RG3", 'B8': 'NIR',
                'B11': 'SWIR_1', 'B12': 'SWIR_2'}

        for band in bands.keys():
            self.image[bands[band]] = os.path.join(reflectance, band + '.img')

        with rasterio.open(self.image['B']) as src:
            self.tif_meta = src.meta
            self.tif_meta.update(
                driver= "GTiff",
                dtype= np.float32,
                nodata= np.nan
            )
            self.res = src.res
            self.scale = src.scales[0]
            self.offset = src.offsets[0]
        
        # #CLOUD REMOVAL
        self.image=f_cloudMask(self.image, s2_mask, s3_mask, self.calculations_dr)

        # ALBEDO TASUMI ET AL. (2008) METHOD WITH KE ET AL. (2016) COEFFICIENTS
        self.image=f_albedo(self.image, self.tif_meta, self.calculations_dr, self.scale, self.offset)


        # METEOROLOGY PARAMETERS
        download_era5_re_hourly(self.image_dir, self.time_start, gdf)
        col_meteorology= get_meteorology(self.image, self.time_start, self.image_dir, self.calculations_dr, self.tif_meta)

        with open(os.path.join('meteorology.pkl'), 'wb') as f:
            pickle.dump(col_meteorology, f)

        #AIR TEMPERATURE [C]
        self.T_air = col_meteorology['AIRT_G']

        #WIND SPEED [M S-1]
        self.ux= col_meteorology['UX_G']

        #RELATIVE HUMIDITY [%]
        self.UR = col_meteorology['RH_G']

        #NET RADIATION 24H [W M-2]
        self.Rn24hobs = col_meteorology['RN24H_G']

        # SPECTRAL IMAGES (NDVI, EVI, SAVI, LAI, T_LST, e_0, e_NB, long, lat)
        self.image=fexp_spec_ind(self.image, self.tif_meta, self.calculations_dr, self.results_dr, self.date_string, self.scale, self.offset)

        with open(os.path.join('image.pkl'), 'wb') as f:
            pickle.dump(self.image, f)

        #COLD PIXEL
        self.d_cold_pixel=fexp_cold_pixel(self.image, self.p_top_NDVI, self.p_coldest_Ts)

        #COLD PIXEL NUMBER
        self.n_Ts_cold = self.d_cold_pixel['temp']

        # INSTANTANEOUS OUTGOING LONG-WAVE RADIATION [W M-2]
        self.image=fexp_radlong_up(self.image, self.calculations_dr, self.tif_meta)

        #INSTANTANEOUS INCOMING SHORT-WAVE RADIATION [W M-2]
        self.image=fexp_radshort_down(self.image,self.z_alt,self.T_air,self.UR, self.sun_elevation, self.time_start,self.tif_meta, self.calculations_dr)

        #INSTANTANEOUS INCOMING LONGWAVE RADIATION [W M-2]
        self.image=fexp_radlong_down(self.image,  self.n_Ts_cold, self.calculations_dr, self.tif_meta)

        #INSTANTANEOUS NET RADIATON BALANCE [W M-2]
        self.image=fexp_radbalance(self.image, self.calculations_dr, self.tif_meta)

        #SOIL HEAT FLUX (G) [W M-2]
        self.image=fexp_soil_heat(self.image, self.calculations_dr, self.tif_meta)

        #HOT PIXEL
        self.d_hot_pixel=fexp_hot_pixel(self.image, self.p_lowest_NDVI, self.p_hottest_Ts)

        #SENSIBLE HEAT FLUX (H) [W M-2]
        self.image=fexp_sensible_heat_flux(self.image, self.ux, self.UR,self.Rn24hobs,self.n_Ts_cold, self.d_hot_pixel, self.calculations_dr, self.tif_meta)

        #DAILY EVAPOTRANSPIRATION (ET_24H) [MM DAY-1]
        self.image=fexp_et(self.image, self.Rn24hobs, self.calculations_dr, self.tif_meta, self.results_dr, self.date_string)

        shutil.rmtree(self.calculations_dr)