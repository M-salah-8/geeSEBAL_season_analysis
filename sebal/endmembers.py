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
import ee
import rasterio
import numpy as np

#A SIMPLIFIED VERSION OF
#CALIBRATION USING INVERSE MODELING AT EXTREME CONDITIONS (CIMEC)
#FROM ALLEN ET AL. (2013) FOR METRIC
#SEE MORE: LAIPELT ET AL. (2020)

#DEFAULT PARAMETERS
#NDVI COLD = 5%
#TS COLD = 20%
#NDVI HOT = 10%
#TS HOT = 20%

#SELECT COLD PIXEL
def fexp_cold_pixel(image, p_top_NDVI, p_coldest_Ts):
  bands = ['NDVI_NEG', 'LST_NW', 'NDVI']
  arrays = {}
  for band in bands:
    src = rasterio.open(image[band])
    array = src.read(1).astype(np.float32)
    array[array == src.nodata] = np.nan
    arrays[band] = array.copy()

  #IDENTIFY THE TOP % NDVI PIXELS
  n_perc_top_NDVI= np.nanpercentile(arrays['NDVI_NEG'], p_top_NDVI)

  #UPDATE MASK WITH NDVI VALUES
  i_top_NDVI= arrays['NDVI_NEG'].copy()
  i_top_NDVI[i_top_NDVI > n_perc_top_NDVI] = np.nan
  #SELECT THE COLDEST TS FROM PREVIOUS NDVI GROUP
  lst_top_NDVI= np.where(np.isnan(i_top_NDVI), np.nan, arrays['LST_NW'])
  n_perc_low_LST= np.nanpercentile(lst_top_NDVI, p_coldest_Ts)

  #UPDATE MASK WITH LST VALUES
  i_cold_lst= lst_top_NDVI.copy()
  i_cold_lst[i_cold_lst > n_perc_low_LST] = np.nan

  #FILTERS    ### ??
  c_lst_cold20 = i_cold_lst.copy()
  c_lst_cold20[c_lst_cold20 < 200] = np.nan
  c_lst_cold20_int=np.round(c_lst_cold20)

  #COUNT NUNMBER OF PIXELS
  n_count_final_cold_pix = np.count_nonzero(~np.isnan(c_lst_cold20_int))

  #SELECT COLD PIXEL RANDOMLY (FROM PREVIOUS SELECTION)
  non_nan_indices = np.where(~np.isnan(c_lst_cold20_int))
  index = np.random.choice(len(non_nan_indices[0]))
  i_0, i_1 = non_nan_indices[0][index], non_nan_indices[1][index]

  n_Ts_cold = arrays['LST_NW'][i_0, i_1]
  # n_long_cold = ee.Number(fc_cold_pix.aggregate_first('longitude'))
  # n_lat_cold = ee.Number(fc_cold_pix.aggregate_first('latitude'))
  n_ndvi_cold = arrays['NDVI'][i_0, i_1]

  #CREATE A DICTIONARY WITH THOSE RESULTS
  d_cold_pixel = {
          'temp': n_Ts_cold,
          'ndvi': n_ndvi_cold,
          'index': [i_0, i_1],
          'sum': n_count_final_cold_pix}

  del(arrays, n_perc_top_NDVI, i_top_NDVI, lst_top_NDVI, n_perc_low_LST, i_cold_lst, c_lst_cold20, n_count_final_cold_pix, non_nan_indices, index, n_ndvi_cold, n_Ts_cold)
  #RETURN DICTIONARY
  return d_cold_pixel

#SELECT HOT PIXEL
def fexp_hot_pixel(image, p_lowest_NDVI, p_hottest_Ts):
  bands = ['POS_NDVI', 'LST_NEG', 'INT', 'G', 'RN', 'NDVI', 'NDVI_NEG', 'LST_NW']
  arrays = {}
  for band in bands:
    src = rasterio.open(image[band])
    array = src.read(1).astype(np.float32)
    array[array == src.nodata] = np.nan
    arrays[band] = array.copy()

  #IDENTIFY THE DOWN % NDVI PIXELS
  n_perc_low_NDVI= np.nanpercentile(arrays['POS_NDVI'], p_lowest_NDVI)

  #UPDATE MASK WITH NDVI VALUES
  i_low_NDVI= arrays['NDVI_NEG'].copy()
  i_low_NDVI[i_low_NDVI > n_perc_low_NDVI] = np.nan

  #SELECT THE HOTTEST TS FROM PREVIOUS NDVI GROUP
  lst_low_NDVI= np.where(np.isnan(i_low_NDVI), np.nan, arrays['LST_NEG'])
  n_perc_top_lst= np.nanpercentile(lst_low_NDVI, p_hottest_Ts)
  
  c_lst_hotpix= lst_low_NDVI.copy()
  c_lst_hotpix[c_lst_hotpix > n_perc_top_lst] = np.nan
  c_lst_hotpix_int=np.round(c_lst_hotpix)
  
  #COUNT NUNMBER OF PIXELS
  n_count_final_hot_pix = np.count_nonzero(~np.isnan(c_lst_hotpix_int))

  #SELECT HOT PIXEL RANDOMLY (FROM PREVIOUS SELECTION)
  non_nan_indices = np.where(~np.isnan(c_lst_hotpix_int))
  index = np.random.choice(len(non_nan_indices[0]))
  i_0, i_1 = non_nan_indices[0][index], non_nan_indices[1][index]

  n_Ts_hot = arrays['LST_NW'][i_0, i_1]
  # n_long_hot = ee.Number(fc_hot_pix.aggregate_first('longitude'))
  # n_lat_hot = ee.Number(fc_hot_pix.aggregate_first('latitude'))  n_ndvi_hot = arrays['NDVI'][i_0, i_1]
  n_ndvi_hot = arrays['NDVI'][i_0, i_1]
  n_Rn_hot = arrays['RN'][i_0, i_1]
  n_G_hot = arrays['G'][i_0, i_1]
  #CREATE A DICTIONARY WITH THOSE RESULTS
  d_hot_pixel = {
        'temp': n_Ts_hot,
        'index': [i_0, i_1],
        'Rn': n_Rn_hot,
        'G': n_G_hot,
        'ndvi': n_ndvi_hot,
        'sum': n_count_final_hot_pix}

  del(arrays, i_low_NDVI, lst_low_NDVI, c_lst_hotpix, c_lst_hotpix_int, non_nan_indices)
  #RETURN DICTIONARY
  return d_hot_pixel