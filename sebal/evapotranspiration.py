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
import rasterio
import numpy as np
import os

def save_data(image, output_dr, meta, array, name):
    with rasterio.open(os.path.join(output_dr, f'{name.lower()}.tif'), 'w', **meta) as dst:
        dst.write(array, 1)
    image[name.upper()] = os.path.join(output_dr, f'{name.lower()}.tif')

def fexp_et(image, Rn24hobs, cal_bands_dr, meta, results_dr, date_string):
    bands = ['RN', 'G', 'T_LST_DEM', 'H']
    arrays = {}
    for band in bands:
        src = rasterio.open(image[band])
        array = src.read(1).astype(np.float32)
        array[array == src.nodata] = np.nan
        arrays[band] = array.copy()
    src = rasterio.open(Rn24hobs)
    array = src.read(1).astype(np.float32)
    array[array == src.nodata] = np.nan
    arrays['RN24H_G'] = array.copy()

    #NET DAILY RADIATION (Rn24h) [W M-2]
    #BRUIN (1982)
    src = rasterio.open(Rn24hobs)
    array = src.read(1).astype(np.float32)
    array[array == src.nodata] = np.nan
    Rn24hobs = array.copy() * 1

    #GET ENERGY FLUXES VARIABLES AND LST
    i_Rn = arrays['RN']
    i_G = arrays['G']
    i_lst = arrays['T_LST_DEM']
    i_H_final = arrays['H']

    #FILTER VALUES
    i_H_final[i_H_final < 0]= 0

    # INSTANTANEOUS LATENT HEAT FLUX (LE) [W M-2]
    #BASTIAANSSEN ET AL. (1998)
    i_lambda_ET = i_Rn-i_G-i_H_final
    save_data(image, cal_bands_dr, meta, i_lambda_ET, 'LE')
    #FILTER
    i_lambda_E= np.where(i_lambda_ET < 0, 0, i_lambda_ET)
    del(i_lambda_E)

    #LATENT HEAT OF VAPORIZATION (LAMBDA) [J KG-1]
    #BISHT ET AL.(2005)
    #LAGOUARDE AND BURNET (1983)
    i_lambda = 2.501-0.002361*(i_lst-273.15)

    #INSTANTANEOUS ET (ET_inst) [MM H-1]
    i_ET_inst = 0.0036 * (i_lambda_ET/i_lambda)
    save_data(image, cal_bands_dr, meta, i_ET_inst, 'ET_INST')

    #EVAPORATIVE FRACTION (EF)
    #CRAGO (1996)
    i_EF = i_lambda_ET/(i_Rn-i_G)
    save_data(image, cal_bands_dr, meta, i_EF, 'EF')

    #DAILY EVAPOTRANSPIRATION (ET_24h) [MM DAY-1]
    i_ET24h_calc = (0.0864 *i_EF * Rn24hobs)/(i_lambda)
    save_data(image, cal_bands_dr, meta, i_ET24h_calc, 'ET_24h')
    output_dr = os.path.join(results_dr, 'eta')
    os.makedirs(output_dr, exist_ok=True)
    save_data(image, output_dr, meta, i_ET24h_calc, f'eta_{date_string}')

    #ADD BANDS
    del(arrays, i_ET24h_calc, i_EF, i_ET_inst, i_lambda, i_lambda_ET, i_H_final, i_lst, Rn24hobs, i_G, i_Rn)
    return image