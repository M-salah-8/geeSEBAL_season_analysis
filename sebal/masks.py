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
import os
import rasterio
import numpy as np

#CLOUD REMOVAL

#FUNCTION TO MASK CLOUDS IN LANDSAT 5 AND 7 FOR SURFACE REFLECTANCE

def f_cloudMaskL457_SR(image):
    quality = image.select('pixel_qa')
    c01 = quality.eq(66)#CLEAR, LOW CONFIDENCE CLOUD
    c02 = quality.eq(68)#WATER, LOW CONFIDENCE CLOUD
    mask = c01.Or(c02)
    return image.updateMask(mask)

#FUNCTION FO MASK CLOUD IN LANDSAT 8 FOR SURFACE REFELCTANCE    ### reduce none values
def f_cloudMaskL8_SR(image):                                    ### check clear\water numbers
    quality = image.select('pixel_qa')
    c01 = quality.eq(322) #CLEAR, LOW CONFIDENCE CLOUD
    c02 = quality.eq(324) #WATER, LOW CONFIDENCE CLOUD
    c03 = quality.eq(1346)#CLEAR TERRAIN
    mask = c01.Or(c02).Or(c03)
    return image.updateMask(mask)

#ALBEDO
#USING TASUMI ET AL. (2008) METHOD FOR LANDSAT 8
#COEFFICIENTS FROM KE ET AL. (2016)
def f_albedoL_8_9(image, meta, cal_bands_dr):
  bands = ['UB', 'B', 'GR', 'R', 'NIR', 'SWIR_1', 'SWIR_2']
  arrays = {}
  for band in bands:
    src = rasterio.open(image[band])
    array = src.read(1).astype(np.float32)
    array[array == src.nodata] = np.nan
    arrays[band] = array.copy() * 0.0000275 - 0.2
  alfa = (0.130*arrays['UB']) + (0.115*arrays['B']) + (0.143*arrays['GR']) + (0.180*arrays['R']) + (0.281*arrays['NIR']) + (0.108*arrays['SWIR_1']) + (0.042*arrays['SWIR_2'])
  with rasterio.open(os.path.join(cal_bands_dr, 'alfa.tif'), 'w', **meta) as dst:
    dst.write(alfa, 1)
  image['ALFA'] = os.path.join(cal_bands_dr, 'alfa.tif')
  del(arrays, alfa)
  return image

if __name__ == "__main__":
    f_albedoL_8_9()
    f_cloudMaskL8_SR()
    f_cloudMaskL457_SR()
