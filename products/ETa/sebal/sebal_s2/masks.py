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
import os, sys
sys.path.append(os.getcwd())
import rasterio
import numpy as np
from funcs import geo


#CLOUD REMOVAL
def f_cloudMask(image, s2_mask, s3_mask, cal_bands_dr):
  geo_in = geo()
  with rasterio.open(s2_mask) as src_s2:
      mask_s2 = src_s2.read(1).astype(bool)
  mask_s3,_ = geo_in.coregister(s2_mask, s3_mask, resampling= 'nearest', dtype=np.int16)
  mask_s3 = mask_s3.astype(bool)[0]
  # make mask_2 and mask_3
  mask = np.logical_and(mask_s2, mask_s3)
  np.save(os.path.join(cal_bands_dr, 'mask.npy'), mask)
  image['MASK'] = os.path.join(cal_bands_dr, 'mask.npy')
  return image

#ALBEDO
#USING TASUMI ET AL. (2008) METHOD FOR LANDSAT 8
#COEFFICIENTS FROM KE ET AL. (2016)
def f_albedo(image, meta, cal_bands_dr, scale, offset):
  image_mask = np.load(image['MASK'])
  bands = ["B", "GR", "R", "RG1", "RG2", "RG3", "NIR", "SWIR_1", "SWIR_2"]
  arrays = {}
  for band in bands:
    with rasterio.open(image[band]) as src:
      array = src.read(1).astype(np.float32)
      array[~image_mask] = np.nan
      arrays[band] = array * scale + offset
  alfa = (0.1324*arrays["B"]) + (0.1269*arrays["GR"]) + (0.1051*arrays["R"]) + (0.0971*arrays["RG1"]) + (0.0890*arrays["RG2"]) + (0.0818*arrays["RG3"]) + (0.0722*arrays["NIR"]) + ( 0.0167*arrays["SWIR_1"]) + ( 0.0002*arrays["SWIR_2"])
  with rasterio.open(os.path.join(cal_bands_dr, 'alfa.tif'), 'w', **meta) as dst:
    dst.write(alfa, 1)
  image['ALFA'] = os.path.join(cal_bands_dr, 'alfa.tif')
  del(arrays, alfa)
  return image

if __name__ == "__main__":
    f_albedo()
    f_cloudMask()