import ee
import requests
import os
import glob
import rasterio
import pandas as pd
import datetime
import rasterio.mask
import geopandas as gpd
from wapordl import wapor_map

def date_dekad(date):
  day = int(date.split("-")[-1])
  if day < 11:
    return date[:-2]+"01"
  elif day < 21:
    return date[:-2]+"11"
  elif day < 32:
    return date[:-2]+"21"
  
def dekad_days(date):
  day = int(date.split("-")[-1])
  if day == 1 or day == 11:
    s_d = datetime.datetime.strptime(date, '%Y-%m-%d')
    e_d = s_d + datetime.timedelta(days=10)
    days = pd.date_range(s_d,e_d - datetime.timedelta(days=1),freq='d').astype(str)
    return days
  elif day == 21:
    s_d = datetime.datetime.strptime(date, '%Y-%m-%d')
    e_d = (s_d.replace(day=1) + datetime.timedelta(days=32)).replace(day=1)
    days = pd.date_range(s_d,e_d - datetime.timedelta(days=1),freq='d').astype(str)
    return days

def export_gee_tifs_localy(image, data, data_name, season_dr, date, aoi):
    date_d = date_dekad(date)
    data_url = {
        'NDVI': image.select('NDVI').getDownloadUrl({
            'bands': 'NDVI',
            'scale': 30,
            'region': aoi,
            'format': 'GEO_TIFF'
            }),
        'ET_24h': image.select('ET_24h').getDownloadUrl({
            'bands': 'ET_24h',
            'scale': 30,
            'region': aoi,
            'format': 'GEO_TIFF'
            }),
        'RGB': image.select(['R','GR','B']).multiply(0.0000275).add(-0.2).getDownloadUrl({
            'bands': ['R','GR','B'],
            'scale': 30,
            'region': aoi,
            'format': 'GEO_TIFF'
            })
    }
    download_dr = os.path.join(season_dr, "dekads", date_d, data_name, 'tifs')
    os.makedirs(download_dr, exist_ok=True)
    response = requests.get(data_url[data])
    with open(os.path.join(download_dr, f'{data_name}_{date}.tif'), 'wb') as fd:
        fd.write(response.content)

def export_image_to_drive(img, description, folder, aoi):
    # Export cloud-optimized GeoTIFF images
    ee.batch.Export.image.toDrive(**{
        'image': img,
        'description': description,
        'scale': 30,
        "folder": folder,
        'region': aoi,
        'fileFormat': 'GeoTIFF',
        'maxPixels': 3784216672400,
        'formatOptions': {
            'cloudOptimized': True
        }
    }).start()

def export_to_drive(result_img, date, folder, aoi):  ### fix date (but outside fun)
    export_image_to_drive(
        result_img.image.select(['R','GR','B']).multiply(0.0000275).add(-0.2),  ### check
        "rgb_" + date,
        folder,
        aoi
    )
    export_image_to_drive(
        result_img.image.select('NDVI'),
        "ndvi_" + date,
        folder,
        aoi
    )
    export_image_to_drive(
        result_img.image.select('ET_24h'),
        "eta_" + date,
        folder,
        aoi
    )

def download_WaPOR(region, variables, period, season_dr, overview = "NONE"):
  for var in variables:
    download_dr = os.path.join(season_dr, 'WaPOR_data', var)
    os.makedirs(download_dr, exist_ok=True)

    if('-E' in var):
      unit = "day"
    elif('-D' in var):
      unit = "dekad"
    elif('-M' in var):
      unit = "month"
    elif ('-A' in var):
      unit = "year"
    else:
      unit = "none"

    wapor_map(region, var, period, download_dr, seperate_unscale = True, unit_conversion = unit)

def kc_tifs(season_dr):
  wapor_pet = sorted(os.listdir(os.path.join(season_dr, "WaPOR_data", "L1-RET-E_resampled")))
  dekads = sorted(os.listdir(season_dr, "dekads"))
  for dekad in dekads:
    tifs_dr = os.path.join(season_dr, "dekads", dekad, "ETa", "tifs")
    if os.path.exists(tifs_dr):
      tifs = glob.glob(tifs_dr + "/*.tif")
      for et_tif in tifs:
        date = os.path.basename(et_tif).split(".")[0].split("_")[1]
        for pet in wapor_pet:
          if date in pet:
            print(pet)
            pet_tif = os.path.join(season_dr, "WaPOR_data", "L1-RET-E_resampled", pet)
            et_src = rasterio.open(et_tif)
            et_arr = et_src.read(1)
            pet_src = rasterio.open(pet_tif)
            pet_arr = pet_src.read(1)
            kc_arr = et_arr / pet_arr
            os.makedirs(os.path.join(season_dr, "dekads", dekad, "kc", 'tifs'), exist_ok=True)
            with rasterio.open(os.path.join(season_dr, "dekads", dekad, "kc", 'tifs', f'kc_{date}.tif'), 'w', **et_src.meta) as dst:
                dst.write(kc_arr, 1)
            break
    else:
      continue

def dekads_tifs(season_dr, project_gdf):
    dekads = sorted(os.listdir(os.path.join(season_dr, "dekads")))
    project_df = pd.read_csv(os.path.join(season_dr, 'sheets', 'daily_data_0.csv'))
    for dekad in dekads:
        tifs_dr = os.path.join(season_dr, "dekads", dekad, "ETa", "tifs")
        if os.path.exists(tifs_dr):
            month_arrays_mean_sum = 0
            tifs = glob.glob(tifs_dr + "/*.tif")
            for eta_tif in tifs:
                date = os.path.basename(eta_tif).split(".")[0].split("_")[1]
                eta_src = rasterio.open(eta_tif)
                eta_array = rasterio.mask.mask(eta_src, project_gdf.geometry, crop = True, nodata= np.nan)[0][0]
                day_mean = project_df.loc[project_df['date'] == date, "ETa"].values[0]
                array = eta_array / day_mean
                month_arrays_mean_sum = month_arrays_mean_sum + array
            month_arrays_mean = month_arrays_mean_sum / len(tifs)
            days = dekad_days(dekad)
            dekad_mean = 0
            for day in days:
                dekad_mean = dekad_mean + project_df.loc[project_df['date'] == day, "ETa"].values[0]
            with rasterio.open(os.path.join(season_dr, "dekads", dekad, "ETa", f'{dekad}.tif'), 'w', **eta_src.meta) as dst:
                dst.write(month_arrays_mean * dekad_mean, 1)
        else:
            continue
        
def export_local_tifs(data_drs, season_dr, project_gdf):
  for data_dr in data_drs:
    tifs_drs = glob.glob(data_dr+'/*.tif')
    for tif in tifs_drs:
      date = os.path.basename(tif).split(".")[0].split("_")[1]
      data_name = os.path.basename(tif).split(".")[0].split("_")[0]
      date_d = date_dekad(date)
      src = rasterio.open(tif)
      pet_array, transform = rasterio.mask.mask(src, project_gdf.to_crs(src.crs).geometry, crop = True)
      meta = src.meta
      meta.update(width= pet_array.shape[-1], height= pet_array.shape[-2], transform= transform)
      download_dr = os.path.join(season_dr, "dekads", date_d, data_name, 'tifs')
      os.makedirs(download_dr, exist_ok=True)
      with rasterio.open(os.path.join(download_dr, f'{data_name}_{date}.tif'), 'w', **meta) as dst:
        dst.write(pet_array)