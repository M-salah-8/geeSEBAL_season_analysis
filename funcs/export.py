import ee
import requests
import os
import glob
import rasterio

def export_tif_localy(image, data, data_name, data_folder, date, aoi):
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
    download_dr = os.path.join(data_folder, 'tifs', data_name)
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
        "RBG_" + date,
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
        "ET_" + date,
        folder,
        aoi
    )

def kc_tifs(et_folder, daily_ET_df):
    ETa_tifs_paths = sorted(glob.glob(os.path.join(et_folder, 'tifs', 'ETa') + '/*.tif'))
    for tif_path in ETa_tifs_paths:
        date = os.path.basename(tif_path).split('.')[0].split('_')[1]
        date_ETp = daily_ET_df.loc[daily_ET_df['date'] == date, 'ETp'].values[0]
        with rasterio.open(tif_path) as src:
            tif_array = src.read(1)
        kc_array = tif_array / date_ETp
        os.makedirs(os.path.join(et_folder, 'tifs', 'kc'), exist_ok= True)
        with rasterio.open(os.path.join(os.path.join(et_folder, 'tifs', 'kc'), f'kc_{date}.tif'), 'w', **src.meta) as dst:
            dst.write(kc_array, 1)

def monthly_seasonal_tifs(data_folder, data_name, season, daily_df, monthly_df):
    # monthly_tifs
    os.makedirs(os.path.join(data_folder, 'tifs', f'{data_name}_month'), exist_ok= True)
    data_tifs_paths = sorted(glob.glob(os.path.join(data_folder, 'tifs', data_name)+'/*.tif'))
    dates = [os.path.basename(i).split('.')[0].split('_')[1] for i in data_tifs_paths]
    months = ['-'.join(i.split('-')[:2]) for i in dates]
    months_data = {}
    for filename, month in zip(data_tifs_paths, months):
        months_data.setdefault(month, []).append(filename)
    for month in months_data:
        month_arrays_sum = 0
        for tif_path in months_data[month]:
            # get image date
            date = os.path.basename(tif_path).split('.')[0].split('_')[1]
            month = '-'.join(date.split("-")[:2])
            with rasterio.open(tif_path) as src:
                array = src.read(1)
            day_mean = daily_df.loc[daily_df['date'] == date, data_name].values[0]
            array = array / day_mean
            month_arrays_sum = month_arrays_sum + array
        month_mean = monthly_df.loc[monthly_df['Month-Year'] == month, 'sum'].values[0]
        month_array = (month_arrays_sum / len(months_data[month])) * month_mean
        with rasterio.open(os.path.join(data_folder, 'tifs', f'{data_name}_month', f'{data_name}_month_{month}.tif'), 'w', **src.meta) as dst:
            dst.write(month_array, 1)
    # seasonal tif
    monthly_tifs = sorted(glob.glob(os.path.join(data_folder, 'tifs', f'{data_name}_month')+'/*.tif'))
    season_array = 0
    for monthly_tif in monthly_tifs:
        with rasterio.open(monthly_tif) as src:
            array = src.read(1)
        season_array += array
        os.makedirs(os.path.join(data_folder, 'tifs', f'{data_name}_season'),exist_ok= True)
        with rasterio.open(os.path.join(data_folder, 'tifs', f'{data_name}_season', f'{data_name}_season_{season}.tif'), 'w', **src.meta) as dst:
            dst.write(season_array, 1)