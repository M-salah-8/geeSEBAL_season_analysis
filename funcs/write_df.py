import pandas as pd
import os
import glob
import numpy as np
import rasterio
import rasterio.mask
import scipy.stats
   
def fill_missing(df, daily_df, data_name):
  daily_df.at[0,data_name] = 0
  daily_df[data_name] = daily_df[data_name].interpolate(method='polynomial', order=3)
  date_2, date_1= [pd.to_datetime(i) for i in df['date'].tolist()[-2:]]
  last_dates_df = daily_df[(daily_df['date'] >= date_2) & (daily_df['date'] <= date_1)]
  dates_len = len(last_dates_df)
  x = range(1, dates_len+1)
  slope, intercept, _, _, _ = scipy.stats.linregress(x, last_dates_df[data_name])
  x = dates_len
  for index in daily_df[(daily_df['date'] > date_1)].index:
    x += 1
    daily_df.loc[index, data_name] = slope*x + intercept
    
def write_dfs(season_dr, daily_ET_df, field_shp):
    ETa_df = pd.DataFrame({'date':[], 'min':[], 'median':[], 'mean':[], 'max':[], 'stdDev':[]})
    kc_df = pd.DataFrame({'date':[], 'min':[], 'median':[], 'mean':[], 'max':[], 'stdDev':[]})
    ndvi_df = pd.DataFrame({'date':[], 'min':[], 'median':[], 'mean':[], 'max':[], 'stdDev':[]})
    daily_ndvi_df = daily_ET_df[['date']].copy()
    et_folder = os.path.join(season_dr, 'ET')
    ndvi_folder = os.path.join(season_dr, 'ndvi')
    ETa_tifs_paths = sorted(glob.glob(os.path.join(et_folder, 'tifs', 'ETa') + '/*.tif'))
    kc_tifs_paths = sorted(glob.glob(os.path.join(et_folder, 'tifs', 'kc') + '/*.tif'))
    ndvi_tifs_paths = sorted(glob.glob(os.path.join(ndvi_folder, 'tifs', 'ndvi') + '/*.tif'))
    
    for ETa_tif_path, kc_tif_path, ndvi_tif_path in zip(ETa_tifs_paths, kc_tifs_paths, ndvi_tifs_paths):
        date = os.path.basename(ETa_tif_path).split('.')[0].split('_')[1]
        with rasterio.open(ETa_tif_path) as src:
            ETa_array = rasterio.mask.mask(src, field_shp.to_crs(src.crs).geometry, crop = True, nodata= np.nan)[0][0]
        with rasterio.open(kc_tif_path) as src:
            kc_array = rasterio.mask.mask(src, field_shp.to_crs(src.crs).geometry, crop = True, nodata= np.nan)[0][0]
        with rasterio.open(ndvi_tif_path) as src:
            ndvi_array = rasterio.mask.mask(src, field_shp.to_crs(src.crs).geometry, crop = True, nodata= np.nan)[0][0]

        ETa_row = {'date': date, 'min': np.nanmin(ETa_array), 'median': np.nanmedian(ETa_array), 'mean': np.nanmean(ETa_array), 'max': np.nanmax(ETa_array), 'stdDev': np.nanstd(ETa_array)}
        kc_row = {'date': date, 'min': np.nanmin(kc_array), 'median': np.nanmedian(kc_array), 'mean': np.nanmean(kc_array), 'max': np.nanmax(kc_array), 'stdDev': np.nanstd(kc_array)}
        ndvi_row = {'date': date, 'min': np.nanmin(ndvi_array), 'median': np.nanmedian(ndvi_array), 'mean': np.nanmean(ndvi_array), 'max': np.nanmax(ndvi_array), 'stdDev': np.nanstd(ndvi_array)}
        ETa_df.loc[len(ETa_df)] = ETa_row
        kc_df.loc[len(kc_df)] = kc_row
        ndvi_df.loc[len(ndvi_df)] = ndvi_row
        # add to daily dfs
        daily_ET_df.loc[daily_ET_df['date'] == date, 'ETa'] = ETa_row['mean']
        daily_ET_df.loc[daily_ET_df['date'] == date, 'kc'] = kc_row['mean']
        daily_ndvi_df.loc[daily_ndvi_df['date'] == date, 'ndvi'] = ndvi_row['mean']

    # complite the daily dfs
    # daily_ET_df.at[0,'kc'] = daily_ET_df.at[daily_ET_df['kc'].first_valid_index(), 'kc']
    daily_ET_df.at[0,'kc'] = 0
    fill_missing(kc_df, daily_ET_df, 'kc')
    daily_ET_df['ETa'] = daily_ET_df['ETp'] * daily_ET_df['kc']
    fill_missing(ndvi_df, daily_ndvi_df, 'ndvi')
    # monthly
    monthly_ET_df = daily_ET_df.copy()
    monthly_ET_df['Month-Year'] = monthly_ET_df['date'].dt.to_period('M')
    monthly_ET_df = monthly_ET_df.groupby('Month-Year')['ETa'].agg(['sum', 'mean']).reset_index()
    monthly_ndvi_df = daily_ndvi_df.copy()
    monthly_ndvi_df['Month-Year'] = monthly_ndvi_df['date'].dt.to_period('M')
    monthly_ndvi_df = monthly_ndvi_df.groupby('Month-Year')['ndvi'].agg(['sum', 'mean']).reset_index()
    # export dataframes
    et_sheets_dr = os.path.join(et_folder,'datasheets')
    os.makedirs(et_sheets_dr, exist_ok= True)
    ETa_df.to_csv(os.path.join(et_sheets_dr,'ETa.csv'), index=False)
    kc_df.to_csv(os.path.join(et_sheets_dr,'kc.csv'), index=False)
    daily_ET_df.to_csv(os.path.join(et_sheets_dr,'daily_ET.csv'), index=False)
    monthly_ET_df.to_csv(os.path.join(et_sheets_dr,'monthly_ET.csv'), index=False)
    ndvi_sheet_dr = os.path.join(ndvi_folder, 'datasheets')
    os.makedirs(ndvi_sheet_dr, exist_ok= True)
    ndvi_df.to_csv(os.path.join(ndvi_sheet_dr, 'ndvi.csv'), index=False)
    daily_ndvi_df.to_csv(os.path.join(ndvi_sheet_dr, 'daily_ndvi.csv'), index=False)
    monthly_ndvi_df.to_csv(os.path.join(ndvi_sheet_dr,'monthly_ndvi.csv'), index=False)
    return daily_ET_df, monthly_ET_df, ETa_df, kc_df, daily_ndvi_df, monthly_ndvi_df, ndvi_df