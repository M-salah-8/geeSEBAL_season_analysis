import pandas as pd
import os
import glob
import numpy as np
import rasterio
import rasterio.mask
import datetime

def fill_missing(daily_df, data_name):
  daily_df.at[0,data_name] = daily_df.at[daily_df[data_name].first_valid_index(), data_name]
  # daily_df.at[0,data_name] = 0
  daily_df[data_name] = daily_df[data_name].interpolate('linear')
  # daily_df[data_name] = daily_df[data_name].interpolate(method='polynomial', order=3)
  # daily_df.fillna(method='ffill')

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

def write_dfs(season_dr, field_gdf):
  dekads = sorted(os.listdir(os.path.join(season_dr, "dekads")))
  os.makedirs(os.path.join(season_dr, 'sheets'), exist_ok= True)
  for index, row in field_gdf.iterrows():
    daily_df = pd.DataFrame({'date':[], 'ETp':[]})
    wapor_pet = sorted(glob.glob(os.path.join(season_dr, 'WaPOR_data', "L1-RET-E_resampled")+"/*.tif"))
    for pet in wapor_pet:
      date = os.path.basename(pet).split('.')[0].split('_')[-1]
      with rasterio.open(pet) as src:
        pet_array = rasterio.mask.mask(src, [row.geometry], crop = True, nodata= np.nan)[0][0]
      etp_row = {'date': date, 'ETp': np.nanmean(pet_array)}
      daily_df.loc[len(daily_df)] = etp_row
    # kc
    for dekad in dekads:
      tifs_dr = os.path.join(season_dr, "dekads", dekad, "kc", "tifs")
      if os.path.exists(tifs_dr):
        tifs = glob.glob(tifs_dr + "/*.tif")
        for kc in tifs:
          date = os.path.basename(kc).split(".")[0].split("_")[1]
          kc_src = rasterio.open(kc)
          kc_array = rasterio.mask.mask(kc_src, [row.geometry], crop = True, nodata= np.nan)[0][0]
          daily_df.loc[daily_df['date'] == date, 'kc'] = np.nanmean(kc_array)
      else:
        tifs_dr = os.path.join(season_dr, "dekads", dekad, "ETa")
        et_tif = glob.glob(tifs_dr + "/*.tif")[0]
        date = os.path.basename(et_tif).split(".")[0]
        et_src = rasterio.open(et_tif)
        et_array = rasterio.mask.mask(et_src, [row.geometry], crop = True, nodata= np.nan)[0][0]
        et_mean = np.nanmean(et_array)
        days = dekad_days(date)
        # et_day_mean = et_mean / len(days)
        pet_dekad = 0
        for day in days:
          pet = daily_df.loc[daily_df['date'] == day, 'ETp'].values[0]  # check
          pet_dekad = pet_dekad + pet
        kc_mean = et_mean / pet_dekad
        for day in days:
          daily_df.loc[daily_df['date'] == day, 'kc'] = kc_mean
      # biomass
      tifs_dr = os.path.join(season_dr, "dekads", dekad, "biomass")
      if os.path.exists(tifs_dr):
        biomass_tif = glob.glob(tifs_dr + "/*.tif")[0]
        date = os.path.basename(biomass_tif).split(".")[0]
        biomass_src = rasterio.open(biomass_tif)
        biomass_array = rasterio.mask.mask(biomass_src, [row.geometry], crop = True, nodata= np.nan)[0][0]
        biomass_mean = np.nanmean(biomass_array)
        days = dekad_days(date)
        biomass_day_mean = biomass_mean / len(days)
        for day in days:
          daily_df.loc[daily_df['date'] == day, 'biomass'] = biomass_day_mean
      else:
        continue
    # complite the daily dfs
    fill_missing(daily_df, 'kc')
    daily_df['ETa'] = daily_df['ETp'] * daily_df['kc']
    daily_df.to_csv(os.path.join(season_dr, 'sheets', f'daily_data_{str(row.id)}.csv'), index=False)