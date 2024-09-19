# Preparation
import os, glob
import datetime
import geopandas as gpd
import pandas as pd
# read inputs from spreadsheet
working_dir = os.getcwd()
inputs_dr = os.path.join(working_dir, 'inputs.csv')
local_images_dr = os.path.join(working_dir, 'local_data')
inputs_df = pd.read_csv(inputs_dr)
# choose a season

# by index
index = 0
row = inputs_df.iloc[index]

# # by name
# name = '2024_winter'
# row = inputs_df.loc[inputs_df['season'] == name].iloc[0]

season = row['season']

# for local images
from products import Image_local

output_dr = os.path.join(working_dir, 'local_data', 'products_outputs', season, 'SEBAL')

local_ls_images = [os.path.split(i)[0] for i in glob.glob(os.path.join(local_images_dr, 'landsat', 'images', '*','*/'))]
for local_image in local_ls_images:
  landsat_img = Image_local(local_image, local_images_dr, output_dr,
                            NDVI_cold=1,
                            Ts_cold=1,
                            NDVI_hot=1,
                            Ts_hot=1)