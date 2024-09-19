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
project_area_dr = row['project_area']
field_area_dr = row['field_area']
PET_area_dr = row['PET_area']
start_date = row['start_date']
end_date = row['end_date']

project_gdf = gpd.read_file(project_area_dr)
field_gdf = gpd.read_file(field_area_dr)
PET_gdf = gpd.read_file(PET_area_dr)

from products import sen_et
# inputs
gpt = "/home/msalah/esa-snap/bin/gpt"		# path to gpt executable
# path to python environment configured to be used with sen-et
py39 = "/home/msalah/.pyenv/versions/3.9.19/bin/python3"
s2_dir = os.path.join(working_dir, 'local_data', 'sentinel_2')
s3_dir = os.path.join(working_dir, 'local_data', 'sentinel_3')
lcc_tif = os.path.join(working_dir, 'local_data', "LCC.tif") # land cover tif
c_l = 42 # crop land in lcc
outputs_dir = os.path.join(working_dir, 'local_data', 'products_outputs', season, 'sen-et')

s2_images = [ 
 '/home/msalah/Desktop/m salah/crop monitoring/crop_monitoring_code/local_data/sentinel_2/images/2024_01_23/S2B_MSIL2A_20240123T083149_N0510_R021_T35QRA_20240123T100647.SAFE.zip',
 
 ]
s3_images = {}
for i in sorted(glob.glob(os.path.join(s3_dir, 'images', '*', "*.SEN3"))):
	s3_images[os.path.basename(os.path.dirname(i))] = i

for s2_image in s2_images:
  date = os.path.basename(os.path.dirname(s2_image))
  s3_image = s3_images[date]
  # at this stage the climate data must be downloaded pefore running sen-et
  # and a template of the climate data must be made
  ecmwf_ERA5_dir = os.path.join(os.path.dirname(s3_image), 'ecmwf_ERA5')
  sen_et(gpt, py39, s2_dir, s2_image, s3_image, project_gdf, outputs_dir, lcc_tif, c_l, ecmwf_ERA5_dir)