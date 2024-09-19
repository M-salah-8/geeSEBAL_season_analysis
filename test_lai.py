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

s2_images = [ 
 '/home/msalah/Desktop/m salah/crop monitoring/crop_monitoring_code/local_data/sentinel_2/images/2023_12_14/S2B_MSIL2A_20231214T083249_N0510_R021_T35QRA_20231214T111928.SAFE.zip',
 '/home/msalah/Desktop/m salah/crop monitoring/crop_monitoring_code/local_data/sentinel_2/images/2023_12_19/S2A_MSIL2A_20231219T083351_N0510_R021_T35QRA_20231219T113423.SAFE.zip',
 '/home/msalah/Desktop/m salah/crop monitoring/crop_monitoring_code/local_data/sentinel_2/images/2023_12_24/S2B_MSIL2A_20231224T083259_N0510_R021_T35QRA_20231224T100541.SAFE.zip',
 '/home/msalah/Desktop/m salah/crop monitoring/crop_monitoring_code/local_data/sentinel_2/images/2023_12_29/S2A_MSIL2A_20231229T083351_N0510_R021_T35QRA_20231229T112451.SAFE.zip',
 '/home/msalah/Desktop/m salah/crop monitoring/crop_monitoring_code/local_data/sentinel_2/images/2024_01_03/S2B_MSIL2A_20240103T083249_N0510_R021_T35QRA_20240103T100512.SAFE.zip',
 '/home/msalah/Desktop/m salah/crop monitoring/crop_monitoring_code/local_data/sentinel_2/images/2024_01_08/S2A_MSIL2A_20240108T083331_N0510_R021_T35QRA_20240108T114557.SAFE.zip',
 '/home/msalah/Desktop/m salah/crop monitoring/crop_monitoring_code/local_data/sentinel_2/images/2024_01_13/S2B_MSIL2A_20240113T083219_N0510_R021_T35QRA_20240113T102012.SAFE.zip',
 '/home/msalah/Desktop/m salah/crop monitoring/crop_monitoring_code/local_data/sentinel_2/images/2024_01_18/S2A_MSIL2A_20240118T083251_N0510_R021_T35QRA_20240118T114555.SAFE.zip',
 '/home/msalah/Desktop/m salah/crop monitoring/crop_monitoring_code/local_data/sentinel_2/images/2024_01_23/S2B_MSIL2A_20240123T083149_N0510_R021_T35QRA_20240123T100647.SAFE.zip',
 '/home/msalah/Desktop/m salah/crop monitoring/crop_monitoring_code/local_data/sentinel_2/images/2024_01_28/S2A_MSIL2A_20240128T083211_N0510_R021_T35QRA_20240128T113252.SAFE.zip',
 '/home/msalah/Desktop/m salah/crop monitoring/crop_monitoring_code/local_data/sentinel_2/images/2024_02_02/S2B_MSIL2A_20240202T083059_N0510_R021_T35QRA_20240202T103958.SAFE.zip',
 ]

from products import s2_biophysical_processor

gpt = "/home/msalah/esa-snap/bin/gpt"		# path to gpt executable
output_dr = os.path.join(working_dir, 'local_data', 'products_outputs', season, 'growth')
s2_biophysical_processor(gpt, s2_images, output_dr, project_gdf)