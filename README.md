## What is SEBAL?
Surface Energy Balance Algorithm for Land (SEBAL) was developed and validated by Bastiaanssen (Bastiaanssen, 1995; Bastiaanssen et al., 1998a, 1998b) to 
estimate evapotranspiration (ET) from energy balance equation (Rn – G = LE + H), where LE, Rn, G and H are Latent Heat Flux, Net Radiation, Soil Heat Flux and Sensible Heat Flux, respectively.

# Google Earth Engine SEBAL Model
geeSEBAL is an open-source implementation of Surface Energy Balance Algorithm for Land (SEBAL) using Google Earth Engine (GEE). geeSEBAL is available in both Javascript and Python API.\

Original geeSEBAL: (https://github.com/gee-hydro/geeSEBAL)

## Edits on the original code
1. Replaced landsat collection 1 with landsat collection 2.
2. add landsat 9 collection 2.
3. Modification to hot and cold pixels selection. (more on that later)
4. Add Crop coefficient (kc) data.
5. Add Functions for visualization and data extraction.
<!-- 3. Add potential evapotranspiration data. -->

## References
 [Laipelt et al. (2021)] Long-term monitoring of evapotranspiration using the SEBAL algorithm and Google Earth Engine cloud computing. (https://doi.org/10.1016/j.isprsjprs.2021.05.018)

# Local SEBAL Model
* Converted from the geeSEBAL model.
* Analyze locally downloaded Landsat data.
* Easy to customize, debug, and excellent for testing specific model sections.
* Landsat data can be processed as soon as level 2 products are available (3 days for Landsat 9 and 4-11 days for Landsat 8) when real-time weather data is provided. And with ERA5 weather data, it's only 5 days behind the present.

## How To Use the Local Model
#### Data acquisition:
1. Download landsat and digital elevation model (DEM) images from earth explorer. (https://earthexplorer.usgs.gov/)
2. Download era5-hourly data. For each landsat image get 10 metre u wind component, 10 metre v wind component, 2 metre dewpoint temperature, 2 metre temperature and Surface solar radiation downwards for the day of the image and the day before it (48 hours in total). (https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=form)
#### Folders structure:
1. in the main directory add "local data" folder.
2. Put the DEM image inside the folder and rename it "strm_30m.tif".\
local data/strm_30m.tif
3. Make "landsat images" folder.
4. For every landsat images folder _that represents one day _ make a folder and put the landsat images folder inside it, alongside with the era5 grib file.\
e.g.:\
local data/landsat images/1/LC08_L2SP_173048_20240711_20240719_02_T1/\
local data/landsat images/1/era5-hourly.grib	(rename era5 grib to "era5-hourly.grib")

### contact:
Mohammed Salah (mohammedalmak98@gmail.com)  