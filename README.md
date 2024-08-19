geeSEBAL is a open-source implementation of Surface Energy Balance Algorithm for Land (SEBAL) using Google Earth Engine (GEE). geeSEBAL is available in both Javascript and Python API.\
\
This code is an edited version of geeSEBAL to achive the following objectives:
1. Replace outdated landsat collection 1 with landsat collection 2.
2. Add data and functions to perform additional analysis and visualization.

\
Original geeSEBAL: (https://github.com/gee-hydro/geeSEBAL)
## How to use Google Earth Engine?

You need an account in GEE (https://earthengine.google.com/).

### Python API

Using pip to install earthengine-api

pip install earthengine-api
Authenticate Earth Engine library
import ee; ee.Authenticate()
## What is SEBAL?

Surface Energy Balance Algorithm for Land (SEBAL) was developed and validated by Bastiaanssen (Bastiaanssen, 1995; Bastiaanssen et al., 1998a, 1998b) to 
estimate evapotranspiration (ET) from energy balance equation (Rn – G = LE + H), where LE, Rn, G and H are Latent Heat Flux, Net Radiation, Soil Heat Flux and Sensible Heat Flux, respectively.

## Edits on the original code
1. Replaced landsat collection 1 with landsat collection 2.
2. Add potential evapotranspiration data.
3. Add kc (actual evapotranspiration / potential evapotranspiration) data.
4. add Functions for visualization and data extraction.

## References
 [Laipelt et al. (2021)] Long-term monitoring of evapotranspiration using the SEBAL algorithm and Google Earth Engine cloud computing. (https://doi.org/10.1016/j.isprsjprs.2021.05.018)

## To Use the Local Model:
* Download landsat and dem images from earth explorer (https://earthexplorer.usgs.gov/)
* Download era5-hourly data for each landsat image (for each landsat image get era5-hourly data (10 metre u wind component, 10 metre v wind component, 2 metre dewpoint temperature, 2 metre temperature, Surface solar radiation downwards) for the day of the image and the day before it (48 hour)) (https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=form)
\
* Add "local data" folder
* Put dem image inside the folder and rename "strm_30m.tif".
* Make "landsat images" folder
* For every landsat images folder that represents one day, make a folder and put the landsat images folder longside with the era5 grib file
e.g.:\
local data/landsat images/1/LC08_L2SP_173048_20240711_20240719_02_T1/\
local data/landsat images/1/era5-hourly.grib	(rename era5 grib to "era5-hourly.grib")