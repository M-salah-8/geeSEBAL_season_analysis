# SEBAL and geeSEBAL

## outputs:
* ETa = actual evapotranspiration
* NDVI = Normalized Difference Vegetation Index
* RGB = red,green,blue image

## main satellite
* Landsat 8 and 9

## What is SEBAL?
Surface Energy Balance Algorithm for Land (SEBAL) was developed and validated by Bastiaanssen (Bastiaanssen, 1995; Bastiaanssen et al., 1998a, 1998b) to 
estimate evapotranspiration (ET) from energy balance equation (Rn – G = LE + H), where LE, Rn, G and H are Latent Heat Flux, Net Radiation, Soil Heat Flux and Sensible Heat Flux, respectively.

# Google Earth Engine SEBAL Model
[geeSEBAL](https://github.com/gee-hydro/geeSEBAL) is an open-source implementation of Surface Energy Balance Algorithm for Land (SEBAL) using Google Earth Engine (GEE). geeSEBAL is available in both Javascript and Python API.

## Edits on the original code
1. Replaced landsat collection 1 with landsat collection 2.
2. add landsat 9 collection 2.
3. Modification to hot and cold pixels selection. (more on that later)
4. Add Crop coefficient (kc) data.
5. Add Functions for visualization and data extraction.

## How To Use geeSEBAL
* Fill inputs.csv
### Make (GEE) account
from [here](https://earthengine.google.com/).

### Folders layout
```bash
┌─── project_folder ───┐
│   ├── area
│   │   └── project.geojson
|   │   
│  ...
│   └── README.md
└─────────────────────┘
```

## References
 [Laipelt et al. (2021)] Long-term monitoring of evapotranspiration using the SEBAL algorithm and Google Earth Engine cloud computing. (https://doi.org/10.1016/j.isprsjprs.2021.05.018)

# SEBAL Model (local)
* Converted from the geeSEBAL model.
* Analyze locally downloaded Landsat data.
* Easy to customize, debug, and excellent for testing specific model sections.
* Landsat data can be processed as soon as level 2 products are available (3 days for Landsat 9 and 4-11 days for Landsat 8) when real-time weather data is provided. And with ERA5 weather data, it's only 5 days behind the present.

## How To Use the Local Model
* Fill inputs.csv
### Data acquisition:
1. Download landsat and digital elevation model (DEM) images for the area of interest, from [earth explorer](https://earthexplorer.usgs.gov/).
2. ERA5-hourly data will be downloaded automticaly, but if you choose to download it manually then, for each landsat image get 10 metre u wind component, 10 metre v wind component, 2 metre dewpoint temperature, 2 metre temperature and Surface solar radiation downwards for the day of the image and the day before it (48 hours in total). Use this [link](https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=form).

### Folders layout
```bash
┌─── project_folder ───┐
│   ├── area
│   │   └── project.geojson
│   ├── local_data
│   │   ├── landsat
│   │   │   └── images
|   │   │       ├── 2024-06-17
|   │   │       |   ├── C09_L2SP_173048_20240617_20240618_02_T1
|   │   │       |   └── era-5-hourly.grib
|   │   │       ├── 2024-07-11
|   │   │       |   ├── C08_L2SP_173048_20240711_20240719_02_T1
|   │   │       |   └── era-5-hourly.grib
|   │   │      ...
|   │   │      
|   │   └── strm_30m.tif
|   │   
│  ...
│   └── README.md
└─────────────────────┘
```