# Growth Monitoring
Using Sentinel Application Platform (SNAP) command line feature to estimate leaf area index (LAI).

## outputs
* LAI = leaf area index
* NDVI = Normalized Difference Vegetation Index

## main satellite
* Sentinel 2

## How To Use
* Fill inputs.csv
### Data acquisition:
1. sentinel 2 images for the area of interest, from [sentinel copernicus browser](https://browser.dataspace.copernicus.eu/).


### Folders layout
```bash
┌─── project_folder ───┐
│   ├── area
│   │   └── project.geojson
│   ├── local_data
│   │   └── sentinel_2
|   │       ├── 2023_11_04
|   │       |   └── S2B_MSIL2A_20231104T083029_N0509_R021_T35QRA_20231104T101901.SAFE.zip
|   │       ├── 2023_11_14
|   │       |   └── S2B_MSIL2A_20231114T083119_N0509_R021_T35QRA_20231114T102315.SAFE.zip
|   │      ...
|   │   
│  ...
│   └── README.md
└─────────────────────┘
```

# References
NASA ARSET [Crop-Specific Time Series Analysis for Growth Monitoring](https://appliedsciences.nasa.gov/get-involved/training/english/arset-mapping-crops-and-their-biophysical-characteristics)