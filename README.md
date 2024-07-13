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
