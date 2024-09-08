from wapordl import wapor_map
import os

def download_WaPOR(region, variables, period, output_dr, overview = "NONE"):
  for var in variables:
    download_dr = os.path.join(output_dr, 'WaPOR_data', var)
    os.makedirs(download_dr, exist_ok=True)

    if('-E' in var):
      unit = "day"
    elif('-D' in var):
      unit = "dekad"
    elif('-M' in var):
      unit = "month"
    elif ('-A' in var):
      unit = "year"
    else:
      unit = "none"

    wapor_map(region, var, period, download_dr, seperate_unscale = True, unit_conversion = unit)