import numpy as np
import json
import os

import pyTSEB.net_radiation as rad
import snappy_utils as su


with open(os.path.join(os.path.dirname(__file__), "longwave_irradiance_parameters.json")) as f:
    parameters = json.load(f)
meteo_product = parameters["meteo_product"]
at_band = parameters["at_band"]
vp_band = parameters["vp_band"]
ap_band = parameters["ap_band"]
at_height = parameters["at_height"]
output_file = parameters["output_file"]

at, geo_coding = su.read_snappy_product(meteo_product, at_band)
at = at.astype(np.float32)
vp = su.read_snappy_product(meteo_product, vp_band)[0].astype(np.float32)
ap = su.read_snappy_product(meteo_product, ap_band)[0].astype(np.float32)

irrad = rad.calc_longwave_irradiance(vp, at, ap, at_height)

band_data = [
        {'band_name': 'longwave_irradiance', 'band_data': irrad}
]

su.write_snappy_product(output_file, band_data, 'longwaveIrradiance', geo_coding)