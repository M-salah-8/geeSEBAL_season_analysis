import numpy as np
import json
import os

from pyTSEB import TSEB

import snappy_utils as su
# import snappy


def _estimate_param_value(landcover, lut, band): 
    param_value = np.ones(landcover.shape) + np.nan

    for lc_class in np.unique(landcover[~np.isnan(landcover)]):
        lc_pixels = np.where(landcover == lc_class)
        lc_index = lut['landcover_class'].index(lc_class)
        param_value[lc_pixels] = lut[band][lc_index]
    return param_value

with open(os.path.join(os.path.dirname(__file__), "structural_params_parameters.json")) as f:
    parameters = json.load(f)
landcover_map = parameters["landcover_map"]
lai_map = parameters["lai_map"]
fgv_map = parameters["fgv_map"]
landcover_band = parameters["landcover_band"]
lookup_table = parameters["lookup_table"]
produce_vh = parameters["produce_vh"]
produce_fc = parameters["produce_fc"]
produce_chwr = parameters["produce_chwr"]
produce_lw = parameters["produce_lw"]
produce_lid = parameters["produce_lid"]
produce_igbp = parameters["produce_igbp"]
output_file = parameters["output_file"]

# Read the required data
PARAMS = ['veg_height', 'lai_max', 'is_herbaceous', 'veg_fractional_cover',
            'veg_height_width_ratio', 'veg_leaf_width', 'veg_inclination_distribution',
            'igbp_classification'
            ]

landcover, geo_coding = su.read_snappy_product(landcover_map, landcover_band)
landcover = landcover.astype(np.float32)
lai = su.read_snappy_product(lai_map, 'lai')[0].astype(np.float32)
fg = su.read_snappy_product(fgv_map, 'frac_green')[0].astype(np.float32)
with open(lookup_table, 'r') as fp:
    lines = fp.readlines()
headers = lines[0].rstrip().split(';')
values = [x.rstrip().split(';') for x in lines[1:]]
lut = {key: [float(x[idx]) for x in values if len(x) == len(headers)]
        for idx, key in enumerate(headers)}

for param in PARAMS:
    if param not in lut.keys():
        print(f'Error: Missing {param} in the look-up table')
        raise ValueError(f'Missing {param} in the look-up table')


band_data = []
param_value = np.ones(landcover.shape, np.float32) + np.nan

if produce_vh:
    for lc_class in np.unique(landcover[~np.isnan(landcover)]):
        lc_pixels = np.where(landcover == lc_class)
        lc_index = lut["landcover_class"].index(lc_class)
        param_value[lc_pixels] = lut['veg_height'][lc_index]

        # Vegetation height in herbaceous vegetation depends on plant area index
        if lut["is_herbaceous"][lc_index] == 1:
            pai = lai / fg
            pai = pai[lc_pixels]
            param_value[lc_pixels] = \
                0.1 * param_value[lc_pixels] + 0.9 * param_value[lc_pixels] *\
                np.minimum((pai / lut['veg_height'][lc_index])**3.0, 1.0)
    band_data.append({'band_name': 'veg_height', 'band_data': param_value})

if produce_fc:
    band_name = 'veg_fractional_cover'
    param_value = _estimate_param_value(landcover, lut, band_name)
    band_data.append({'band_name': band_name, 'band_data': param_value})

if produce_chwr:
    band_name = 'veg_height_width_ratio'
    param_value = _estimate_param_value(landcover, lut, band_name)
    band_data.append({'band_name': band_name, 'band_data': param_value})

if produce_lw:
    band_name = 'veg_leaf_width'
    param_value = _estimate_param_value(landcover, lut, band_name)
    band_data.append({'band_name': band_name, 'band_data': param_value})

if produce_lid:
    band_name = 'veg_inclination_distribution'
    param_value = _estimate_param_value(landcover, lut, band_name)
    band_data.append({'band_name': band_name, 'band_data': param_value})

if produce_igbp:
    band_name = 'igbp_classification'
    param_value = _estimate_param_value(landcover, lut, band_name)
    band_data.append({'band_name': band_name, 'band_data': param_value})

su.write_snappy_product(output_file, band_data, 'landcoverParams', geo_coding)
