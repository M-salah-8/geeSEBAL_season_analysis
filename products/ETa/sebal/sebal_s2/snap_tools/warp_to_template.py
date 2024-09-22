import tempfile
import os
import json

import gdal_utils as gu
import snappy_utils as su


with open(os.path.join(os.path.dirname(__file__), "warp_to_template_parameters.json")) as f:
    parameters = json.load(f)

source = parameters["source"]
template = parameters["template"]
output = parameters["output"]
resample_algorithm = parameters["resample_algorithm"]

# Save source and template to GeoTIFF becasue it will need to be read by GDAL
temp_file = tempfile.NamedTemporaryFile(suffix=".tif", delete=False)
temp_source_path = temp_file.name
temp_file.close()
su.copy_bands_to_file(source, temp_source_path)
temp_file = tempfile.NamedTemporaryFile(suffix=".tif", delete=False)
temp_template_path = temp_file.name
temp_file.close()
su.copy_bands_to_file(template, temp_template_path)

# Wrap the source based on tamplate
wraped = gu.resample_with_gdalwarp(temp_source_path, temp_template_path, resample_algorithm)

# Save with snappy
name, geo_coding = su.get_product_info(template)[0:2]
bands = su.get_bands_info(source)
for i, band in enumerate(bands):
    band['band_data'] = wraped.GetRasterBand(i+1).ReadAsArray()
su.write_snappy_product(output, bands, name, geo_coding)

# Clean up
try:
    os.remove(temp_source_path)
    os.remove(temp_template_path)
except Exception:
    pass

