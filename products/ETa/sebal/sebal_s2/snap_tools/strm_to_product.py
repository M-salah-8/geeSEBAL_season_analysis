import os
import json
from snappy_utils import read_snappy_product, write_snappy_product


with open (os.path.join(os.path.dirname(__file__), "strm_to_product_parameters.json")) as f:
    parameters = json.load(f)

strm_tif = parameters["strm_tif"]
calculations_dir = parameters["calculations_dir"]
data, geo_coding = read_snappy_product(strm_tif)
write_snappy_product(os.path.join(calculations_dir, "!OUTPUT_SRTM_elevation!.dim"), [{'band_name': 'elevation', 'band_data': data}], '!OUTPUT_SRTM_elevation!', geo_coding)

