import os
import shutil
import json
from .snap_funcs import *

with open(os.path.join(os.path.dirname(__file__), "s2_parameters.json")) as f:
    inputs = json.load(f)

def run_sen_et(gpt, py39, s2_file, s3_file, calculations_dir, strm_tif, gdf):
    graphs_dir = os.path.join(os.path.dirname(__file__), "snap_graphs")
    snap_tools = os.path.join(os.path.dirname(__file__), "snap_tools")
    graphs_output_dir = os.path.join(calculations_dir, "graphs")
    os.makedirs(graphs_output_dir, exist_ok=True)
    wkt_data = gdf.boundary.to_wkt().iloc[0]

    # sentinel 2
    reflectance, s2_mask, lai = s2_preprocessing(graphs_dir, s2_file, gpt, wkt_data, calculations_dir, graphs_output_dir)
    co_dem = add_elevation(strm_tif, calculations_dir, snap_tools, py39)
    print("-------------------------end of s2-------------------------")

    # sentinel 3
    lst, s3_mask = s3_preprocessing(graphs_dir, s3_file, gpt, wkt_data, calculations_dir, graphs_output_dir)

    # warp to template projection, resolution, and extent
    warp_to_template(calculations_dir, snap_tools, py39, inputs["warp_to_template"])

    # sharpen LST with data mining sharpener
    sharpened_LST = data_mining_sharpener(s3_file, calculations_dir, snap_tools, py39, inputs["data_mining_sharpener"])
    print("-------------------------end of s3-------------------------")

    # clean up
    # shutil.rmtree(calculations_dir)

    return reflectance, s2_mask, lai, lst, s3_mask, sharpened_LST, co_dem
