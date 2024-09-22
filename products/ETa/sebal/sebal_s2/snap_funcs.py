import sys, os
sys.path.append(os.getcwd())
import xml.etree.ElementTree as ET
import rasterio
import subprocess
import numpy as np
from funcs import geo
import json

# sentinel 2
def s2_preprocessing(graphs_dir, s2_file, gpt, wkt_data, calculations_dir, graphs_output_dir):
    graph_dir = os.path.join(graphs_dir, "sentinel_2_pre_processing.xml")
    tree = ET.parse(graph_dir)
    root = tree.getroot()
    geo_region = root.findall('.//geoRegion')
    for i in geo_region:
        i.text = wkt_data
    file_elements = root.findall('.//file')
    input_file = file_elements[0]
    input_file.text = s2_file
    output_files = file_elements[1:]
    for output_file in output_files:
        output_file.text = os.path.join(calculations_dir, os.path.basename(output_file.text))
    graph_output = os.path.join(graphs_output_dir, os.path.basename(graph_dir))
    tree.write(graph_output)
    command = [gpt, graph_output]
    subprocess.run(command)
    reflectance = os.path.join(calculations_dir, "!OUTPUT_reflectance!.data")
    lai = os.path.join(calculations_dir, "!OUTPUT_biophysical!.data", "lai.img")
    s2_mask = os.path.join(calculations_dir, "!OUTPUT_mask!.data", "mask.img")
    return reflectance, s2_mask, lai

def add_elevation(strm_tif, calculations_dir, snap_tools, py39):
    source = os.path.join(calculations_dir, "!OUTPUT_mask!.data", "mask.img")
    geo_in = geo()
    array, meta = geo_in.coregister(source,strm_tif,dtype=np.int16,resampling='bilinear')
    tif_save_name = os.path.join(calculations_dir, "co_dem")
    geo_in.save_tif(array, meta, tif_save_name)
    co_dem = tif_save_name+'.tif'
    strm_to_product = os.path.join(snap_tools, "strm_to_product.py")
    strm_to_product_parameters = os.path.join(snap_tools, "strm_to_product_parameters.json")
    with open(strm_to_product_parameters) as f:
        parameters = json.load(f)
    parameters["strm_tif"] = co_dem
    parameters["calculations_dir"] = calculations_dir
    with open(strm_to_product_parameters, 'w') as f:
        json.dump(parameters, f)
    command = [py39, strm_to_product]
    subprocess.run(command)
    strm_output = os.path.join(calculations_dir, "!OUTPUT_SRTM_elevation!.dim", "elevation.img")
    return co_dem

def s3_preprocessing(graphs_dir, s3_file, gpt, wkt_data, calculations_dir, graphs_output_dir):
    source = os.path.join(calculations_dir, "!OUTPUT_mask!.data", "mask.img")
    graph_dir = os.path.join(graphs_dir, "sentinel_3_pre_processing.xml")
    tree = ET.parse(graph_dir)
    root = tree.getroot()
    geo_regions = root.findall('.//geoRegion')
    for geo_region in geo_regions:
        geo_region.text = wkt_data
    file_elements = root.findall('.//file')
    input_file = file_elements[0]
    input_file.text = s3_file
    output_files = file_elements[1:]
    for output_file in output_files:
        output_file.text = os.path.join(calculations_dir, os.path.basename(output_file.text))
    src = rasterio.open(source)
    crs_element = root.find('.//crs')
    crs_element.text = src.crs.to_wkt()
    graph_output = os.path.join(graphs_output_dir, os.path.basename(graph_dir))
    tree.write(graph_output)
    command = [gpt, graph_output]
    subprocess.run(command)
    lst = os.path.join(calculations_dir, "!OUTPUT_LST!.data", "LST.img")
    s3_mask = os.path.join(calculations_dir, "!OUTPUT_S3_mask!.data", "mask.img")
    return lst, s3_mask


def warp_to_template(calculations_dir, snap_tools, py39, inputs):
    source = os.path.join(calculations_dir, "!OUTPUT_observation_geometry!.dim")
    template = os.path.join(calculations_dir, "!OUTPUT_reflectance!.dim")
    resample_algorithm = inputs['resample_algorithm'] # default "cubicspline"
    output = os.path.join(calculations_dir, "!OUTPUT_observation_geometry_s2wrap!.dim")
    warp_to_template = os.path.join(snap_tools, "warp_to_template.py")
    warp_to_template_parameters = os.path.join(snap_tools, "warp_to_template_parameters.json")
    with open(warp_to_template_parameters) as f:
        parameters = json.load(f)
    parameters["source"] = source
    parameters["template"] = template
    parameters["resample_algorithm"] = resample_algorithm
    parameters["output"] = output
    with open(warp_to_template_parameters, 'w') as f:
        json.dump(parameters, f)
    command = [py39, warp_to_template]
    subprocess.run(command)

def data_mining_sharpener(s3_file, calculations_dir, snap_tools, py39, inputs):
    dt = os.path.basename(s3_file).split("_")[-11]
    dt = dt[:4] + "-" + dt[4:6] + "-" + dt[6:8] + " " + dt[9:11] + ":" + dt[11:13]
    sentinel_2_reflectance = os.path.join(calculations_dir, "!OUTPUT_reflectance!.dim")
    sentinel_3_lst = os.path.join(calculations_dir, "!OUTPUT_LST!.dim")
    high_res_dem = os.path.join(calculations_dir, "!OUTPUT_SRTM_elevation!.dim")
    high_res_geom = os.path.join(calculations_dir, "!OUTPUT_observation_geometry_s2wrap!.dim")
    lst_quality_mask = os.path.join(calculations_dir, "!OUTPUT_S3_mask!.dim")
    date_time_utc = dt
    elevation_band = inputs["elevation_band"]                      #default "elevation"
    lst_good_quality_flags = inputs["lst_good_quality_flags"]      #default "1"
    cv_homogeneity_threshold = inputs["cv_homogeneity_threshold"]  #default 0.0
    moving_window_size = inputs["moving_window_size"]              #default 3
    parallel_jobs = inputs["parallel_jobs"]                        #default 1
    output = os.path.join(calculations_dir, "!OUTPUT_SharpLST!.dim")

    data_mining_sharpener = os.path.join(snap_tools, "data_mining_sharpener.py")
    data_mining_sharpener_parameters = os.path.join(snap_tools, "data_mining_sharpener_parameters.json")
    with open(data_mining_sharpener_parameters) as f:
        parameters = json.load(f)
    parameters["sentinel_2_reflectance"] = sentinel_2_reflectance
    parameters["sentinel_3_lst"] = sentinel_3_lst
    parameters["high_res_dem"] = high_res_dem
    parameters["high_res_geom"] = high_res_geom
    parameters["lst_quality_mask"] = lst_quality_mask
    parameters["date_time_utc"] = date_time_utc
    parameters["elevation_band"] = elevation_band
    parameters["lst_good_quality_flags"] = lst_good_quality_flags
    parameters["cv_homogeneity_threshold"] = cv_homogeneity_threshold
    parameters["moving_window_size"] = moving_window_size
    parameters["parallel_jobs"] = parallel_jobs
    parameters["output"] = output
    with open(data_mining_sharpener_parameters, 'w') as f:
        json.dump(parameters, f)
    command = [py39, data_mining_sharpener]
    subprocess.run(command)
    output_image = os.path.join(output.replace('.dim', '.data'), 'sharpened_LST.img')
    return output_image