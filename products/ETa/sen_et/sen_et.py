import sys, os, glob
sys.path.append(os.getcwd())
import xml.etree.ElementTree as ET
import rasterio
import subprocess
from funcs import geo
import json
import shutil

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

def add_elevation(s2_dir, s2_file, calculations_dir, graphs_dir, gpt, wkt_data, graphs_output_dir):
    graph_dir = os.path.join(graphs_dir,"add_elevation.xml")
    tree = ET.parse(graph_dir)
    root = tree.getroot()
    file_elements = root.findall('.//file')
    input_file = file_elements[0]
    input_file.text = s2_file
    output_file = file_elements[1]
    output_file.text = os.path.join(s2_dir, os.path.basename(output_file.text))
    if os.path.exists(output_file.text):
        print("elevation already exists")
        elevation_files = glob.glob(os.path.join(s2_dir, "*!OUTPUT_SRTM_elevation!.*"))
        for file in elevation_files:
            if os.path.isfile(file):
                if not os.path.exists(os.path.join(calculations_dir, os.path.basename(file))):
                    shutil.copy(file, os.path.join(calculations_dir, os.path.basename(file)))
                else:
                    print("-----------elevation already exists in calculations----------")
            else:
                if not os.path.exists(os.path.join(calculations_dir, os.path.basename(file))):
                    shutil.copytree(file, os.path.join(calculations_dir, os.path.basename(file)))
                else:
                    print("-----------elevation already exists in calculations----------")
    else:
        geo_region = root.find('.//geoRegion')
        geo_region.text = wkt_data
        graph_output = os.path.join(graphs_output_dir, os.path.basename(graph_dir))
        tree.write(graph_output)
        command = [gpt, graph_output]
        subprocess.run(command)
        elevation_files = glob.glob(os.path.join(s2_dir, "*!OUTPUT_SRTM_elevation!.*"))
        for file in elevation_files:
            if os.path.isfile(file):
                shutil.copy(file, os.path.join(calculations_dir, os.path.basename(file)))
            else:
                shutil.copytree(file, os.path.join(calculations_dir, os.path.basename(file)))

def lc_preprocessing(lc_tif, c_l, graphs_dir, gpt, calculations_dir, graphs_output_dir):
    geo_in = geo()
    source = os.path.join(calculations_dir, "!OUTPUT_mask!.data", "mask.img")
    lc_array, meta = geo_in.coregister(source, lc_tif)
    lai_img = os.path.join(calculations_dir, "!OUTPUT_biophysical!.data", "lai.img")
    lai_array, _ = geo_in.coregister(source, lai_img)
    lc_array[lai_array > 0.1] = c_l
    tif_out_dir = os.path.join(lc_tif.split(".")[0] + "_co")
    geo_in.save_tif(lc_array, meta, tif_out_dir)
    graph_dir = os.path.join(graphs_dir,"tif_to_dim.xml")
    tree = ET.parse(graph_dir)
    root = tree.getroot()
    file_elements = root.findall('.//file')
    input_file = file_elements[0]
    input_file.text = tif_out_dir + ".tif"
    output_file = file_elements[1]
    output_file.text = os.path.join(calculations_dir, "!OUTPUT_CCI_landcover!.dim")
    band_new_name = root.find('.//name')
    band_new_name.text = "CCILandCover-2015"
    graph_output = os.path.join(graphs_output_dir, os.path.basename(graph_dir))
    tree.write(graph_output)
    command = [gpt, graph_output]
    subprocess.run(command)

def leaf_spectra(calculations_dir, sen_et_tools_dir, py39):
    biophysical_file = os.path.join(calculations_dir, '!OUTPUT_biophysical!.dim')
    output_file = os.path.join(calculations_dir, '!OUTPUT_LR&T!.dim')
    leaf_spectra = os.path.join(sen_et_tools_dir, "leaf_spectra.py")
    leaf_spectra_parameters = os.path.join(sen_et_tools_dir, "leaf_spectra_parameters.json")
    with open(leaf_spectra_parameters) as f:
        parameters = json.load(f)
    parameters['biophysical_file']= biophysical_file
    parameters['output_file']= output_file
    with open(leaf_spectra_parameters, 'w') as f:
        json.dump(parameters, f)
    command = [py39, leaf_spectra]
    subprocess.run(command)

def frac_green(calculations_dir, sen_et_tools_dir, py39, inputs):
    sza_file = os.path.join(calculations_dir, '!OUTPUT_sun_zenith_angle!.dim')
    biophysical_file = os.path.join(calculations_dir, '!OUTPUT_biophysical!.dim')
    min_frac_green = inputs['min_frac_green']
    output_file = os.path.join(calculations_dir, '!OUTPUT_FVG!.dim')
    frac_green = os.path.join(sen_et_tools_dir, "frac_green.py")
    frac_green_parameters = os.path.join(sen_et_tools_dir, "frac_green_parameters.json")
    with open(frac_green_parameters) as f:
        parameters = json.load(f)
    parameters['sza_file']= sza_file
    parameters['biophysical_file']= biophysical_file
    parameters['min_frac_green']= min_frac_green
    parameters['output_file']= output_file
    with open(frac_green_parameters, 'w') as f:
        json.dump(parameters, f)
    command = [py39, frac_green]
    subprocess.run(command)

def structural_params(calculations_dir, sen_et_tools_dir, py39):
    landcover_map = os.path.join(calculations_dir, "!OUTPUT_CCI_landcover!.dim")
    lai_map = os.path.join(calculations_dir, '!OUTPUT_biophysical!.dim')
    fgv_map = os.path.join(calculations_dir, '!OUTPUT_FVG!.dim')
    lookup_table = os.path.join(sen_et_tools_dir,"auxdata","LUT","WaPOR_CCI_LUT.csv")
    output_file = os.path.join(calculations_dir, '!OUTPUT_VSP!.dim')
    structural_params = os.path.join(sen_et_tools_dir, "structural_params.py")
    structural_params_parameters = os.path.join(sen_et_tools_dir, "structural_params_parameters.json")
    with open(structural_params_parameters) as f:
        parameters = json.load(f)
    parameters["landcover_map"] =  landcover_map
    parameters["lai_map"] =  lai_map
    parameters["fgv_map"] =  fgv_map
    parameters["landcover_band"] =  "CCILandCover-2015"
    parameters["lookup_table"] =  lookup_table
    parameters["produce_vh"] =  True
    parameters["produce_fc"] =  True
    parameters["produce_chwr"] =  True
    parameters["produce_lw"] =  True
    parameters["produce_lid"] =  True
    parameters["produce_igbp"] =  True
    parameters["output_file"] =  output_file
    with open(structural_params_parameters, 'w') as f:
        json.dump(parameters, f)
    command = [py39, structural_params]
    subprocess.run(command)

def aerodynamic_roughness(calculations_dir, sen_et_tools_dir, py39, inputs):
    lai_map = os.path.join(calculations_dir, '!OUTPUT_biophysical!.dim')
    landcover_params_map = os.path.join(calculations_dir, '!OUTPUT_VSP!.dim')
    soil_roughness = inputs['soil_roughness'] # default 0.01
    output_file = os.path.join(calculations_dir, '!OUTPUT_AeroRough!.dim')
    aerodynamic_roughness = os.path.join(sen_et_tools_dir, "aerodynamic_roughness.py")
    aerodynamic_roughness_parameters = os.path.join(sen_et_tools_dir, "aerodynamic_roughness_parameters.json")
    with open(aerodynamic_roughness_parameters) as f:
        parameters = json.load(f)
    parameters["lai_map"] = lai_map
    parameters["landcover_params_map"] = landcover_params_map
    parameters["soil_roughness"] = soil_roughness
    parameters["output_file"] = output_file
    with open(aerodynamic_roughness_parameters, 'w') as f:
        json.dump(parameters, f)
    command = [py39, aerodynamic_roughness]
    subprocess.run(command)

# sentinel 3
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

def warp_to_template(calculations_dir, sen_et_tools_dir, py39, inputs):
    source = os.path.join(calculations_dir, "!OUTPUT_observation_geometry!.dim")
    template = os.path.join(calculations_dir, "!OUTPUT_reflectance!.dim")
    resample_algorithm = inputs['resample_algorithm'] # default "cubicspline"
    output = os.path.join(calculations_dir, "!OUTPUT_observation_geometry_s2wrap!.dim")
    warp_to_template = os.path.join(sen_et_tools_dir, "warp_to_template.py")
    warp_to_template_parameters = os.path.join(sen_et_tools_dir, "warp_to_template_parameters.json")
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

def data_mining_sharpener(s3_file, calculations_dir, sen_et_tools_dir, py39, inputs):
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

    data_mining_sharpener = os.path.join(sen_et_tools_dir, "data_mining_sharpener.py")
    data_mining_sharpener_parameters = os.path.join(sen_et_tools_dir, "data_mining_sharpener_parameters.json")
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

# climate data
def download_ecmwf_data(s3_file, download_path, sen_et_tools_dir, py39, gdf):
    if os.path.exists(download_path):
        print("-------------------------ecmwf data already exists-------------------------")
    else:
        print("-------------------------downloading ecmwf data-------------------------")
        gdf = gdf.to_crs("epsg:4326")
        bbox = gdf.total_bounds
        area = int(bbox[3]+1), int(bbox[0]), int(bbox[1]), int(bbox[2]+1)
        dt = os.path.basename(s3_file).split("_")[-11]

        area = area
        date = dt[:4] + "-" + dt[4:6] + "-" + dt[6:8]
        download_temperature = True
        download_dewpoint = True
        download_pressure = True
        download_wind_speed = True
        download_clear_sky_solar_radiation = True
        download_solar_radiation = True
        overwrite = True

        os.makedirs(download_path, exist_ok=True)
        ecmwf_data_download = os.path.join(sen_et_tools_dir, "ecmwf_data_download.py")
        ecmwf_data_download_parameters = os.path.join(sen_et_tools_dir, "ecmwf_data_download_parameters.json")
        with open(ecmwf_data_download_parameters) as f:
            parameters = json.load(f)
        parameters["area"] = area
        parameters["date"] = date
        parameters["download_path"] = download_path
        parameters["download_temperature"] = download_temperature
        parameters["download_dewpoint"] = download_dewpoint
        parameters["download_pressure"] = download_pressure
        parameters["download_wind_speed"] = download_wind_speed
        parameters["download_clear_sky_solar_radiation"] = download_clear_sky_solar_radiation
        parameters["download_solar_radiation"] = download_solar_radiation
        parameters["overwrite"] = overwrite
        with open(ecmwf_data_download_parameters, 'w') as f:
            json.dump(parameters, f)
        command = [py39, ecmwf_data_download]
        subprocess.run(command)

def prepare_ecmwf_data(s3_file, calculations_dir, sen_et_tools_dir, py39, gdf, ecmwf_ERA5_dir, inputs):
    gdf = gdf.to_crs("EPSG:4326")
    center_lon = gdf.geometry.centroid.x.iloc[0]
    dt = os.path.basename(s3_file).split("_")[-11]
    dt = dt[:4] + "-" + dt[4:6] + "-" + dt[6:8] + " " + dt[9:11] + ":" + dt[11:13]
    elevation_map = os.path.join(calculations_dir, "!OUTPUT_SRTM_elevation!.dim")
    elevation_band = inputs["elevation_band"]
    ecmwf_data_file = os.path.join(ecmwf_ERA5_dir, "era5.nc")
    date_time_utc = dt
    time_zone = int(center_lon/15)
    print(time_zone)
    ecmwf_data_preparation = os.path.join(sen_et_tools_dir, "ecmwf_data_preparation.py")
    ecmwf_data_preparation_parameters = os.path.join(sen_et_tools_dir, "ecmwf_data_preparation_parameters.json")
    output_file = os.path.join(calculations_dir, "!OUTPUT_SurfaceMD!.dim")
    with open(ecmwf_data_preparation_parameters) as f:
        parameters = json.load(f)
    parameters["elevation_map"] = elevation_map
    parameters["elevation_band"] = elevation_band
    parameters["ecmwf_data_file"] = ecmwf_data_file
    parameters["date_time_utc"] = date_time_utc
    parameters["time_zone"] = time_zone
    parameters["prepare_temperature"] = True
    parameters["prepare_vapour_pressure"] = True
    parameters["prepare_air_pressure"] = True
    parameters["prepare_wind_speed"] = True
    parameters["prepare_clear_sky_solar_radiation"] = True
    parameters["prepare_daily_solar_irradiance"] = True
    parameters["output_file"] = output_file
    with open(ecmwf_data_preparation_parameters, 'w') as f:
        json.dump(parameters, f)
    command = [py39, ecmwf_data_preparation]
    subprocess.run(command)

# evapotranspiration estimation
def longwave_irradiance(calculations_dir, sen_et_tools_dir, py39, inputs):
    meteo_product = os.path.join(calculations_dir, "!OUTPUT_SurfaceMD!.dim")
    at_band = inputs["at_band"]
    vp_band = inputs["vp_band"]
    ap_band = inputs["ap_band"]
    at_height = inputs["at_height"]
    output_file = os.path.join(calculations_dir, "!OUTPUT_Atm_LI!.dim")
    longwave_irradiance = os.path.join(sen_et_tools_dir, "longwave_irradiance.py")
    longwave_irradiance_parameters = os.path.join(sen_et_tools_dir, "longwave_irradiance_parameters.json")
    with open(longwave_irradiance_parameters) as f:
        parameters = json.load(f)
    parameters["meteo_product"] = meteo_product
    parameters["at_band"] = at_band
    parameters["vp_band"] = vp_band
    parameters["ap_band"] = ap_band
    parameters["at_height"] = at_height
    parameters["output_file"] = output_file
    with open(longwave_irradiance_parameters, 'w') as f:
        json.dump(parameters, f)
    command = [py39, longwave_irradiance]
    subprocess.run(command)

def net_shortwave_radiation(calculations_dir, sen_et_tools_dir, py39, inputs):
    lsp_product = os.path.join(calculations_dir, "!OUTPUT_LR&T!.dim")
    lai_product = os.path.join(calculations_dir, "!OUTPUT_biophysical!.dim")
    csp_product = os.path.join(calculations_dir, "!OUTPUT_VSP!.dim")
    mi_product = os.path.join(calculations_dir, "!OUTPUT_SurfaceMD!.dim")
    sza_product = os.path.join(calculations_dir, "!OUTPUT_observation_geometry_s2wrap!.dim")
    soil_ref_vis = inputs["soil_ref_vis"]
    soil_ref_nir = inputs["soil_ref_nir"]
    output_file = os.path.join(calculations_dir, "!OUTPUT_NSR!.dim")
    net_shortwave_radiation = os.path.join(sen_et_tools_dir, "net_shortwave_radiation.py")
    net_shortwave_radiation_parameters = os.path.join(sen_et_tools_dir, "net_shortwave_radiation_parameters.json")
    with open(net_shortwave_radiation_parameters) as f:
        parameters = json.load(f)
    parameters["lsp_product"] = lsp_product
    parameters["lai_product"] = lai_product
    parameters["csp_product"] = csp_product
    parameters["mi_product"] = mi_product
    parameters["sza_product"] = sza_product
    parameters["soil_ref_vis"] = soil_ref_vis
    parameters["soil_ref_nir"] = soil_ref_nir
    parameters["output_file"] = output_file
    with open(net_shortwave_radiation_parameters, 'w') as f:
        json.dump(parameters, f)
    command = [py39, net_shortwave_radiation]
    subprocess.run(command)

def energy_fluxes(calculations_dir, sen_et_tools_dir, py39, inputs):
    lst = os.path.join(calculations_dir, "!OUTPUT_SharpLST!.dim")
    lst_vza = os.path.join(calculations_dir, "!OUTPUT_observation_geometry_s2wrap!.dim")
    lai = os.path.join(calculations_dir, "!OUTPUT_biophysical!.dim")
    csp = os.path.join(calculations_dir, "!OUTPUT_VSP!.dim")
    fgv = os.path.join(calculations_dir, "!OUTPUT_FVG!.dim")
    ar = os.path.join(calculations_dir, "!OUTPUT_AeroRough!.dim")
    mi = os.path.join(calculations_dir, "!OUTPUT_SurfaceMD!.dim")
    nsr = os.path.join(calculations_dir, "!OUTPUT_NSR!.dim")
    li = os.path.join(calculations_dir, "!OUTPUT_Atm_LI!.dim")
    mask = os.path.join(calculations_dir, "!OUTPUT_mask!.dim")
    soil_roughness = inputs["soil_roughness"]
    alpha_pt = inputs["alpha_pt"]
    atmospheric_measurement_height = inputs["atmospheric_measurement_height"]
    green_vegetation_emissivity = inputs["green_vegetation_emissivity"]
    soil_emissivity = inputs["soil_emissivity"]
    save_component_fluxes = True
    save_component_temperature = True
    save_aerodynamic_parameters = True
    output_file = os.path.join(calculations_dir, "!OUTPUT_LS_EnergyFluxes!.dim")
    energy_fluxes = os.path.join(sen_et_tools_dir, "energy_fluxes.py")
    energy_fluxes_parameters = os.path.join(sen_et_tools_dir, "energy_fluxes_parameters.json")
    with open(energy_fluxes_parameters) as f:
        parameters = json.load(f)
    parameters["lst"] = lst
    parameters["lst_vza"] = lst_vza
    parameters["lai"] = lai
    parameters["csp"] = csp
    parameters["fgv"] = fgv
    parameters["ar"] = ar
    parameters["mi"] = mi
    parameters["nsr"] = nsr
    parameters["li"] = li
    parameters["mask"] = mask
    parameters["soil_roughness"] = soil_roughness
    parameters["alpha_pt"] = alpha_pt
    parameters["atmospheric_measurement_height"] = atmospheric_measurement_height
    parameters["green_vegetation_emissivity"] = green_vegetation_emissivity
    parameters["soil_emissivity"] = soil_emissivity
    parameters["save_component_fluxes"] = save_component_fluxes
    parameters["save_component_temperature"] = save_component_temperature
    parameters["save_aerodynamic_parameters"] = save_aerodynamic_parameters
    parameters["output_file"] = output_file
    with open(energy_fluxes_parameters, 'w') as f:
        json.dump(parameters, f)
    command = [py39, energy_fluxes]
    subprocess.run(command)

def daily_evapotranspiration(s3_file, calculations_dir, output_dir, sen_et_tools_dir, py39):
    dt = os.path.basename(s3_file).split("_")[-11]
    date = dt[:4] + "_" + dt[4:6] + "_" + dt[6:8]
    ief_file =os.path.join(calculations_dir, "!OUTPUT_LS_EnergyFluxes!.dim")
    mi_file = os.path.join(calculations_dir, "!OUTPUT_SurfaceMD!.dim")
    os.makedirs(os.path.join(output_dir, "eta", date), exist_ok=True)
    output_file = os.path.join(output_dir, "eta", date, "!OUTPUT_eta!.dim")
    daily_evapotranspiration = os.path.join(sen_et_tools_dir, "daily_evapotranspiration.py")
    daily_evapotranspiration_parameters = os.path.join(sen_et_tools_dir, "daily_evapotranspiration_parameters.json")
    with open(daily_evapotranspiration_parameters) as f:
        parameters = json.load(f)
    parameters["ief_file"] = ief_file
    parameters["mi_file"] = mi_file
    parameters["output_file"] = output_file
    with open(daily_evapotranspiration_parameters, 'w') as f:
        json.dump(parameters, f)
    command = [py39, daily_evapotranspiration]
    subprocess.run(command)