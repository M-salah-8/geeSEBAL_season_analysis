import os
import shutil
import json
from .sen_et import *

with open(os.path.join(os.path.dirname(__file__), "parameters.json")) as f:
    inputs = json.load(f)

def run_sen_et(gpt, py39, s2_dir, s2_file, s3_file, gdf, outputs_dir, lc_tif, c_l, ecmwf_ERA5_dir):
    sen_et_dir = os.path.dirname(__file__)
    graphs_dir = os.path.join(sen_et_dir, "sen-et_tools", "auxdata", "graphs")
    sen_et_tools_dir = os.path.join(sen_et_dir, "sen-et_tools")
    calculations_dir = os.path.join(os.path.dirname(s2_file), "calculations")
    os.makedirs(calculations_dir, exist_ok=True)
    graphs_output_dir = os.path.join(calculations_dir, "graphs")
    os.makedirs(graphs_output_dir, exist_ok=True)
    wkt_data = gdf.boundary.to_wkt().iloc[0]

    # sentinel 2
    s2_preprocessing(graphs_dir, s2_file, gpt, wkt_data, calculations_dir, graphs_output_dir)
    add_elevation(s2_dir, s2_file, calculations_dir, graphs_dir, gpt, wkt_data, graphs_output_dir)
    lc_preprocessing(lc_tif, c_l, graphs_dir, gpt, calculations_dir, graphs_output_dir)

    # estimate leaf reflectance and transmittance
    leaf_spectra(calculations_dir, sen_et_tools_dir, py39)

    # estimation fraction of vegetation which is green
    frac_green(calculations_dir, sen_et_tools_dir, py39, inputs["frac_green"])

    # produce maps of vegetation structural parameters
    structural_params(calculations_dir, sen_et_tools_dir, py39)

    # estimate aerodynamic roughness
    aerodynamic_roughness(calculations_dir, sen_et_tools_dir, py39, inputs["aerodynamic_roughness"])
    print("-------------------------end of s2-------------------------")

    # sentinel 3
    s3_preprocessing(graphs_dir, s3_file, gpt, wkt_data, calculations_dir, graphs_output_dir)

    # warp to template projection, resolution, and extent
    warp_to_template(calculations_dir, sen_et_tools_dir, py39, inputs["warp_to_template"])

    # sharpen LST with data mining sharpener
    data_mining_sharpener(s3_file, calculations_dir, sen_et_tools_dir, py39, inputs["data_mining_sharpener"])
    print("-------------------------end of s3-------------------------")

    # download ECMWF data
    download_ecmwf_data(s3_file, ecmwf_ERA5_dir, sen_et_tools_dir, py39, gdf)

    # prepare ERA5 data
    prepare_ecmwf_data(s3_file, calculations_dir, sen_et_tools_dir, py39, gdf, ecmwf_ERA5_dir, inputs["prepare_ecmwf_data"])
    print("-------------------------end of ecmwf-------------------------")

    # estimate atmosphere longwave irradiance
    longwave_irradiance(calculations_dir, sen_et_tools_dir, py39, inputs["longwave_irradiance"])

    # estimate net shortwave radiation
    net_shortwave_radiation(calculations_dir, sen_et_tools_dir, py39, inputs["net_shortwave_radiation"])

    # estimate land surface energy fluxes
    energy_fluxes(calculations_dir, sen_et_tools_dir, py39, inputs["energy_fluxes"])
    print("-------------------------end of energy fluxes-------------------------")

    # estimate daily evapotranspiration
    daily_evapotranspiration(s3_file, calculations_dir, outputs_dir, sen_et_tools_dir, py39)

    # clean up
    shutil.rmtree(calculations_dir)
