import numpy as np
import json
import os

import pyTSEB.net_radiation as rad
import pyTSEB.clumping_index as ci

import snappy_utils as su


with open(os.path.join(os.path.dirname(__file__), "net_shortwave_radiation_parameters.json")) as f:
    parameters = json.load(f)

lsp_product = parameters["lsp_product"]
lai_product = parameters["lai_product"]
csp_product = parameters["csp_product"]
mi_product = parameters["mi_product"]
sza_product = parameters["sza_product"]
soil_ref_vis = parameters["soil_ref_vis"]
soil_ref_nir = parameters["soil_ref_nir"]
output_file = parameters["output_file"]
    
refl_vis_c, geo_coding = su.read_snappy_product(lsp_product, 'refl_vis_c')
refl_vis_c = refl_vis_c.astype(np.float32)
refl_nir_c = su.read_snappy_product(lsp_product, 'refl_nir_c')[0].astype(np.float32)
trans_vis_c = su.read_snappy_product(lsp_product, 'trans_vis_c')[0].astype(np.float32)
trans_nir_c = su.read_snappy_product(lsp_product, 'trans_nir_c')[0].astype(np.float32)

lai = su.read_snappy_product(lai_product, 'lai')[0].astype(np.float32)


lad = su.read_snappy_product(csp_product, 'veg_inclination_distribution')[0].astype(np.float32)
frac_cover = su.read_snappy_product(csp_product, 'veg_fractional_cover')[0].astype(np.float32)
hw_ratio = su.read_snappy_product(csp_product, 'veg_height_width_ratio')[0].astype(np.float32)


p = su.read_snappy_product(mi_product, 'air_pressure')[0].astype(np.float32)
irradiance = su.read_snappy_product(mi_product, 'clear_sky_solar_radiation')[0].astype(np.float32)

sza = su.read_snappy_product(sza_product, 'solar_zenith_tn')[0].astype(np.float32)

net_rad_c = np.zeros(lai.shape, np.float32)
net_rad_s = np.zeros(lai.shape, np.float32)
soil_ref_vis = np.full(lai.shape, soil_ref_vis, np.float32)
soil_ref_nir = np.full(lai.shape, soil_ref_nir, np.float32)

#Estimate diffuse and direct irradiance
difvis, difnir, fvis, fnir = rad.calc_difuse_ratio(irradiance, sza, p)
skyl = difvis * fvis + difnir * fnir
irradiance_dir = irradiance * (1.0 - skyl)
irradiance_dif = irradiance * skyl

# Net shortwave radition for bare soil
i = lai <= 0
spectra_soil = fvis[i] * soil_ref_vis[i] + fnir[i] * soil_ref_nir[i]
net_rad_s[i] = (1. - spectra_soil) * (irradiance_dir[i] + irradiance_dif[i])

# Net shortwave radiation for vegetated areas
i = lai > 0
F = lai[i] / frac_cover[i] 
# Clumping index
omega0 = ci.calc_omega0_Kustas(lai[i], frac_cover[i], lad[i], isLAIeff=True)
omega = ci.calc_omega_Kustas(omega0, sza[i], hw_ratio[i])
lai_eff = F * omega
[net_rad_c[i], net_rad_s[i]] = rad.calc_Sn_Campbell(lai[i],
                                                    sza[i],
                                                    irradiance_dir[i],
                                                    irradiance_dif[i],
                                                    fvis[i],
                                                    fnir[i],
                                                    refl_vis_c[i],
                                                    trans_vis_c[i],
                                                    refl_nir_c[i],
                                                    trans_nir_c[i],
                                                    soil_ref_vis[i],
                                                    soil_ref_nir[i],
                                                    lad[i],
                                                    lai_eff
                                                    )


band_data = [
        {'band_name': 'net_shortwave_radiation_canopy', 'band_data': net_rad_c},
        {'band_name': 'net_shortwave_radiation_soil', 'band_data': net_rad_s}
]

su.write_snappy_product(output_file, band_data, 'netShortwaveRadiation', geo_coding)