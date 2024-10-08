import numpy as np
import json
import os

from pyTSEB import TSEB

import snappy_utils as su


with open(os.path.join(os.path.dirname(__file__), "energy_fluxes_parameters.json")) as f:
    parameters = json.load(f)

lst = parameters["lst"]
lst_vza = parameters["lst_vza"]
lai = parameters["lai"]
csp = parameters["csp"]
fgv = parameters["fgv"]
ar = parameters["ar"]
mi = parameters["mi"]
nsr = parameters["nsr"]
li = parameters["li"]
mask = parameters["mask"]
soil_roughness = parameters["soil_roughness"]
alpha_pt = parameters["alpha_pt"]
atmospheric_measurement_height = parameters["atmospheric_measurement_height"]
green_vegetation_emissivity = parameters["green_vegetation_emissivity"]
soil_emissivity = parameters["soil_emissivity"]
save_component_fluxes = parameters["save_component_fluxes"]
save_component_temperature = parameters["save_component_temperature"]
save_aerodynamic_parameters = parameters["save_aerodynamic_parameters"]
output_file = parameters["output_file"]

# Read the required data
lst = su.read_snappy_product(lst, 'sharpened_LST')[0].astype(np.float32)
vza = su.read_snappy_product(lst_vza, 'sat_zenith_tn')[0].astype(np.float32)
lai, geo_coding = su.read_snappy_product(lai, 'lai')
lai = lai.astype(np.float32)
lad = su.read_snappy_product(csp, 'veg_inclination_distribution')[0].astype(np.float32)
frac_cover = su.read_snappy_product(csp, 'veg_fractional_cover')[0].astype(np.float32)
h_w_ratio = su.read_snappy_product(csp, 'veg_height_width_ratio')[0].astype(np.float32)
leaf_width = su.read_snappy_product(csp, 'veg_leaf_width')[0].astype(np.float32)
veg_height = su.read_snappy_product(csp, 'veg_height')[0].astype(np.float32)
landcover_band = su.read_snappy_product(csp, 'igbp_classification')[0].astype(np.float32)
frac_green = su.read_snappy_product(fgv, 'frac_green')[0].astype(np.float32)
z_0M = su.read_snappy_product(ar, 'roughness_length')[0].astype(np.float32)
d_0 = su.read_snappy_product(ar, 'zero_plane_displacement')[0].astype(np.float32)
ta = su.read_snappy_product(mi, 'air_temperature')[0].astype(np.float32)
u = su.read_snappy_product(mi, 'wind_speed')[0].astype(np.float32)
ea = su.read_snappy_product(mi, 'vapour_pressure')[0].astype(np.float32)
p = su.read_snappy_product(mi, 'air_pressure')[0].astype(np.float32)
shortwave_rad_c = su.read_snappy_product(nsr, 'net_shortwave_radiation_canopy')[0].astype(np.float32)
shortwave_rad_s = su.read_snappy_product(nsr, 'net_shortwave_radiation_soil')[0].astype(np.float32)
longwave_irrad = su.read_snappy_product(li, 'longwave_irradiance')[0].astype(np.float32)
mask = su.read_snappy_product(mask, 'mask')[0].astype(np.float32)

# Model outputs
t_s = np.full(lai.shape, np.nan, np.float32)
t_c = np.full(lai.shape, np.nan, np.float32)
t_ac = np.full(lai.shape, np.nan, np.float32)
h_s = np.full(lai.shape, np.nan, np.float32)
h_c = np.full(lai.shape, np.nan, np.float32)
le_s = np.full(lai.shape, np.nan, np.float32)
le_c = np.full(lai.shape, np.nan, np.float32)
g = np.full(lai.shape, np.nan, np.float32)
ln_s = np.full(lai.shape, np.nan, np.float32)
ln_c = np.full(lai.shape, np.nan, np.float32)
r_s = np.full(lai.shape, np.nan, np.float32)
r_x = np.full(lai.shape, np.nan, np.float32)
r_a = np.full(lai.shape, np.nan, np.float32)
u_friction = np.full(lai.shape, np.nan, np.float32)
mol = np.full(lai.shape, np.nan, np.float32)
n_iterations = np.full(lai.shape, np.nan, np.float32)
flag = np.full(lai.shape, 255)
# ======================================
# First process bare soil cases
i = np.logical_and(np.isin(landcover_band, [13, 16]), mask == 1)
t_s[i] = lst[i]

# Calculate soil fluxes
[flag[i], ln_s[i], le_s[i], h_s[i], g[i], r_a[i], u_friction[i], mol[i],
n_iterations[i]] = TSEB.OSEB(lst[i],
                                ta[i],
                                u[i],
                                ea[i],
                                p[i],
                                shortwave_rad_s[i],
                                longwave_irrad[i],
                                soil_emissivity,
                                z_0M[i],
                                d_0[i],
                                atmospheric_measurement_height,
                                atmospheric_measurement_height,
                                calcG_params=[[1], 0.35])

# Set canopy fluxes to 0
ln_c[i] = 0.0
le_c[i] = 0.0
h_c[i] = 0.0
# ======================================
# Then process vegetated cases
i = np.logical_and(~np.isin(landcover_band, [13, 16, 17]), mask == 1)
# Emissivity of canopy containing green and non-green elements.
emissivity_veg = green_vegetation_emissivity * frac_green[i] + 0.91 * (1 - frac_green[i])
# Caculate component fluxes
[flag[i], t_s[i], t_c[i], t_ac[i], ln_s[i], ln_c[i], le_c[i], h_c[i], le_s[i], h_s[i],
g[i], r_s[i], r_x[i], r_a[i], u_friction[i], mol[i],
n_iterations[i]] = TSEB.TSEB_PT(lst[i],
                                vza[i],
                                ta[i],
                                u[i],
                                ea[i],
                                p[i],
                                shortwave_rad_c[i],
                                shortwave_rad_s[i],
                                longwave_irrad[i],
                                lai[i],
                                veg_height[i],
                                emissivity_veg,
                                soil_emissivity,
                                z_0M[i],
                                d_0[i],
                                atmospheric_measurement_height,
                                atmospheric_measurement_height,
                                f_c=frac_cover[i],
                                f_g=frac_green[i],
                                w_C=h_w_ratio[i],
                                leaf_width=leaf_width[i],
                                z0_soil=soil_roughness,
                                alpha_PT=alpha_pt,
                                x_LAD=lad[i],
                                calcG_params=[[1], 0.35],
                                resistance_form=[0, {}])

# Calculate the bulk fluxes
le = le_c + le_s
h = h_c + h_s
r_ns = shortwave_rad_c + shortwave_rad_s
r_nl = ln_c + ln_s
r_n = r_ns + r_nl
band_data = [
        {'band_name': 'sensible_heat_flux', 'band_data': h},
        {'band_name': 'latent_heat_flux', 'band_data': le},
        {'band_name': 'ground_heat_flux', 'band_data': g},
        {'band_name': 'net_radiation', 'band_data': r_n},
        {'band_name': 'quality_flag', 'band_data': flag}
        ]

if save_component_fluxes:
    band_data.extend(
            [
                {'band_name': 'sensible_heat_flux_canopy', 'band_data': h_c},
                {'band_name': 'sensible_heat_flux_soil', 'band_data': h_s},
                {'band_name': 'latent_heat_flux_canopy', 'band_data': le_c},
                {'band_name': 'latent_heat_flux_soil', 'band_data': le_s},
                {'band_name': 'net_longwave_radiation_canopy', 'band_data': ln_c},
                {'band_name': 'net_longwave_radiation_soil', 'band_data': ln_s}
            ]
    )
if save_component_temperature:
    band_data.extend(
            [
                {'band_name': 'temperature_canopy', 'band_data': t_c}, 
                {'band_name': 'temperature_soil', 'band_data': t_s},
                {'band_name': 'temperature_canopy_air', 'band_data': t_ac}
            ]
    )
if save_aerodynamic_parameters:
    band_data.extend(
            [ 
                {'band_name': 'resistance_surface', 'band_data': r_a}, 
                {'band_name': 'resistance_canopy', 'band_data': r_x},
                {'band_name': 'resistance_soil', 'band_data': r_s},
                {'band_name': 'friction_velocity', 'band_data': u_friction},
                {'band_name': 'monin_obukhov_length', 'band_data': mol}
            ]
    )

su.write_snappy_product(output_file, band_data, 'turbulentFluxes', geo_coding)