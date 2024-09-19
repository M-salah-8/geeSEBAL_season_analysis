import os
import subprocess
import xml.etree.ElementTree as ET

def s2_biophysical_processor(gpt, images_drs, output_dr, gdf):
	xml_file = os.path.join("products", "growth", "s2_bp", "BP_s2A_10m.xml")
	gdf = gdf.to_crs('epsg:4326')
	wkt_data = gdf.boundary.to_wkt().iloc[0]
	tree = ET.parse(xml_file)
	root = tree.getroot()
	geo_region = root.find('.//geoRegion')
	geo_region.text = wkt_data
	file_elements = root.findall('.//file')
	input_file = file_elements[0]
	output_file = file_elements[1]

	for image_dr in images_drs:
		date = os.path.basename(image_dr).split('_')[-1].split('.')[0][0:8]
		date = f"{date[0:4]}_{date[4:6]}_{date[6:8]}"
		os.makedirs(os.path.join(output_dr, "lai", date),exist_ok=True)
		tif_output_dr = os.path.join(output_dr, "lai", date, f"{date}.dim")
		input_file.text = image_dr
		output_file.text = tif_output_dr
		os.makedirs(os.path.join(output_dr, "graphs"), exist_ok=True)
		xml_dr = os.path.join(output_dr, "graphs", f"{date}.xml")
		tree.write(xml_dr)
		command = [gpt, xml_dr]
		subprocess.run(command)