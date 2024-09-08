import ee
from .image import Image_gee

def process_images(start_date, end_date, aoi, hc_limits = {'NDVI_cold': 1, 'Ts_cold': 1, 'NDVI_hot': 1, 'Ts_hot': 1}, ls8=True, ls9=True):
    result_imgs = []
    print('processing images')
    if ls8:
        ls8_imgs = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')\
            .filterBounds(aoi)\
            .filter(ee.Filter.date(start_date, end_date)).distinct(['DATE_ACQUIRED'])   ### check
        ls8_list_size = ls8_imgs.size().getInfo()
        print(f"{ls8_list_size} landsat 8 images found")
        ls8_list = ls8_imgs.toList(ls8_list_size)
        print('processing')
        for i in range(ls8_list_size):
            img = ee.Image_gee(ls8_list.get(i))
            landsat_img = Image_gee(img.get('system:id').getInfo(),
                            NDVI_cold= hc_limits['NDVI_cold'],
                            Ts_cold= hc_limits['Ts_cold'],
                            NDVI_hot= hc_limits['NDVI_hot'],
                            Ts_hot= hc_limits['Ts_hot'])
            result_imgs.append(landsat_img)
            print(f"{i+1} / {ls8_list_size} done")
        print('landsate 8 images done')
    
    if ls9:
        ls9_imgs = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2')\
            .filterBounds(aoi)\
            .filter(ee.Filter.date(start_date, end_date)).distinct(['DATE_ACQUIRED'])
        ls9_list_size = ls9_imgs.size().getInfo()
        print(f"{ls9_list_size} landsat 9 images found")
        ls9_list = ls9_imgs.toList(ls9_list_size)
        print('processing')
        for i in range(ls9_list_size):
            img = ee.Image_gee(ls9_list.get(i))
            landsat_img = Image_gee(img.get('system:id').getInfo(),
                            NDVI_cold= hc_limits['NDVI_cold'],
                            Ts_cold= hc_limits['Ts_cold'],
                            NDVI_hot= hc_limits['NDVI_hot'],
                            Ts_hot= hc_limits['Ts_hot'])
            result_imgs.append(landsat_img)
            print(f"{i+1} / {ls9_list_size} done")
        print('landsate 9 images done')

    return result_imgs