import rasterio
import os
from rasterio.warp import reproject, Resampling
import numpy as np

class GeoFunctions:
    def __init__(self):
        self.resampling_method = {
            'nearest': Resampling.nearest,
            'bilinear': Resampling.bilinear,
            'cubic': Resampling.cubic,
            'cubic_spline': Resampling.cubic_spline,
            'lanczos': Resampling.lanczos
        }

    def get_tif_props(self, tif):
        src = rasterio.open(tif)
        tif_props = {
            'crs': src.crs,
            'res': src.res,
            'nodata': src.nodata,
            'transform': src.transform,
            'shape': src.shape,
            'count': src.count,
            'meta': src.meta
        }
        return tif_props

    def save_tif(self, array, meta, name, output_dr = os.getcwd()):
        if len(array.shape) == 2:
            array = np.expand_dims(array, axis=0)
        with rasterio.open(os.path.join(output_dr, f'{name}.tif'), 'w', **meta) as dst:
            dst.write(array)

    def coregister(self, source, target, target_array = None, resampling = 'nearest', dtype = np.float32):
        # check if source is a dic
        if isinstance(source, dict):
            source_props = source
        elif isinstance(source, str):
            source_props = self.get_tif_props(source)
        else:
            raise TypeError('source must be a str or a dict')

        if isinstance(target, dict):
            target_props = target
            if target_array is None:
                raise ValueError('target_array must be provided if target is a dict')
        elif isinstance(target, str):
            target_props = self.get_tif_props(target)
            if target_array is None:
                target_array = rasterio.open(target).read()
        else:
            raise TypeError('target must be a str or a dict')
        
        array, transform = reproject(target_array,
            destination= np.zeros(source_props['shape'], dtype=dtype),
            src_transform= target_props['transform'],
            dst_transform=source_props['transform'],
            src_crs=target_props['crs'], 
            dst_crs=source_props['crs'],
            src_nodata=target_props['nodata'],
            dst_nodata=source_props['nodata'],
            dst_resolution=source_props['res'],
            resampling=self.resampling_method[resampling],
            )
        meta = target_props['meta']
        meta.update(
            width = array.shape[-1],
            height = array.shape[-2],
            transform = transform,
            dtype = dtype,
        )
        return array, meta

    def stack_images(self, tif_files):
        # get tifs properties from the first file
        source_tif = tif_files[0]
        source_props = self.get_tif_props(source_tif)
        count = source_props['count']
        meta = source_props['meta']
        meta.update(count= len(tif_files)*count)
        # stack images
        arrays_stack = np.zeros([len(tif_files)*count,meta['height'], meta['width']])
        for index,tif in enumerate(tif_files):
            array,_ = self.coregister(source_props, tif, 'nearest', dtype= meta['dtype'])
            bands_start = index*count
            bands_end = bands_start + count
            arrays_stack[bands_start:bands_end] = array

        return arrays_stack, meta
    
