from osgeo import gdal, osr
import os 
import subprocess, sys
import numpy as np 
from util import execute, numpy_array_to_raster, cropRaster


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Output geo grid')
    parser.add_argument('-d', '--direc', dest='path', type=str, required=True,
            help='Directory path')
    parser.add_argument('-f', '--factor', dest='factor', type=int, default=None,
            help='Scalar Factor')
    parser.add_argument('-shp', '--glacier_shp', dest='shpfile', type=str, default=None,
            help='Glacier Shapefile path')
    args = parser.parse_args()
    
    # Set the resolution of both incidence and velocity
    if os.path.exists('velocity.tif') and not os.path.exists('velocity_res.tif'):
        execute('gdalwarp -tr 30 -30 velocity.tif velocity_res.tif')
    else:
        print(f"WARNING: File named {'velocity.tif'} doesn't exists")
        
    # Set the resolution of both incidence and velocity
    if os.path.exists('snr.tif'):
        execute('gdalwarp -tr 30 -30 snr.tif snr_res.tif')
    else:
        print(f"WARNING: File named {'snr.tif'} doesn't exists")
    
    if os.path.exists('incid_crop.tif'):
        execute('gdalwarp -tr 30 -30 incid_crop.tif incid_res.tif')
    else:
        print(f"WARNING: File named {'incid_crop.tif'} doesn't exists")
        
    # Clip incid to velocity raster extent
    data = gdal.Open('velocity_res.tif')
    geoTransform = data.GetGeoTransform()
    projs = data.GetProjection()
    velocity = data.ReadAsArray()/args.factor
    minx = geoTransform[0]
    maxy = geoTransform[3]
    maxx = minx + geoTransform[1] * data.RasterXSize
    miny = maxy + geoTransform[5] * data.RasterYSize
    data = None
    
    if args.factor is not None:
        ds = numpy_array_to_raster('velocity_res.tif', velocity, projs, geoTransform, nband=2)
        ds.FlushCache()
        ds = None
    
    if os.path.exists('incid_res.tif'):
        execute('gdal_translate -projwin ' + ' '.join([str(x) for x in [minx, maxy, maxx, miny]]) + ' -of GTiff incid_res.tif incid_ang.tif')
        execute('rm incid_res.tif incid_crop.tif')
    
    # Cropping the velocity map to glacier boundary
    if args.shpfile is not None:
        cropRaster(args.shpfile, 'velocity_res.tif', 'velocity_res_crop.tif')
        cropRaster(args.shpfile, 'snr_res.tif', 'snr_res_crop.tif')
        