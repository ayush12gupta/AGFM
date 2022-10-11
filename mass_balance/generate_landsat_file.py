import os, glob
import argparse
import numpy as np
from osgeo import gdal

from utils import read_raster, numpy_array_to_raster
# from geogrid_autorift.util import numpy_array_to_raster


parser = argparse.ArgumentParser()
parser.add_argument('--out_filename', type=str, default="optical.tif", help="Path to the output optical merged image")
parser.add_argument('--landsat_dir', type=str, required=True, help="directory in which Landsat file are saved")
 
args = parser.parse_args()


def generate_indices(landsat_dir, out_filename):
    
    print(landsat_dir+'/*_B3.TIF')
    green = read_raster(glob.glob(landsat_dir+'/*_B3.TIF')[0])
    # red = read_raster(glob.glob(landsat_dir+'/*_B4.TIF')[0])
    nir = read_raster(glob.glob(landsat_dir+'/*_B5.TIF')[0])
    swir = read_raster(glob.glob(landsat_dir+'/*_B6.TIF')[0])
    nodata = (green==-32767)

    NDSI = (green-swir)/(green+swir)
    NDWI = (green-nir)/(green+nir)
    NDSI[nodata] = -32767
    NDWI[nodata] = -32767

    NSWIR = nir/swir
    NSWIR[nodata] = -32767

    # RSWIR = red/swir
    # RSWIR[nodata] = -32767

    ds = gdal.Open(glob.glob(landsat_dir+'/*_B3.TIF')[0])
    projs = ds.GetProjection()
    geo = ds.GetGeoTransform()
    ds = None

    ds = numpy_array_to_raster(f'{landsat_dir}/{out_filename}', np.array([NSWIR, NDSI, NDWI]), projs, geo, nband=3)
    ds.FlushCache()
    ds = None


if __name__=='__main__':
    print(args.landsat_dir, args.out_filename)
    generate_indices(args.landsat_dir, args.out_filename)
