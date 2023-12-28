import os, glob
import argparse
import numpy as np
from osgeo import gdal
from math import radians, sin

from utils import read_raster, numpy_array_to_raster
# from geogrid_autorift.util import numpy_array_to_raster


parser = argparse.ArgumentParser()
parser.add_argument('--out_filename', type=str, default="optical.tif", help="Path to the output optical merged image")
parser.add_argument('--sentinel_dir', type=str, required=True, help="directory in which Landsat file are saved")
 
args = parser.parse_args()


def read_meta(metadata):
    # metadata = glob.glob(f'{filepath}/*_MTL.txt')[0]
#     '{filepath}/LC08_L1TP_147037_20161206_20200905_02_T1_MTL.txt'
    fh = open(metadata)

    # These two lists will contain the rescaling parameters.
    mult_term = []
    add_term = []
    for line in fh:
        # Read the file line-by-line looking for the reflectance transformation parameters
        if "REFLECTANCE_MULT_BAND_" in line:
            mult_term.append(float(line.split("=")[1].strip()))
        elif "REFLECTANCE_ADD_BAND_" in line:
            add_term.append(float(line.split("=")[1].strip()))
        elif "SUN_ELEVATION" in line:
            # We're also getting the sun elevation from the metadata. It has
            # to be converted to a float and radians.
            sun_elevation = radians(float(line.split("=")[1].strip()))
    fh.close()  # Be sure to close an open file
    return mult_term, add_term, sun_elevation


def convert2refl(arr, band, mult_term, add_term, sun_elevation):
    refl = (arr*mult_term[band-1] + add_term[band-1])/sin(sun_elevation)
    return refl


def generate_indices(sentinel_dir, out_filename, reflectance=False):
    
    green = read_raster(glob.glob(sentinel_dir+'/*_B03_20m.jp2')[0])
    nir = read_raster(glob.glob(sentinel_dir+'/*_B8A_20m.jp2')[0])
    swir = read_raster(glob.glob(sentinel_dir+'/*_B11_20m.jp2')[0])
    # qa = read_raster(glob.glob(sentinel_dir+'/*_QA_PIXEL.TIF')[0])
    # mult_term, add_term, sun_elevation = read_meta(glob.glob(f'{sentinel_dir}/*_MTL.txt')[0])
    # if reflectance:
    #     print("Applying reflectance")
    #     green = convert2refl(green, 3, mult_term, add_term, sun_elevation)
    #     nir_ref = convert2refl(nir, 5, mult_term, add_term, sun_elevation)
    #     swir_ref = convert2refl(swir, 6, mult_term, add_term, sun_elevation)

    nodata = (green==-32767)
    if reflectance:
        NDSI = (green-swir_ref)/(green+swir_ref)
        NDWI = (green-nir_ref)/(green+nir_ref)
        swir_ref = None
        nir_ref = None
        green = None
    else:
        NDSI = (green-swir)/(green+swir)
        NDWI = (green-nir)/(green+nir)
        
    NDSI[nodata] = -32767
    NDWI[nodata] = -32767
    # qa[nodata] = -32767

    NSWIR = nir/swir
    mask = np.isnan(NSWIR)
    NSWIR[nodata] = -32767
    NSWIR[mask] = -32767
    print(NSWIR.max())
    swir = None
    nir = None
    # RSWIR = red/swir
    # RSWIR[nodata] = -32767

    ds = gdal.Open(glob.glob(sentinel_dir+'/*_B03_20m.jp2')[0])
    projs = ds.GetProjection()
    geo = ds.GetGeoTransform()
    ds = None

    ds = numpy_array_to_raster(f'{sentinel_dir}/{out_filename}', np.array([NSWIR, NDSI, NDWI]), projs, geo, nband=3)
    ds.FlushCache()
    ds = None


if __name__=='__main__':
    print(args.sentinel_dir, args.out_filename)
    generate_indices(args.sentinel_dir, args.out_filename, False)
