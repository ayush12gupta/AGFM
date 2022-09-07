import subprocess, sys
from pyproj import Transformer
from osgeo import gdal, osr
import numpy as np

from geogrid_autorift.util import numpy_array_to_raster


def execute(cmd):
    subprocess.check_call(cmd, shell=True, stdout=sys.stdout, stderr=subprocess.STDOUT)


def horn_gradient(z, geo):
    """calculate x and y gradients according to Horn (1981) method"""

    nrows, ncols = z.shape

    z_00 = z[0:nrows-2, 0:ncols-2]
    z_10 = z[1:nrows-1, 0:ncols-2]
    z_20 = z[2:nrows, 0:ncols-2]

    z_02 = z[0:nrows-2, 2:ncols]
    z_12 = z[1:nrows-1, 2:ncols]
    z_22 = z[2:nrows, 2:ncols]
    
    dz_dx = (z_02+2.0*z_12+z_22-z_00-2.0*z_10-z_20)/8.0
    
    z_00 = z[0:nrows-2, 0:ncols-2]
    z_01 = z[0:nrows-2, 1:ncols-1]
    z_02 = z[0:nrows-2, 2:ncols]

    z_20 = z[2:nrows, 0:ncols-2]
    z_21 = z[2:nrows, 1:ncols-1]
    z_22 = z[2:nrows, 2:ncols]
    
    dz_dy = (z_20+2.0*z_21+z_22-z_00-2.0*z_01-z_02)/8.0	

    dz_dy = dz_dy/geo[1]
    dz_dx = dz_dx/geo[5]
    
    return (dz_dy, dz_dx)


def generate_dem_products(dem_dir, bbox):

    bbox = np.array(bbox[1:-1].replace(' ','').split(',')).astype('float')
    zone = round((180+bbox[2])/6)
    execute(f'gdalwarp -s_srs "EPSG:4326" -t_srs "+proj=utm +zone={zone} +datum=WGS84 +units=m +no_defs" -of GTIFF {dem_dir} dem.tif')

    transformer = Transformer.from_crs("epsg:4326","epsg:32643")
    xl,yt = bbox[0], bbox[2]
    xr, yb = bbox[1], bbox[3]
    x2,y2 = transformer.transform(xl,yt)
    x3,y3 = transformer.transform(xr,yb)

    execute(f'gdalwarp -te {x2} {y2} {x3} {y3} dem.tif dem_crop.tif')
    execute('rm dem.tif')

    demvel = gdal.Open('dem_crop.tif')
    dem = demvel.GetRasterBand(1).ReadAsArray().astype(float)
    projs = demvel.GetProjection()
    geo = demvel.GetGeoTransform()
    demvel = None

    dz_dy, dz_dx = horn_gradient(dem, geo)

    # Saving slope along X and Y directions
    ds = numpy_array_to_raster('dem_x.tif', np.expand_dims(dz_dx, 0), projs, geo, nband=1)
    ds.FlushCache()
    ds = None 
    ds = numpy_array_to_raster('dem_y.tif', np.expand_dims(dz_dy, 0), projs, geo, nband=1)
    ds.FlushCache()
    ds = None 
