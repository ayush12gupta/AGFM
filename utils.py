import os
import subprocess, sys
from pyproj import Transformer
from osgeo import gdal, osr
import numpy as np
import datetime
import xml.etree.ElementTree as ET

from geogrid_autorift.util import numpy_array_to_raster


def execute(cmd):
    subprocess.check_call(cmd, shell=True, stdout=sys.stdout, stderr=subprocess.STDOUT)


def read_vals(fn, nodat=None, band=1):
    ds = gdal.Open(fn)
    disp = ds.GetRasterBand(band).ReadAsArray()
    nodata = np.isnan(disp)
    disp[nodata] = -32767
    if nodat is not None:
        nodata = (nodat|nodata)
        disp[nodat] = -32767
    ds = None
    return disp, nodata


def get_DT(date1, date2):
    date1, date2 = str(date1), str(date2)
    date1 = datetime.date(int(date1[:4]), int(date1[4:6]), int(date1[6:]))
    date2 = datetime.date(int(date2[:4]), int(date2[4:6]), int(date2[6:]))
    return (date2-date1).days


def get_deltaT(dates):
    deltaT = []
    for i in range(len(dates)-1):
        deltaT.append(get_DT(dates[i], dates[i+1]))

    return deltaT


def getazimuthAngle(filename):
    from zipfile import ZipFile
    print(filename)
    with ZipFile(filename, 'r') as zipObj:
        listOfiles = zipObj.namelist()
        for l in listOfiles:
            if l.endswith('xml') and l.split('/')[-2]=='annotation':
                f = zipObj.open(l)
                fl = str(f.read())
                azimuth_angle = float(fl.split('platformHeading>')[1][:-2])
                break
    return azimuth_angle


def horn_gradient(z, geo):
    """calculate x and y gradients according to Horn (1981) method"""

    nrows, ncols = z.shape

    z_00 = z[0:nrows-2, 0:ncols-2]
    z_10 = z[1:nrows-1, 0:ncols-2]
    z_20 = z[2:nrows, 0:ncols-2]

    z_02 = z[0:nrows-2, 2:ncols]
    z_12 = z[1:nrows-1, 2:ncols]
    z_22 = z[2:nrows, 2:ncols]
    
    dz_dx_tmp = (z_02+2.0*z_12+z_22-z_00-2.0*z_10-z_20)/8.0
    dz_dx = np.zeros((nrows, ncols))
    dz_dx[1:-1,1:-1] = dz_dx_tmp
    
    z_00 = z[0:nrows-2, 0:ncols-2]
    z_01 = z[0:nrows-2, 1:ncols-1]
    z_02 = z[0:nrows-2, 2:ncols]

    z_20 = z[2:nrows, 0:ncols-2]
    z_21 = z[2:nrows, 1:ncols-1]
    z_22 = z[2:nrows, 2:ncols]
    
    dz_dy_tmp = (z_20+2.0*z_21+z_22-z_00-2.0*z_01-z_02)/8.0	
    dz_dy = np.zeros((nrows, ncols))
    dz_dy[1:-1,1:-1] = dz_dy_tmp

    dz_dy = dz_dy/geo[1]
    dz_dx = dz_dx/geo[5]
    
    return (dz_dy, dz_dx)


def directional_slope(slopex, slopey, angle):
    rad = (angle*np.pi)/180
    dslope = (slopex*np.sin(rad)) + (slopey*np.cos(rad))
    return dslope


def get_espg(lat, zone):
    epsg_code = 32600
    epsg_code += int(zone)
    if (lat < 0): # South
        epsg_code += 100    
    return epsg_code


def generate_dem_products(dem_dir, bbox, config=None):

    bbox = np.array(bbox[1:-1].replace(' ','').split(',')).astype('float')
    zone = round((180+bbox[2])/6)
    if not os.path.exists('dem.tif'):
        execute(f'gdalwarp -s_srs "EPSG:4326" -t_srs "+proj=utm +zone={zone} +datum=WGS84 +units=m +no_defs" -of GTIFF {dem_dir} dem.tif')

    # Assuming the entire region lies in the same UTM zone
    epsg = get_espg(bbox[0], zone)
    transformer = Transformer.from_crs("epsg:4326",f"epsg:{epsg}")
    xl,yt = bbox[0], bbox[2]
    xr, yb = bbox[1], bbox[3]
    x2,y2 = transformer.transform(xl,yt)
    x3,y3 = transformer.transform(xr,yb)

    if not os.path.exists('dem_crop.tif'):
        execute(f'gdalwarp -te {x2} {y2} {x3} {y3} dem.tif dem_crop.tif')
    execute('rm dem.tif')

    demvel = gdal.Open('dem_crop.tif')
    dem = demvel.GetRasterBand(1).ReadAsArray().astype(float)
    projs = demvel.GetProjection()
    geo = demvel.GetGeoTransform()
    demvel = None

    # Get azimuth angle
    try: 
        tree = ET.parse('reference.xml')
        ref_dir = str(tree.getroot().findall(".//*[@name='safe']")[0].text[1:-1])
        azimuth_angle = getazimuthAngle(ref_dir[1:-1])
    except:
        azimuth_angle = config['azimuthAngle']

    dz_dy, dz_dx = horn_gradient(dem, geo)
    slopey = directional_slope(dz_dx, dz_dy, float(azimuth_angle))
    slopey = np.arctan(slopey)*(180.0)/np.pi
    slopex = directional_slope(dz_dx, dz_dy, (float(azimuth_angle)+90))
    slopex = np.arctan(slopex)*(180.0)/np.pi
    slope = np.sqrt(dz_dx**2 + dz_dy**2)
    D = np.arctan(slope)*(180.0)/np.pi

    # Saving slope along X and Y directions
    ds = numpy_array_to_raster('dem_x.tif', np.expand_dims(slopex, 0), projs, geo, nband=1)
    ds.FlushCache()
    ds = None 
    ds = numpy_array_to_raster('dem_y.tif', np.expand_dims(slopey, 0), projs, geo, nband=1)
    ds.FlushCache()
    ds = None 
    ds = numpy_array_to_raster('dem_slope.tif', np.expand_dims(D, 0), projs, geo, nband=1)
    ds.FlushCache()
    ds = None
