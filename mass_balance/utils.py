import cv2
import numpy as np
from osgeo import gdal
from pyproj import Transformer


def get_espg(lat, zone):
    epsg_code = 32600
    epsg_code += int(zone)
    if (lat < 0): # South
        epsg_code += 100    
    return epsg_code


def read_raster(fn, band=1):
    ds = gdal.Open(fn)
    disp = ds.GetRasterBand(band).ReadAsArray()
    nodat = ds.GetRasterBand(1).GetNoDataValue()
    nandata = np.isnan(disp)
    disp[disp==nodat] = -32767
    disp[nandata] = -32767
    ds = None
    return disp


def get_bbox(bbox):
    zone = round((180+bbox[2])/6)
    epsg = get_espg(bbox[0], zone)
    transformer = Transformer.from_crs("epsg:4326",f"epsg:{epsg}")
    xl,yt = bbox[0], bbox[2]
    xr, yb = bbox[1], bbox[3]
    x2,y2 = transformer.transform(xl,yt)
    x3,y3 = transformer.transform(xr,yb)
    return x2, y2, x3, y3


def horn_gradient(z, geo):
    """
    calculate x and y gradients according to Horn (1981) method
    """
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
    slope = np.sqrt(dz_dx**2 + dz_dy**2)
    D = np.arctan(slope)*(180.0)/np.pi
    
    return D


def clean_mask(mask, area_thresh=250):
    mask2 = mask.copy()
    mask2 -= 255
    _, contour, _ = cv2.findContours(mask2, cv2.RETR_TREE, cv2.CHAIN_APPROX_NONE)

    # get the biggest contour # returns _, contours, _ if using OpenCV 3
    area = []
    contours = []
    for con in contour:
        areac = cv2.contourArea(con)
        if areac>area_thresh:
            area.append(areac)
            contours.append(con)

    # fill in the contour
    cv2.drawContours(mask2, contours, -1, 255, -1)
    return mask2
    