import cv2, os
import rasterio
import argparse
import numpy as np
import geopandas as gpd
from osgeo import gdal
from affine import Affine
from rasterio.features import shapes
from shapely.geometry import shape

from utils import *
# from geo.utils import numpy_array_to_raster
# from geogrid_autorift.util import numpy_array_to_raster

parser = argparse.ArgumentParser()
parser.add_argument('--out_tif', type=str, default="./optical.tif", help="Path to the output optical merged image")
parser.add_argument('--roi', type=str, default="[32.06, 32.77, 76.86, 77.82]", help="Region of interest")
parser.add_argument('--dem_path', type=str, required=True, help="directory in which DEM file needs to be saved")
parser.add_argument('--out_shapefile', type=str, default="./", help="Path to the output shapefile")
parser.add_argument('--landsat_files', nargs='+', help='List of landsat index files needed for merging', required=True)

args = parser.parse_args()


def stitch_images(images, dem_path, band):
    N = len(images)
    ds = gdal.Open(dem_path)
    dem = ds.GetRasterBand(1).ReadAsArray()
    h, w = dem.shape
    ds = None
    comb = np.zeros((h, w))
    region = (comb==1)
    for i in range(N):
        img = read_raster(images[i], band)
        img = cv2.resize(img, (w, h))
        valid = (img>-2000)
        new_region = (region|valid)
        valid = (np.logical_not(region)&(new_region))
        comb[valid] = img[valid]
        region = new_region
    
    return comb
    
    
# Assumption the entire region lies in the same UTM zone
def preprocess_optical(opt_name, img_paths, roi, dem_path):
    if len(img_paths)<1:
        print("No Optical Image provided. Exiting...")
        return
    
    ds = gdal.Open(dem_path)
    geo = ds.GetGeoTransform()
    ds = None
    bbox = get_bbox(roi)
    dem_out = os.path.join(*dem_path.split('/')[:-1])+'/dem_crop.tif'
    zone = round((180+roi[2])/6)
    if not os.path.exists(dem_out):
        os.system(f'gdalwarp -s_srs "EPSG:4326" -t_srs "+proj=utm +zone={zone} +datum=WGS84 +units=m +no_defs" -te {bbox[0]} {bbox[1]} {bbox[2]} {bbox[3]} -of GTIFF {dem_path} {dem_out}')

    ds = gdal.Open(dem_out)
    dem = ds.GetRasterBand(1).ReadAsArray()
    geo = ds.GetGeoTransform()
    projs = ds.GetProjection()
    ds = None
    D = horn_gradient(dem, geo)
    dem_slope = dem_path[:-4]+'_slope.tif'
    ds = numpy_array_to_raster(dem_slope, np.array([D]), projs, geo, nband=1)
    ds.FlushCache()
    ds = None
    
    crop_img = []
    for img_path in img_paths:
        out_img = img_path[:-4] + '_cropped.tif'
#         os.system
        if not os.path.exists(out_img):
            os.system(f'gdalwarp -tr {geo[1]} {geo[5]} -te {bbox[0]} {bbox[1]} {bbox[2]} {bbox[3]} {img_path} {out_img}')
        if os.path.exists(out_img):
            crop_img.append(out_img)
    
    if len(crop_img)==0:
        print("No cropped image generated")
        return 
        
    combined = []
    for i in range(2):
        comb = stitch_images(crop_img, dem_slope, band=(i+1))
        combined.append(comb)
    
    ds = gdal.Open(dem_slope)
    projs = ds.GetProjection()
    geo = ds.GetGeoTransform()
    ds = None
    
    ds = numpy_array_to_raster(opt_name, np.array(combined), projs, geo, nband=2)
    ds.FlushCache()
    ds = None
    

def generate_shapefile(optical_data, dem_slope_path, out_shp, crs):  
    dem_slope = read_raster(dem_slope_path)
    dem = dem_slope.astype('float32')/(dem_slope.max())
    
    ds = gdal.Open(optical_data)
    geo = ds.GetGeoTransform()
    nswir = ds.GetRasterBand(1).ReadAsArray()
    ndsi = ds.GetRasterBand(2).ReadAsArray()
    valid = (nswir>0)
    valid2 = (ndsi>0)
    nodata = (nswir<0)
    nodata2 = (ndsi<0)
    nswir[nodata] = np.median(nswir[valid])
    ndsi[nodata2] = 0.1
    ds = None

    # Segmentation
    img = nswir.astype('float32')/(nswir[valid].max())
    img2 = ndsi.astype('float32')/(ndsi[valid2].max())
    ret, thresh1 = cv2.threshold((img*255).astype('uint8'), 120, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
    ret2, thresh2 = cv2.threshold((img2*255).astype('uint8'), 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)

    im = ((img**0.1)*img2)
    ret, thresh3 = cv2.threshold((im*255).astype('uint8'), 50, 255, cv2.THRESH_BINARY)  # Hard coded
    th3 = cv2.adaptiveThreshold((im*255).astype('uint8'),255,cv2.ADAPTIVE_THRESH_GAUSSIAN_C, cv2.THRESH_BINARY,51,2)
    th3 = cv2.GaussianBlur(th3,(7,7),2)/255

    # Get thresh3 using DEM
    mas = ((th3*im)**0.5*thresh3)
        
    # Try 1
    # i2 = mas*(1-dem)**2
    # i2 = cv2.GaussianBlur(i2.astype('uint8'),(7,7),7)
    # ret, thresh5 = cv2.threshold((i2).astype('uint8'), 140, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)

    # Try 2
    ii = ((th3*im)**0.5)*(1-dem)   #mas*((1-dem1)**1)
    ret, i3 = cv2.threshold((ii*255).astype('uint8'), 120, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
    gaus_img = cv2.GaussianBlur(i3.astype('uint8'),(7,7),7)
    ret, mask = cv2.threshold((gaus_img).astype('uint8'), 140, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)

    # Removing small contours
    mask_clean = clean_mask(mask, area_thresh=250)
    cv2.imwrite('test.jpg', np.hstack([mask, mask_clean]))
    glaciers = []
    area = []
    affine = Affine(geo[1], geo[2], geo[0], geo[4], geo[5], geo[3])
    
    for shp, val in shapes(mask_clean.astype('float32'), transform=affine):
        if val>0:
            glaciers.append(shape(shp))
            area.append(shape(shp).area)

    # print(glaciers)
    gdf3 = gpd.GeoDataFrame(geometry=glaciers, crs=crs)  # Note GeoDataFrame geometry requires a list
    gdf3.to_file(filename=out_shp, driver='ESRI Shapefile')


if __name__=='__main__':
    roi = np.array(args.roi[1:-1].replace(' ','').split(',')).astype('float')
    zone = round((180+roi[2])/6)
    epsg = get_espg(roi[0], zone)
    crs = f'EPSG:{epsg}'
    dem_slope_path = args.dem_path[:-4]+'_slope.tif'
    print(args.landsat_files, crs)
    # img_paths = ['../2022_1/landsat2022_1.tif','../2022_2/landsat2022_2.tif']
    preprocess_optical(args.out_tif, args.landsat_files, roi, args.dem_path)
    generate_shapefile(args.out_tif, dem_slope_path, args.out_shapefile, crs)
