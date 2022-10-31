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
parser.add_argument('--dem_path', type=str, default=None, help="directory in which DEM file needs to be saved")
parser.add_argument('--out_shapefile', type=str, default="./", help="Path to the output shapefile")
parser.add_argument('--ref_shapefile', type=str, default=None, help="Path to the reference shapefile")
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
        print(img.shape)
        img = cv2.resize(img, (w, h))
        valid = (img>-2000)
        print(img)
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
    # print(bbox)
    dem_out = os.path.join(*dem_path.split('/')[:-1])+'/dem_crop.tif'
    zone = round((180+roi[2])/6)

    # Input DEM mush be in WGS Projection
    if os.path.exists(dem_out):
        os.remove(dem_out)

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
    print(geo, dem_out)
    crop_img = []
    for img_path in img_paths:
        out_img = img_path[:-4] + '_cropped.tif'
#         os.system
        print(f'gdalwarp -tr {geo[1]} {geo[5]} -te {bbox[0]} {bbox[1]} {bbox[2]} {bbox[3]} {img_path} {out_img}')
        if not os.path.exists(out_img):
            os.system(f'gdalwarp -tr {geo[1]} {geo[5]} -te {bbox[0]} {bbox[1]} {bbox[2]} {bbox[3]} {img_path} {out_img}')
        else:
            os.remove(out_img)
            os.system(f'gdalwarp -tr {geo[1]} {geo[5]} -te {bbox[0]} {bbox[1]} {bbox[2]} {bbox[3]} {img_path} {out_img}')
        if os.path.exists(out_img):
            crop_img.append(out_img)
    
    if len(crop_img)==0:
        print("No cropped image generated")
        return 
        
    print(len(crop_img))
    combined = []
    for i in range(3):   # 3 Bands
        comb = stitch_images(crop_img, dem_slope, band=(i+1))
        combined.append(comb)
    
    ds = gdal.Open(dem_slope)
    projs = ds.GetProjection()
    geo = ds.GetGeoTransform()
    ds = None
    
    print(geo)
    ds = numpy_array_to_raster(opt_name, np.array(combined), projs, geo, nband=3)
    ds.FlushCache()
    ds = None
    print("Pre-processing completed")


def generate_shapefile(optical_data, dem_slope_path, out_shp, ref_shp, crs):  
    dem_slope = read_raster(dem_slope_path)
    dem = dem_slope.astype('float32')/(dem_slope.max())
    dem = cv2.GaussianBlur(dem,(51,51),10)
    ds = gdal.Open(optical_data)
    geo = ds.GetGeoTransform()
    nswir = ds.GetRasterBand(1).ReadAsArray()
    ndsi = ds.GetRasterBand(2).ReadAsArray()
    ndwi = ds.GetRasterBand(3).ReadAsArray()
    valid = (nswir>0)
    valid2 = (ndsi>0)&(ndsi<=1)
    valid3 = (ndwi>0)&(ndwi<=1)
    nodata = (nswir<0)
    nodata2 = (ndsi<0)|(ndsi>1)
    nodata3 = (ndwi<0)|(ndwi>1)
    nswir[nodata] = np.median(nswir[valid])
    ndsi[nodata2] = 0.1
    ndwi[nodata3] = np.median(ndwi[valid3])
    ds = None

    # Segmentation
    img = nswir.astype('float32')/(nswir[valid].max())
    img2 = ndsi.astype('float32')/(ndsi[valid2].max())
    img3 = ndwi.astype('float32')/(ndwi[valid3].max())
    # ret, thresh1 = cv2.threshold((img*255).astype('uint8'), 120, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
    # ret2, thresh2 = cv2.threshold((img2*255).astype('uint8'), 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)

    im = ((img**0.1)*img2)
    ret, thresh3 = cv2.threshold((im*255).astype('uint8'), 50, 255, cv2.THRESH_BINARY)  # Hard coded
    th3 = cv2.adaptiveThreshold((im*255).astype('uint8'),255,cv2.ADAPTIVE_THRESH_GAUSSIAN_C, cv2.THRESH_BINARY,51,2)
    th3 = cv2.GaussianBlur(th3,(7,7),2)/255

    # Get thresh3 using DEM
    # mas = ((th3*im)**0.5*thresh3)

    ret, _ = cv2.threshold((img3*255).astype('uint8'), 100, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
    ret = max(ret, 100)
    im3 = (img3*255).astype('uint8')
    glacial_lake = clean_mask((im3>ret).astype('uint8')*255, 10)
    # mas[glacial_lake!=0] = 0
    # cl = clean_mask(mas.astype('uint8'))
        
    # Try 1
    # i2 = mas*(1-dem)**2
    # i2 = cv2.GaussianBlur(i2.astype('uint8'),(7,7),7)
    # ret, thresh5 = cv2.threshold((i2).astype('uint8'), 140, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)

    # Try 2
    # ii = ((th3*im)**0.5)*(1-dem)   #mas*((1-dem1)**1)
    # ret, i3 = cv2.threshold((ii*255).astype('uint8'), 120, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
    # # ret, i3 = cv2.threshold((ii*255).astype('uint8'), 120, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
    # # gaus_img = cv2.GaussianBlur(i3.astype('uint8'),(7,7),7)
    # gaus_img = cv2.GaussianBlur(i3.astype('uint8'),(5,5),0)
    # ret, mask = cv2.threshold((gaus_img).astype('uint8'), 140, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
    # mask[glacial_lake!=0] = 0

    # Try 3
    dem = (dem*255).astype('uint8')
    dem_bg = cv2.medianBlur(dem, 9) # Best 9
    imm = (im*255).astype('uint8')
    mask = (imm**(1-(dem_bg/255)**1.5)>40).astype('uint8')
    mask[glacial_lake!=0] = 0
    mask = clean_mask(mask, 10)

    # Removing small contours
    mask_clean = clean_mask(mask, area_thresh=0)
    glaciers = []
    area = []
    glims_id = []
    affine = Affine(geo[1], geo[2], geo[0], geo[4], geo[5], geo[3])
    unit_area = abs(geo[1]*geo[5])
    inter_thresh = unit_area*10
    
    # for shp, val in shapes(mask_clean.astype('float32'), transform=affine):
    #     if val>0:
    #         glaciers.append(shape(shp))
    #         area.append(shape(shp).area)

    if ref_shp is None:
        gips = None
    else:
        gips = gpd.read_file(ref_shp).to_crs(crs)
    print("Generating Shapefile .... ")
    for shp, val in shapes(mask_clean.astype('float32'), transform=affine):
        if val>0:
            if ref_shp is not None:
                gid = None
                max_ar = inter_thresh
                for j, row in gips.iterrows():
                    poly = shape(shp)
                    gip = row['geometry']
                    ar = poly.intersection(gip).area
                    if ar>max_ar:
                        gid = row['GLIMS_ID']
                        max_ar = ar
                if gid!=None:
                    glaciers.append(shape(shp))
                    area.append(shape(shp).area/1e6)
                    glims_id.append(gid)
            else:
                glaciers.append(shape(shp))
                area.append(shape(shp).area)

    if ref_shp is not None:
        d = {"GLIMS_ID":glims_id, "Area":area, "geometry": glaciers}
    else:
        d = {"Area":area, "geometry": glaciers}
    gdf3 = gpd.GeoDataFrame(d, crs=crs)  # Note GeoDataFrame geometry requires a list
    gdf3.to_file(filename=out_shp, driver='ESRI Shapefile')


if __name__=='__main__':
    roi = np.array(args.roi[1:-1].replace(' ','').split(',')).astype('float')
    zone = round((180+roi[2:].mean())/6)
    print(roi, zone)
    epsg = get_espg(roi[0], zone)
    crs = f'EPSG:{epsg}'
    # Checking for DEM file
    if args.dem_path is None:
        dem_path = os.path.join(*args.out_tif.split('/')[:-1],'dem_tmp.tif')
        # os.remove(dem_path)
    else:
        dem_path = args.dem_path
    
    if not os.path.exists(dem_path):
        print("Downloading DEM to", dem_path)
        download_DEM(roi, dem_path)
    
    dem_slope_path = dem_path[:-4]+'_slope.tif'
    print(args.landsat_files, crs)
    # img_paths = ['../2022_1/landsat2022_1.tif','../2022_2/landsat2022_2.tif']
    preprocess_optical(args.out_tif, args.landsat_files, roi, dem_path)
    generate_shapefile(args.out_tif, dem_slope_path, args.out_shapefile, args.ref_shapefile, crs)
