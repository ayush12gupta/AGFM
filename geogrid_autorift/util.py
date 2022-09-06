#!/usr/bin/env python3

from shutil import copyfile, move # Utilities for copying and moving files
from osgeo import gdal, osr
import os 
import numpy as np  
import fiona
import netCDF4
import rasterio
from rasterio.mask import mask


# config ### HARD CODED
GDAL_DATA_TYPE = gdal.GDT_Float64 
GEOTIFF_DRIVER_NAME = "GTiff"
NO_DATA = -32767
SPATIAL_REFERENCE_SYSTEM_WKID = 32767


def create_raster(output_path,
                  columns,
                  rows,
                  nband = 1,
                  gdal_data_type = GDAL_DATA_TYPE,
                  driver = GEOTIFF_DRIVER_NAME):

    driver = gdal.GetDriverByName(driver)

    output_raster = driver.Create(output_path, int(columns), int(rows), nband, eType = gdal_data_type)    
    return output_raster


def numpy_array_to_raster(output_path,
                          numpy_array,
                          proj,
                          trans,
                          nband = 1,
                          no_data = NO_DATA,
                          gdal_data_type = GDAL_DATA_TYPE,
                          spatial_reference_system_wkid = SPATIAL_REFERENCE_SYSTEM_WKID,
                          driver = GEOTIFF_DRIVER_NAME):

    shp = numpy_array.shape
    rows, columns = shp[1], shp[2]
    output_raster = create_raster(output_path, int(columns), int(rows), nband, gdal_data_type) 
    geotransform = trans

    spatial_reference = osr.SpatialReference()
    spatial_reference.ImportFromWkt(proj)
    output_raster.SetProjection(spatial_reference.ExportToWkt())
    output_raster.SetGeoTransform(geotransform)
    for i in range(1, nband+1):
        output_band = output_raster.GetRasterBand(i)
        output_band.SetNoDataValue(no_data)
        output_band.WriteArray(numpy_array[i-1])          
        output_band.FlushCache()
    
    return  output_raster


def cropRaster(shape_filenm, raster_filenm, cropped_filenm):
    
    # Read Shape file
    with fiona.open(shape_filenm, "r") as shapefile:
        shapes = [feature["geometry"] for feature in shapefile]

    # read imagery file
    with rasterio.open(raster_filenm) as src:
        out_image, out_transform = mask(src, shapes, crop=True)
        out_meta = src.meta

    # Save clipped imagery
    out_meta.update({"driver": "GTiff",
                     "height": out_image.shape[1],
                     "width": out_image.shape[2],
                     "transform": out_transform})

    with rasterio.open(cropped_filenm, "w", **out_meta) as dest:
        dest.write(out_image)


def generateCrops_NC(fn):

    file2read = netCDF4.Dataset(fn,'r')
    nm = fn.split('_')[0]
    #-----Params for velocity processing----------------
    shapefile_dir = '../data/' ## Hard Coded
    incidenceAngle = file2read['img_pair_info'].incidence_angle
    deg2rad = np.pi/180
    sin_incid = np.sin(incidenceAngle*deg2rad)
    cos_incid = np.cos(incidenceAngle*deg2rad)
    #----------- Angle u-----------------------
    demsp = gdal.Open('../data/demroi_2.tif')   ## Hard Coded
    projs = demsp.GetProjection()
    geo = demsp.GetGeoTransform()
    nodata = demsp.GetRasterBand(1).GetNoDataValue()
    angle_u = demsp.GetRasterBand(1).ReadAsArray()*deg2rad
    angle_u_nodt = (angle_u==nodata)
    demsp = None
    #----------- Compute angle V------------------------
    vr = np.array(file2read['vr']).astype(float)
    vr_nodt = (vr==-32767)
    vr[vr_nodt] = 0  #np.nan
    va = np.array(file2read['va']).astype(float)
    va_nodt = (va==-32767)
    va[va_nodt] = 0  #np.nan
    nodata = (va!=0)
    tmp = np.zeros_like(vr)
    tmp[nodata] = (vr*sin_incid)[nodata]/va[nodata]
    angle_v = np.arctan(tmp)
    h, w = angle_v.shape
    angle_u = angle_u[:h,:w]
    angle_u_nodt = angle_u_nodt[:h,:w]

    ## VF
    vf = vr*(((np.cos(angle_u)* np.sin(angle_v))*sin_incid) - (cos_incid*np.sin(angle_u))) + va*(np.cos(angle_u)*np.cos(angle_v))
    vf[angle_u_nodt] = 0
    vf[va_nodt] = -32767
    vf[vr_nodt] = -32767
    
    ## VN
    vn = vr*(((np.sin(angle_u)* np.sin(angle_v))*sin_incid) - (cos_incid*np.cos(angle_u))) + va*(np.sin(angle_u)*np.cos(angle_v))
    vn[angle_u_nodt] = 0
    vn[va_nodt] = -32767
    vn[vr_nodt] = -32767
    ds = numpy_array_to_raster('tmp_V.tif', np.array([vf, vn]), projs, geo, nband=2)
    ds.FlushCache()
    ds = None

    cropRaster(shapefile_dir + 'CBbasin.shp', 'tmp_V.tif', nm+'_V_basin.tiff')
    cropRaster(shapefile_dir + 'CBbasin_full.shp', 'tmp_V.tif', nm+'_V_full.tiff')
    os.remove('tmp_V.tif')


def generateCrops(fn, incidenceAngle, azimuthAngle, config):

    file2read = netCDF4.Dataset(fn,'r')
    nm = fn.split('_')[0]
    #-----Params for velocity processing----------------
    shapefile_dir = config['shapefile_dir']    
    # incidenceAngle = file2read['img_pair_info'].incidence_angle
    deg2rad = np.pi/180
    sin_incid = np.sin(incidenceAngle*deg2rad)
    cos_incid = np.cos(incidenceAngle*deg2rad)
    sin_az = np.sin(azimuthAngle*deg2rad)
    cos_az = np.cos(azimuthAngle*deg2rad)
    #----------- Angle u-----------------------
    demsp = gdal.Open(config['dem']['Slope'])   
    projs = demsp.GetProjection()
    geo = demsp.GetGeoTransform()
    nodata = demsp.GetRasterBand(1).GetNoDataValue()
    angle_u = demsp.GetRasterBand(1).ReadAsArray()*deg2rad
    angle_u_nodt = (angle_u==nodata)
    demsp = None
    #----------- Compute angle Va and Vr------------------------
    demvel = gdal.Open('velocity.tif')
    vx = demvel.GetRasterBand(1).ReadAsArray().astype(float)
    vx_nodt = (vx==-32767)
    vx[vx_nodt] = 0  #np.nan
    vy = demvel.GetRasterBand(2).ReadAsArray().astype(float)
    vy_nodt = (vy==-32767)
    vy[vy_nodt] = 0  #np.nan

    va = (vy*cos_az) - (vx*sin_az)
    vr = (vx*cos_az) + (vy*sin_az)

    nodata = (va!=0)
    tmp = np.zeros_like(vr)
    tmp[nodata] = (vr*sin_incid)[nodata]/va[nodata]
    angle_v = np.arctan(tmp)

    ## VF
    vf = vr*(((np.cos(angle_u)* np.sin(angle_v))*sin_incid) - (cos_incid*np.sin(angle_u))) + va*(np.cos(angle_u)*np.cos(angle_v))
    vf[angle_u_nodt] = 0
    vf[vx_nodt] = -32767
    vf[vy_nodt] = -32767
    ## VN
    vn = vr*(((np.sin(angle_u)* np.sin(angle_v))*sin_incid) - (cos_incid*np.cos(angle_u))) + va*(np.sin(angle_u)*np.cos(angle_v))
    vn[angle_u_nodt] = 0
    vn[vx_nodt] = -32767
    vn[vy_nodt] = -32767
    ds = numpy_array_to_raster('tmp_V.tif', np.array([vf, vn]), projs, geo, nband=2)
    ds.FlushCache()
    ds = None
    cropRaster(shapefile_dir + 'CBbasin.shp', 'tmp_V.tif', nm+'_V_basin.tiff')
    cropRaster(shapefile_dir + 'CBbasin_full.shp', 'tmp_V.tif', nm+'_V_full.tiff')
    os.remove('tmp_V.tif')
    os.chdir('..')


def postprocess(filename, kwargs, config):
    '''
    Custom postprocessing function
    '''
    generateCrops(filename, kwargs['incidence_angle'], kwargs['azimuth_angle'], config)