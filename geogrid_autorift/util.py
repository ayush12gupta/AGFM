#!/usr/bin/env python3

from shutil import copyfile, move # Utilities for copying and moving files
from osgeo import gdal, osr
import os 
import subprocess, sys
import numpy as np  
import fiona
import netCDF4
import datetime
import rasterio
from rasterio.mask import mask
from pyproj import Transformer
import xml.etree.ElementTree as ET
# from geogrid_autorift.util import numpy_array_to_raster


# config ### HARD CODED
GDAL_DATA_TYPE = gdal.GDT_Float32
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


def read_velocity(path, xoff, yoff, xsize, ysize):
    ds = gdal.Open(path)
    V_R_asc = ds.GetRasterBand(1).ReadAsArray(xoff, yoff, xsize, ysize)
    V_A_asc = ds.GetRasterBand(2).ReadAsArray(xoff, yoff, xsize, ysize)
    ds = None
    return V_R_asc, V_A_asc


def get_flightang(path):
    with open(path, 'r') as fp:
        for l_no, line in enumerate(fp):
            # search string
            if 'Azimuth angle' in line:
                # print(float(line.split(":")[1][1:]))
                return float(line.split(":")[1][1:])


def get_incang(path):
    with open(path, 'r') as fp:
        for l_no, line in enumerate(fp):
            # search string
            if 'Incidence Angle' in line:
                return float(line.split(":")[1][1:])
            

def get_Va_mat(phi):
    return np.array([np.cos(phi), np.sin(phi), 0])


def get_Vr_mat(phi, incid):
    return np.array([(-np.sin(phi)*np.sin(incid)), (np.sin(incid)*np.cos(phi)), np.cos(incid)])


def get_matA_3D(phi_asc, incid_asc, phi_des, incid_des):
    A = np.zeros((4, 3))
    A[0] = get_Vr_mat(phi_asc, incid_asc)
    A[1] = get_Va_mat(phi_asc)
    A[2] = get_Vr_mat(phi_des, incid_des)
    A[3] = get_Va_mat(phi_des)
    return A
    
def get_matB(velR_asc, velA_asc, velR_des, velA_des):
    b = np.zeros((4, 1))
    b[0] = velR_asc
    b[1] = velA_asc
    b[2] = velR_des
    b[3] = velA_des
    return b


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
    demsp = gdal.Open('../dem_slope.tif')   
    projs = demsp.GetProjection()
    geo = demsp.GetGeoTransform()
    nodata = demsp.GetRasterBand(1).GetNoDataValue()
    angle_u = demsp.GetRasterBand(1).ReadAsArray()*deg2rad
    # angle_u = angle_u[1:,1:]
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
    h, w = angle_v.shape
    angle_u = angle_u[:h,:w]
    angle_u_nodt = angle_u_nodt[:h,:w]

    ## VF
    print(vr.shape, angle_u.shape)
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


def generate_dem_products(dem_dir, bbox, do_incid=True, config=None, stackfile=None):

    bbox = np.array(bbox[1:-1].replace(' ','').split(',')).astype('float')
    zone = round((180+bbox[2])/6)
    if not os.path.exists('dem_crop.tif'):
        execute(f'gdalwarp -s_srs "EPSG:4326" -t_srs "+proj=utm +zone={zone} +datum=WGS84 +units=m +no_defs" -of GTIFF {dem_dir} dem.tif')
    
    if os.path.exists('../merged/los.rdr.geo.vrt') and do_incid:
         execute(f'gdalwarp -s_srs "EPSG:4326" -t_srs "+proj=utm +zone={zone} +datum=WGS84 +units=m +no_defs" -of GTIFF ../merged/los.rdr.geo.vrt incid.tif')

    # Assuming the entire region lies in the same UTM zone
    epsg = get_espg(bbox[0], zone)
    transformer = Transformer.from_crs("epsg:4326",f"epsg:{epsg}")
    xl,yt = bbox[0], bbox[2]
    xr, yb = bbox[1], bbox[3]
    x2,y2 = transformer.transform(xl,yt)
    x3,y3 = transformer.transform(xr,yb)

    if not os.path.exists('dem_crop.tif'):
        execute(f'gdalwarp -te {x2} {y2} {x3} {y3} dem.tif dem_crop.tif')
    if os.path.exists('dem.tif'):
        execute('rm dem.tif')
    
    if os.path.exists('incid.tif') and not os.path.exists('incid_crop.tif'):
        execute(f'gdalwarp -te {x2} {y2} {x3} {y3} incid.tif incid_crop.tif')

    demvel = gdal.Open('dem_crop.tif')
    dem = demvel.GetRasterBand(1).ReadAsArray().astype(float)
    projs = demvel.GetProjection()
    geo = demvel.GetGeoTransform()
    demvel = None

    # Get azimuth angle
    if stackfile is None:
        try: 
            tree = ET.parse('../reference.xml')
            ref_dir = str(tree.getroot().findall(".//*[@name='safe']")[0].text[1:-1])
            azimuth_angle = getazimuthAngle(f'{ref_dir[1:-1]}')
        except:
            azimuth_angle = config['azimuthAngle']
    else:
        try:
            fl = open(stackfile)
            for l_no, line in enumerate(fl):
                if 'dirname' in line:
                    sar_dir = line.split(':')[-1][1:]
            fl.close()
            azimuth_angle = getazimuthAngle(sar_dir)
        except:
            azimuth_angle = 0

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
