import cv2
import numpy as np
import scipy.optimize as opt
from osgeo import gdal, osr
from scipy.interpolate import RectBivariateSpline


def BSpline(z, factor):
    c, r = z.shape
    x = np.arange(0, r, 1)
    y = np.arange(0, c, 1)
    spline = RectBivariateSpline(x, y, z)
    
    new_x = np.arange(0, r - 1e-10, (1/factor))
    new_y = np.arange(0, c - 1e-10, (1/factor))
    z_interp = spline(new_x, new_y)

    return z_interp

def twoD_Gaussian(xy, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    x, y = xy
    xo = float(xo)
    yo = float(yo)
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo)
                            + c*((y-yo)**2)))
    return g.ravel()

def gaussian_fit(res, scale, fit=True):
    res2 = BSpline(res, scale)
    (_, maxVal, _, maxLoc) = cv2.minMaxLoc(res2)
    if not fit:
        return res2
    
    data = res2.ravel()
    
    h, w = res2.shape
    
    x = np.linspace(0, w, w)
    y = np.linspace(0, h, h)
    x, y = np.meshgrid(x, y)
    
    # initial guess of parameters
    initial_guess = []
    initial_guess.append(res2.max())
    initial_guess.extend(list(maxLoc))
    initial_guess.extend([w/4, h/4, 0, res2.mean()])

    # find the optimal Gaussian parameters
    popt, pcov = opt.curve_fit(twoD_Gaussian, (x, y), data, p0=initial_guess)
    
    data_fitted = twoD_Gaussian((x, y), *popt).reshape(h,w)    
    return data_fitted


def upscale(img, scale, type='spline'):
    scr = img.copy()
    rows, cols = map(int, img.shape)
    s = 1

    if scale>1:
        if type=='spline':
            scr = gaussian_fit(scr, scale, fit=False)
        elif type=='gaussian':
            scr = gaussian_fit(scr, scale, fit=False)
        elif type=='pyramid':
            while (s < scale):
                s *= 2
                scr = cv2.pyrUp(scr, dstsize=(s * cols, s * rows))
    
    return scr


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
    return np.array([(-np.sin(phi)*np.sin(incid)), (np.sin(incid)*np.cos(phi)), -np.cos(incid)])


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
