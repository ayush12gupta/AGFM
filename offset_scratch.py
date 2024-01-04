import isce
import cv2, os
import subprocess
import numpy as np
from osgeo import gdal
from imageMath import IML
import matplotlib.pyplot as plt
from scipy.interpolate import RectBivariateSpline
import scipy.optimize as opt


NODATA = -32767

def read_img(filename):
    IMG = IML.mmapFromISCE(filename)
    img = IMG.bands[0]
    return np.absolute(img)

def read_raster(filename):
    ds = gdal.Open(filename, gdal.GA_ReadOnly)
    projs = ds.GetProjection()
    geo = ds.GetGeoTransform()
    v = ds.ReadAsArray()
    ds = None
    return v, projs, geo

def runCmd(cmd):
    out = subprocess.getoutput(cmd)
    return out

def read_config(txt_file):
    chipsizex0 = float(str.split(runCmd(f'fgrep "Smallest Allowable Chip Size in m:" {txt_file}'))[-1])
    gridspacingx = float(str.split(runCmd(f'fgrep "Grid spacing in m:" {txt_file}'))[-1])
#     pixsizex = float(str.split(runCmd(f'fgrep "Ground range pixel size:" {txt_file}'))[-1])
    
    incidenceAngle = float(str.split(runCmd(f'fgrep "Incidence Angle:" {txt_file}'))[-1])
    azimuthAngle = float(str.split(runCmd(f'fgrep "Azimuth angle:" {txt_file}'))[-1])
    rangePixelSize = float(str.split(runCmd(f'fgrep "Ground range pixel size:" {txt_file}'))[4])
    azimuthPixelSize = float(str.split(runCmd(f'fgrep "Azimuth pixel size:" {txt_file}'))[3])
    dt = float(str.split(runCmd(f'fgrep "Repeat Time:" {txt_file}'))[2])
    epsg = float(str.split(runCmd(f'fgrep "EPSG:" {txt_file}'))[1])
    kwargs = {'incidence_angle': incidenceAngle, 'azimuth_angle': azimuthAngle, 'rangePixelSize': rangePixelSize, 
              'azimuthPixelSize': azimuthPixelSize, 'gridspace': gridspacingx, 'chipsize0': chipsizex0, 'repeat_time': dt}
    return kwargs


class OffsetTracking:
    def __init__(self, I1_path, I2_path, po2v_path):

        txt_file_asc = f'{po2v_path}/testGeogrid.txt'
        config = read_config(txt_file_asc)
        
        self.chipsizex0 = config['chipsize0']
        self.pixsizex = config['rangePixelSize']
        self.pixsizey = config['azimuthPixelSize']
        self.ChipSize0X = int(np.ceil(self.chipsizex0/self.pixsizex/4)*4)
        self.ChipSize0Y = int(np.ceil(self.chipsizex0/self.pixsizey/4)*4)
#         self.ScaleChipSizeY = self.ChipSize0Y/self.ChipSize0X
        self.ScaleChipSizeY = 1

        self.sparseSearchSampleRate = 4
        self.chipmin = 240
        self.chipmax = 1920
        self.ChipSizeMaxX = (self.chipmax / self.chipsizex0) * self.ChipSize0X
        self.ChipSizeMinX = (self.chipmin / self.chipsizex0) * self.ChipSize0X
        self.OverSampleRatio = {self.ChipSize0X:32, self.ChipSize0X*2:64,self.ChipSize0X*4:128, self.ChipSize0X*8:128}
        ChipSizeUniX = np.unique(np.append(np.unique(self.ChipSizeMinX), np.unique(self.ChipSizeMaxX)))
        self.ChipSizeUniX = np.delete(ChipSizeUniX,np.where(ChipSizeUniX == 0)[0])
        ChipRangeX = self.ChipSize0X * np.array([1,2,4,8,16,32,64],np.float32)
        self.ChipSizeUniX = ChipRangeX[(ChipRangeX >= np.min(ChipSizeUniX)) & (ChipRangeX <= np.max(ChipSizeUniX))]
        self.maxScale = np.max(ChipSizeUniX) / self.ChipSize0X
        self.ChipSize0_GridSpacing_oversample_ratio = int(self.ChipSize0X / config['gridspace'])

        # Reading files
        gridfile_asc = f'{po2v_path}/window_location.tif'
        self.pos_asc, _, _ = read_raster(gridfile_asc)
        self.I1 = read_img(I1_path)
        self.I2 = read_img(I2_path)
        mask, _, _ = read_raster('/home/pc/DATA/mask.tif')
        self.mask = np.logical_not(mask[0]==-32767)
        self.po2vr, _, _ = read_raster(f'{po2v_path}/window_rdr_off2vel_x_vec.tif')
        self.po2va, _, _ = read_raster(f'{po2v_path}/window_rdr_off2vel_y_vec.tif')
        self.po2vr = self.po2vr[2]
        self.po2va = self.po2va[2]
        
        self.xGrid = self.pos_asc[0] #.astype(np.int32)
        self.yGrid = self.pos_asc[1] #.astype(np.int32)
        self.noDataMask = (self.xGrid == NODATA)
        # generate the nodata mask where offset searching will be skipped based on 1) imported nodata mask and/or 2) zero values in the image
        for ii in range(self.xGrid.shape[0]):
            for jj in range(self.xGrid.shape[1]):
                if (self.yGrid[ii,jj]!= NODATA)&(self.xGrid[ii,jj] != NODATA):
                    if (self.I1[self.yGrid[ii,jj]-1,self.xGrid[ii,jj]-1]==0)|(self.I2[self.yGrid[ii,jj]-1,self.xGrid[ii,jj]-1]==0):
                        self.noDataMask[ii,jj] = True

        self.SearchLimitX = np.ones(self.xGrid.shape, np.float32) * np.round(25)
        self.SearchLimitY = np.ones(self.xGrid.shape, np.float32) * np.round(25)
        
        # Setting values as zero
        self.xGrid[self.noDataMask] = 0
        self.yGrid[self.noDataMask] = 0
        self.SearchLimitX[self.noDataMask] = 0
        self.SearchLimitY[self.noDataMask] = 0
        self.xGrid = self.xGrid.astype(np.float32)
        self.yGrid = self.yGrid.astype(np.float32)
        
        # self.ChipSizeMinX = np.ones(self.xGrid.shape, np.float32) * np.round(self.ChipSizeMinX)
        # self.ChipSizeMaxX = np.ones(self.xGrid.shape, np.float32) * np.round(self.ChipSizeMaxX)
        # self.ChipSizeX = np.zeros(self.xGrid.shape, np.float32)
        # self.InterpMask = np.zeros(self.xGrid.shape, bool)
        # Dx = np.empty(self.xGrid.shape, dtype=np.float32)
        # Dx.fill(np.nan)
        # Dy = np.empty(self.xGrid.shape, dtype=np.float32)
        # Dy.fill(np.nan)
        
    
    def iter_variables(self, itr):
        Scale = self.ChipSize0X / self.ChipSizeUniX[itr]
        if Scale<1:
            dstShape = (int(self.xGrid.shape[0]*Scale),int(self.xGrid.shape[1]*Scale))
            self.xGrid0 = cv2.resize(self.xGrid.astype(np.float32),dstShape[::-1],interpolation=cv2.INTER_AREA)
            self.yGrid0 = cv2.resize(self.yGrid.astype(np.float32),dstShape[::-1],interpolation=cv2.INTER_AREA)
            self.mask0 = cv2.resize(self.mask.astype(np.uint8),dstShape[::-1],interpolation=cv2.INTER_NEAREST).astype(bool)
            self.noDataMask0 = cv2.resize(self.noDataMask.astype(np.uint8),dstShape[::-1],interpolation=cv2.INTER_NEAREST).astype(bool)

            self.po2vr0 = cv2.resize(self.po2vr.astype(np.float32),dstShape[::-1],interpolation=cv2.INTER_AREA)
            self.po2va0 = cv2.resize(self.po2va.astype(np.float32),dstShape[::-1],interpolation=cv2.INTER_AREA)

            self.SearchLimitX0 = np.ceil(cv2.resize(self.SearchLimitX, dstShape[::-1]))
            self.SearchLimitY0 = np.ceil(cv2.resize(self.SearchLimitY, dstShape[::-1]))
        else:
            self.xGrid0 = self.xGrid
            self.yGrid0 = self.yGrid
            self.mask0 = self.mask
            self.noDataMask0 = self.noDataMask
            self.po2vr0 = self.po2vr
            self.po2va0 = self.po2va
            self.SearchLimitX0 = self.SearchLimitX
            self.SearchLimitY0 = self.SearchLimitY
        
        if np.mod(self.ChipSizeUniX[itr],2) == 0:
            self.xGrid0 = np.round(self.xGrid0+0.5)-0.5
            self.yGrid0 = np.round(self.yGrid0+0.5)-0.5

        
        self.overSampleRatio = self.OverSampleRatio[self.ChipSizeUniX[itr]]
        ChipSizeXF = self.ChipSizeUniX[itr]
        ChipSizeYF = np.float32(np.round(ChipSizeXF*self.ScaleChipSizeY/2)*2)
        self.SearchLimitY0 = np.float32(np.round(self.SearchLimitY0*self.ScaleChipSizeY))
        
        self.ChipSizeX0 = np.ones(self.xGrid0.shape, dtype=np.float32) * ChipSizeXF
        self.ChipSizeY0 = np.ones(self.xGrid0.shape, dtype=np.float32) * ChipSizeYF
        # self.SearchLimitX0 = np.ones(self.xGrid0.shape, dtype=np.float32) * self.SearchLimitX0
        # self.SearchLimitY0 = np.ones(self.xGrid0.shape, dtype=np.float32) * self.SearchLimitY0
        
        # adding loop specific variables
        # self.DxF = np.empty(self.xGrid0.shape,dtype=np.float32)
        # self.DxF.fill(np.nan)
        # self.DyF = self.DxF.copy()

        SLx_max = np.max(self.SearchLimitX0)
        Px = int(np.max(ChipSizeXF)/2 + SLx_max + 2)
        SLy_max = np.max(self.SearchLimitY0)
        Py = int(np.max(ChipSizeYF)/2 + SLy_max + 2)
        self.I1F = np.lib.pad(self.I1,((Py,Py),(Px,Px)),'constant')
        self.I2F = np.lib.pad(self.I2,((Py,Py),(Px,Px)),'constant')
        self.xGrid0 += (Px + 0.5)
        self.yGrid0 += (Py + 0.5)
        
    def get_chip_ref(self, ii, jj):
        clx = np.floor(self.ChipSizeX0[ii,jj]/2)
        ChipRangeX = slice(int(-clx + self.xGrid0[ii,jj]) , int(clx + self.xGrid0[ii,jj]))
        cly = np.floor(self.ChipSizeY0[ii,jj]/2)
        ChipRangeY = slice(int(-cly + self.yGrid0[ii,jj]) , int(cly + self.yGrid0[ii,jj]))
        ChipI = self.I2F[ChipRangeY,ChipRangeX]
        
        SearchRangeX = slice(int(-clx - self.SearchLimitX0[ii,jj] + self.xGrid0[ii,jj]) , int(clx + self.SearchLimitX0[ii,jj] - 1 + self.xGrid0[ii,jj]))
        SearchRangeY = slice(int(-cly - self.SearchLimitY0[ii,jj] + self.yGrid0[ii,jj]) , int(cly + self.SearchLimitY0[ii,jj] - 1 + self.yGrid0[ii,jj]))
        RefI = self.I1F[SearchRangeY,SearchRangeX]
        
        minChipI = np.min(ChipI)
        if minChipI < 0:
            ChipI = ChipI - minChipI
        
        minRefI = np.min(RefI)
        if minRefI < 0:
            RefI = RefI - minRefI
        
        return ChipI, RefI
    
    def checkValid(self, ii, jj):
        # Not vald
        sr_check = (self.SearchLimitX0[ii,jj] <= 0) | (self.SearchLimitY0[ii,jj] <= 0)
        mask_check = (self.noDataMask0[ii,jj] == 1) | (self.mask0[ii,jj] == 0)
        return sr_check | mask_check
    
    def compute_vel(self, dx, dy, ii, jj):
        dx = (dx-self.SearchLimitX0[ii,jj])*self.po2vr0[ii,jj]
        dy = (dy-self.SearchLimitY0[ii,jj])*self.po2va0[ii,jj]
        return dx, dy