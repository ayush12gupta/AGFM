#!/usr/bin/env python3

import pdb
import subprocess
import re
import string
import sys
import cv2, time
import numpy as np
from util import upscale


class glacier_autoRIFT:
    '''
    Class for mapping regular geographic grid on radar imagery.
    '''
    
    def preprocess_filt_wal(self):
        '''
        Do the pre processing using wallis filter (10 min vs 15 min in Matlab).
        '''
        
        self.zeroMask = (self.IMG[0] == 0)    
        kernel = np.ones((self.WallisFilterWidth,self.WallisFilterWidth), dtype=np.float32)
        
        for i in range(len(self.IMG)):
            m = cv2.filter2D(self.IMG[i],-1,kernel,borderType=cv2.BORDER_CONSTANT)/np.sum(kernel)
            m2 = (self.IMG[i])**2
            m2 = cv2.filter2D(m2,-1,kernel,borderType=cv2.BORDER_CONSTANT)/np.sum(kernel)
            s = np.sqrt(m2 - m**2) * np.sqrt(np.sum(kernel)/(np.sum(kernel)-1.0))
            self.IMG[i] = (self.IMG[i] - m) / s  
        
    
#    ####   obsolete definition of "preprocess_filt_hps"
#    def preprocess_filt_hps(self):
#        '''
#        Do the pre processing using (orig - low-pass filter) = high-pass filter filter (3.9/5.3 min).
#        '''
#        import cv2
#        import numpy as np
#
#        if self.zeroMask is not None:
#            self.zeroMask = (self.I1 == 0)
#
#        kernel = np.ones((self.WallisFilterWidth,self.WallisFilterWidth), dtype=np.float32)
#
#        lp = cv2.filter2D(self.I1,-1,kernel,borderType=cv2.BORDER_CONSTANT)/np.sum(kernel)
#
#        self.I1 = (self.I1 - lp)
#
#        lp = cv2.filter2D(self.I2,-1,kernel,borderType=cv2.BORDER_CONSTANT)/np.sum(kernel)
#
#        self.I2 = (self.I2 - lp)


    def preprocess_filt_hps(self):
        '''
        Do the pre processing using (orig - low-pass filter) = high-pass filter filter (3.9/5.3 min).
        '''

        kernel = -np.ones((self.WallisFilterWidth,self.WallisFilterWidth), dtype=np.float32)
        kernel[int((self.WallisFilterWidth-1)/2),int((self.WallisFilterWidth-1)/2)] = kernel.size - 1
        kernel = kernel / kernel.size

        for i in range(len(self.IMG)):
            self.IMG[i] = cv2.filter2D(self.IMG[i],-1,kernel,borderType=cv2.BORDER_CONSTANT)
            
    
    def preprocess_db(self):
        '''
        Do the pre processing using db scale (4 min).
        '''

        self.zeroMask = (self.IMG[0] == 0)
        self.IMG = 20.0 * np.log10(self.IMG)

        
    def preprocess_filt_sob(self):
        '''
        Do the pre processing using sobel filter (4.5/5.8 min).
        '''
        
        sobelx = cv2.getDerivKernels(1,0,self.WallisFilterWidth)
        kernelx = np.outer(sobelx[0],sobelx[1])
        sobely = cv2.getDerivKernels(0,1,self.WallisFilterWidth)
        kernely = np.outer(sobely[0],sobely[1])
        
        kernel = kernelx + kernely
        
        for i in range(len(self.IMG)):
            self.IMG[i] = cv2.filter2D(self.IMG[i],-1,kernel,borderType=cv2.BORDER_CONSTANT)
    

    def preprocess_filt_lap(self):
        '''
        Do the pre processing using Laplacian filter (2.5 min / 4 min).
        '''

        self.zeroMask = (self.IMG[0] == 0)

        for i in range(len(self.IMG)):
            self.IMG[i] = 20.0 * np.log10(self.IMG[i])
            self.IMG[i] = cv2.Laplacian(self.IMG[i],-1,ksize=self.WallisFilterWidth,borderType=cv2.BORDER_CONSTANT)
                       
    
    def uniform_data_type(self):
        
        if self.DataType == 0:
            
            for i in range(len(self.IMG)):
                if self.zeroMask is not None:
                    validData = np.isfinite(self.IMG[i])
                    temp = self.IMG[i][validData]
                else:
                    temp = self.IMG[i]
                    
                S1 = np.std(temp)*np.sqrt(temp.size/(temp.size-1.0))
                M1 = np.mean(temp)
                self.IMG[i] = (self.IMG[i] - (M1 - 3*S1)) / (6*S1) * (2**8 - 0)
                self.IMG[i] = np.round(np.clip(self.IMG[i], 0, 255)).astype(np.uint8)
        
            if self.zeroMask is not None:
                # self.IMG[:,self.zeroMask] = 0
                self.IMG[:,self.zeroMask] = 0
                self.zeroMask = None
        
        elif self.DataType == 1:
            
            if self.zeroMask is not None:
                self.IMG[np.logical_not(np.isfinite(self.IMG))] = 0
                # self.I1[np.logical_not(np.isfinite(self.I1))] = 0
                # self.I2[np.logical_not(np.isfinite(self.I2))] = 0
            
            self.IMG = self.IMG.astype(np.float32)
            # self.I1 = self.I1.astype(np.float32)
            # self.I2 = self.I2.astype(np.float32)
            
            if self.zeroMask is not None:
                self.IMG[:,self.zeroMask] = 0
                # self.I1[:,self.zeroMask] = 0
                # self.I2[:,self.zeroMask] = 0
                self.zeroMask = None

        else:
            sys.exit('invalid data type for the image pair which must be unsigned integer 8 or 32-bit float')


    def glacier_autorift(self):
        '''
        Do the actual processing.
        '''
        from scipy import ndimage

        ChipSizeUniX = np.unique(np.append(np.unique(self.ChipSizeMinX), np.unique(self.ChipSizeMaxX)))
        ChipSizeUniX = np.delete(ChipSizeUniX,np.where(ChipSizeUniX == 0)[0])
        # print(ChipSizeUniX, self.ChipSize0X, 'chipsz')
        if np.any(np.mod(ChipSizeUniX,self.ChipSize0X) != 0):
            sys.exit('chip sizes must be even integers of ChipSize0')

        ChipRangeX = self.ChipSize0X * np.array([1,2,4,8,16,32,64],np.float32)
#        ChipRangeX = ChipRangeX[ChipRangeX < (2**8 - 1)]
        if np.max(ChipSizeUniX) > np.max(ChipRangeX):
            sys.exit('max each chip size is out of range')

        ChipSizeUniX = ChipRangeX[(ChipRangeX >= np.min(ChipSizeUniX)) & (ChipRangeX <= np.max(ChipSizeUniX))]

        maxScale = np.max(ChipSizeUniX) / self.ChipSize0X

        if (np.mod(self.xGrid.shape[0],maxScale) != 0)|(np.mod(self.xGrid.shape[1],maxScale) != 0):
            message = 'xgrid and ygrid have an incorect size ' + str(self.xGrid.shape) + ' for nested search, they must have dimensions that an interger multiple of ' + str(maxScale)
            sys.exit(message)
        
        self.xGrid = self.xGrid.astype(np.float32)
        self.yGrid = self.yGrid.astype(np.float32)
        
        if np.size(self.SearchLimitX) == 1:
            self.SearchLimitX = np.ones(self.xGrid.shape, np.float32) * np.round(self.SearchLimitX)
        else:
            self.SearchLimitX = self.SearchLimitX.astype(np.float32)
        if np.size(self.SearchLimitY) == 1:
            self.SearchLimitY = np.ones(self.xGrid.shape, np.float32) * np.round(self.SearchLimitY)
        else:   
            self.SearchLimitY = self.SearchLimitY.astype(np.float32)
        if np.size(self.ChipSizeMinX) == 1:
            self.ChipSizeMinX = np.ones(self.xGrid.shape, np.float32) * np.round(self.ChipSizeMinX)
        else:
            self.ChipSizeMinX = self.ChipSizeMinX.astype(np.float32)
        if np.size(self.ChipSizeMaxX) == 1:
            self.ChipSizeMaxX = np.ones(self.xGrid.shape, np.float32) * np.round(self.ChipSizeMaxX)
        else:
            self.ChipSizeMaxX = self.ChipSizeMaxX.astype(np.float32)
        if self.mask is None:
            self.mask = np.ones(self.xGrid.shape, bool)

        ChipSizeX = np.zeros(self.xGrid.shape, np.float32)
        InterpMask = np.zeros(self.xGrid.shape, bool)
        Dx = np.empty(self.xGrid.shape, dtype=np.float32)
        Dx.fill(np.nan)
        Dy = np.empty(self.xGrid.shape, dtype=np.float32)
        Dy.fill(np.nan)
        SNR = np.empty(self.xGrid.shape, dtype=np.float32)
        SNR.fill(0)

        Flag = 3  
        
        if self.ChipSize0X > self.GridSpacingX:
            if np.mod(self.ChipSize0X,self.GridSpacingX) != 0:
                sys.exit('when GridSpacing < smallest allowable chip size (ChipSize0), ChipSize0 must be integer multiples of GridSpacing')
            else:
                ChipSize0_GridSpacing_oversample_ratio = int(self.ChipSize0X / self.GridSpacingX)
        else:
            ChipSize0_GridSpacing_oversample_ratio = 1
                
        DispFiltC = DISP_FILT()
        overlap_c = np.max((1 - self.sparseSearchSampleRate / ChipSize0_GridSpacing_oversample_ratio,0))
        DispFiltC.FracValid = self.FracValid * (1 - overlap_c) + overlap_c**2
        DispFiltC.FracSearch = self.FracSearch
        DispFiltC.FiltWidth = (self.FiltWidth - 1) * ChipSize0_GridSpacing_oversample_ratio + 1
        DispFiltC.Iter = self.Iter - 1
        DispFiltC.MadScalar = self.MadScalar
        DispFiltC.colfiltChunkSize = self.colfiltChunkSize

        DispFiltF = DISP_FILT()
        overlap_f = 1 - 1 / ChipSize0_GridSpacing_oversample_ratio
        DispFiltF.FracValid = self.FracValid * (1 - overlap_f) + overlap_f**2
        DispFiltF.FracSearch = self.FracSearch
        DispFiltF.FiltWidth = (self.FiltWidth - 1) * ChipSize0_GridSpacing_oversample_ratio + 1
        DispFiltF.Iter = self.Iter
        DispFiltF.MadScalar = self.MadScalar
        DispFiltF.colfiltChunkSize = self.colfiltChunkSize

        for i in range(ChipSizeUniX.__len__()-1, -1, -1):
            print("Chipsize", ChipSizeUniX[i])
            
            # Nested grid setup: chip size being ChipSize0X no need to resize, otherwise has to resize the arrays
            if self.ChipSize0X != ChipSizeUniX[i]:
                Scale = self.ChipSize0X / ChipSizeUniX[i]
                dstShape = (int(self.xGrid.shape[0]*Scale),int(self.xGrid.shape[1]*Scale))
                xGrid0 = cv2.resize(self.xGrid.astype(np.float32),dstShape[::-1],interpolation=cv2.INTER_AREA)
                yGrid0 = cv2.resize(self.yGrid.astype(np.float32),dstShape[::-1],interpolation=cv2.INTER_AREA)
                
                if np.mod(ChipSizeUniX[i],2) == 0:
                    xGrid0 = np.round(xGrid0+0.5)-0.5
                    yGrid0 = np.round(yGrid0+0.5)-0.5
                else:
                    xGrid0 = np.round(xGrid0)
                    yGrid0 = np.round(yGrid0)

                # M0 = (ChipSizeX != 0) & (self.ChipSizeMinX <= ChipSizeUniX[i]) & (self.ChipSizeMaxX >= ChipSizeUniX[i])
                # M0 = colfilt(M0.copy(), (int(1/Scale*6), int(1/Scale*6)), 0, self.colfiltChunkSize)
                # M0 = cv2.resize(np.logical_not(M0).astype(np.uint8),dstShape[::-1],interpolation=cv2.INTER_NEAREST).astype(bool)
                mask0 = cv2.resize(self.mask.astype(np.uint8),dstShape[::-1],interpolation=cv2.INTER_NEAREST).astype(bool)

                SearchLimitX0 = colfilt(self.SearchLimitX.copy(), (int(1/Scale), int(1/Scale)), 0, self.colfiltChunkSize)
                SearchLimitY0 = colfilt(self.SearchLimitY.copy(), (int(1/Scale), int(1/Scale)), 0, self.colfiltChunkSize)

                SearchLimitX0 = np.ceil(cv2.resize(SearchLimitX0,dstShape[::-1]))
                SearchLimitY0 = np.ceil(cv2.resize(SearchLimitY0,dstShape[::-1]))
                # SearchLimitX0[M0] = 0
                # SearchLimitY0[M0] = 0
            else:
                SearchLimitX0 = self.SearchLimitX.copy()
                SearchLimitY0 = self.SearchLimitY.copy()
                xGrid0 = self.xGrid.copy()
                yGrid0 = self.yGrid.copy()
                mask0 = self.mask.copy()
#                M0 = (ChipSizeX == 0) & (self.ChipSizeMinX <= ChipSizeUniX[i]) & (self.ChipSizeMaxX >= ChipSizeUniX[i])
#                SearchLimitX0[np.logical_not(M0)] = 0
#                SearchLimitY0[np.logical_not(M0)] = 0

            if np.logical_not(np.any(SearchLimitX0 != 0)):
                continue
            
            start_time = time.time()
            idxZero = (SearchLimitX0 <= 0) | (SearchLimitY0 <= 0) | (mask0 == 0)
            SearchLimitX0[idxZero] = 0
            SearchLimitY0[idxZero] = 0
            SearchLimitX0[(np.logical_not(idxZero)) & (SearchLimitX0 < self.minSearch)] = self.minSearch
            SearchLimitY0[(np.logical_not(idxZero)) & (SearchLimitY0 < self.minSearch)] = self.minSearch

            if ((xGrid0.shape[0] - 2)/(self.sparseSearchSampleRate * ChipSize0_GridSpacing_oversample_ratio) < 5) | ((xGrid0.shape[1] - 2)/(self.sparseSearchSampleRate * ChipSize0_GridSpacing_oversample_ratio) < 5):
                Flag = 2
                return Flag
        
            # Setup for coarse search: sparse sampling / resize
            rIdxC = slice((self.sparseSearchSampleRate * ChipSize0_GridSpacing_oversample_ratio)-1,xGrid0.shape[0],(self.sparseSearchSampleRate * ChipSize0_GridSpacing_oversample_ratio))
            cIdxC = slice((self.sparseSearchSampleRate * ChipSize0_GridSpacing_oversample_ratio)-1,xGrid0.shape[1],(self.sparseSearchSampleRate * ChipSize0_GridSpacing_oversample_ratio))
            xGrid0C = xGrid0[rIdxC,cIdxC]
            yGrid0C = yGrid0[rIdxC,cIdxC]
            
            if np.remainder((self.sparseSearchSampleRate * ChipSize0_GridSpacing_oversample_ratio),2) == 0:
                filtWidth = (self.sparseSearchSampleRate * ChipSize0_GridSpacing_oversample_ratio) + 1
            else:
                filtWidth = (self.sparseSearchSampleRate * ChipSize0_GridSpacing_oversample_ratio)

            SearchLimitX0C = colfilt(SearchLimitX0.copy(), (int(filtWidth), int(filtWidth)), 0, self.colfiltChunkSize)
            SearchLimitY0C = colfilt(SearchLimitY0.copy(), (int(filtWidth), int(filtWidth)), 0, self.colfiltChunkSize)
            SearchLimitX0C = SearchLimitX0C[rIdxC,cIdxC]
            SearchLimitY0C = SearchLimitY0C[rIdxC,cIdxC]

            # Coarse search
            SubPixFlag = False
            ChipSizeXC = ChipSizeUniX[i]
            ChipSizeYC = np.float32(np.round(ChipSizeXC*self.ScaleChipSizeY/2)*2)
            # Change
            SearchLimitY0C = np.float32(np.round(SearchLimitY0C*self.ScaleChipSizeY))
            
            if type(self.OverSampleRatio) is dict:
                overSampleRatio = self.OverSampleRatio[ChipSizeUniX[i]]
            else:
                overSampleRatio = self.OverSampleRatio

            if self.IMG.dtype == np.uint8:
                DxC, DyC, SNR_C = arImgDisp_u(self.IMG, xGrid0C.copy(), yGrid0C.copy(), ChipSizeXC, ChipSizeYC, SearchLimitX0C, SearchLimitY0C, self.pad_img, SubPixFlag, overSampleRatio, self.StepSize, self.MultiThread)
            elif self.IMG.dtype == np.float32:
                DxC, DyC, SNR_C = arImgDisp_s(self.IMG, xGrid0C.copy(), yGrid0C.copy(), ChipSizeXC, ChipSizeYC, SearchLimitX0C, SearchLimitY0C, self.pad_img, SubPixFlag, overSampleRatio, self.StepSize, self.MultiThread)
            else:
                sys.exit('invalid data type for the image pair which must be unsigned integer 8 or 32-bit float')
            
            chpt1 = time.time()
            print("Time-taken for coarse search:", chpt1-start_time)
            
            # M0C is the mask for reliable estimates after coarse search, MC is the mask after disparity filtering, MC2 is the mask after area closing for fine search
            M0C = np.logical_not(np.isnan(DxC))

            MC = DispFiltC.filtDisp(DxC.copy(), DyC.copy(), SearchLimitX0C.copy(), SearchLimitY0C.copy(), M0C.copy(), overSampleRatio)

            MC[np.logical_not(M0C)] = False
    
            ROIC = (SearchLimitX0C > 0)
            CoarseCorValidFac = np.sum(MC[ROIC]) / np.sum(M0C[ROIC])
            if (CoarseCorValidFac < self.CoarseCorCutoff):
                continue
            
            MC2 = ndimage.distance_transform_edt(np.logical_not(MC)) < self.BuffDistanceC
            dstShape = (int(MC2.shape[0]*(self.sparseSearchSampleRate*ChipSize0_GridSpacing_oversample_ratio)),int(MC2.shape[1]*(self.sparseSearchSampleRate*ChipSize0_GridSpacing_oversample_ratio)))

            MC2 = cv2.resize(MC2.astype(np.uint8),dstShape[::-1],interpolation=cv2.INTER_NEAREST).astype(bool)
            if np.logical_not(np.all(MC2.shape == SearchLimitX0.shape)):
                rowAdd = SearchLimitX0.shape[0] - MC2.shape[0]
                colAdd = SearchLimitX0.shape[1] - MC2.shape[1]
                if rowAdd>0:
                    MC2 = np.append(MC2,MC2[-rowAdd:,:],axis=0)
                if colAdd>0:
                    MC2 = np.append(MC2,MC2[:,-colAdd:],axis=1)

            SearchLimitX0[np.logical_not(MC2)] = 0
            SearchLimitY0[np.logical_not(MC2)] = 0

            # Fine Search
            SubPixFlag = True
            ChipSizeXF = ChipSizeUniX[i]
            ChipSizeYF = np.float32(np.round(ChipSizeXF*self.ScaleChipSizeY/2)*2)
            SearchLimitY0 = np.float32(np.round(SearchLimitY0*self.ScaleChipSizeY))

            print("Starting Fine Search")
            if self.IMG.dtype == np.uint8:
                DxF, DyF, SNR_F = arImgDisp_u(self.IMG, xGrid0.copy(), yGrid0.copy(), ChipSizeXF, ChipSizeYF, SearchLimitX0, SearchLimitY0, self.pad_img, SubPixFlag, overSampleRatio, self.StepSize, self.MultiThread)
            elif self.IMG.dtype == np.float32:
                DxF, DyF, SNR_F = arImgDisp_s(self.IMG, xGrid0.copy(), yGrid0.copy(), ChipSizeXF, ChipSizeYF, SearchLimitX0, SearchLimitY0, self.pad_img, SubPixFlag, overSampleRatio, self.StepSize, self.MultiThread)
            else:
                sys.exit('invalid data type for the image pair which must be unsigned integer 8 or 32-bit float')
            
            M0 = DispFiltF.filtDisp(DxF.copy(), DyF.copy(), SearchLimitX0.copy(), SearchLimitY0.copy(), np.logical_not(np.isnan(DxF)), overSampleRatio)
            DxF[np.logical_not(M0)] = np.nan
            DyF[np.logical_not(M0)] = np.nan
            SNR_F[np.logical_not(M0)] = 0
            
            # Light interpolation with median filtered values: DxFM (filtered) and DxF (unfiltered)
            DxFM = colfilt(DxF.copy(), (self.fillFiltWidth, self.fillFiltWidth), 3, self.colfiltChunkSize)
            DyFM = colfilt(DyF.copy(), (self.fillFiltWidth, self.fillFiltWidth), 3, self.colfiltChunkSize)
            SNR_FM = colfilt(SNR_F.copy(), (self.fillFiltWidth, self.fillFiltWidth), 3, self.colfiltChunkSize)
            
            # M0 is mask for original valid estimates, MF is mask for filled ones, MM is mask where filtered ones exist for filling
            MF = np.zeros(M0.shape, dtype=bool)
            MM = np.logical_not(np.isnan(DxFM))

            for j in range(3):
                foo = MF | M0   # initial valid estimates
                foo1 = (cv2.filter2D(foo.astype(np.float32),-1,np.ones((3,3)),borderType=cv2.BORDER_CONSTANT) >= 6) | foo     # 1st area closing followed by the 2nd (part of the next line calling OpenCV)
                fillIdx = np.logical_not(bwareaopen(np.logical_not(foo1).astype(np.uint8), 5)) & np.logical_not(foo) & MM
                MF[fillIdx] = True
                DxF[fillIdx] = DxFM[fillIdx]
                DyF[fillIdx] = DyFM[fillIdx]
                SNR_F[fillIdx] = SNR_FM[fillIdx]
            
            # Below is for replacing the valid estimates with the bicubic filtered values for robust and accurate estimation
            if self.ChipSize0X == ChipSizeUniX[i]:
                idx = np.isnan(Dx) & np.logical_not(DxF)
                snr_gain_idx = (SNR_F >= SNR)
                Dx[snr_gain_idx|idx] = DxF[snr_gain_idx|idx]
                Dy[snr_gain_idx|idx] = DyF[snr_gain_idx|idx]
                SNR[snr_gain_idx|idx] = SNR_F[snr_gain_idx|idx]
                ChipSizeX[(M0|MF)&(snr_gain_idx|idx)] = ChipSizeUniX[i]
                InterpMask[MF&(snr_gain_idx|idx)] = True
            else:
                Scale = ChipSizeUniX[i] / self.ChipSize0X
                dstShape = (int(Dx.shape[0]/Scale),int(Dx.shape[1]/Scale))
            
                # DxF0 (filtered) / Dx (unfiltered) is the result from earlier iterations, DxFM (filtered) / DxF (unfiltered) is that of the current iteration
                # first colfilt nans within 2-by-2 area (otherwise 1 nan will contaminate all 4 points)
                DxF0 = colfilt(Dx.copy(),(int(Scale+1),int(Scale+1)),2, self.colfiltChunkSize)
                # then resize to half size using area (similar to averaging) to match the current iteration
                DxF0 = cv2.resize(DxF0,dstShape[::-1],interpolation=cv2.INTER_AREA)
                DyF0 = colfilt(Dy.copy(),(int(Scale+1),int(Scale+1)),2, self.colfiltChunkSize)
                DyF0 = cv2.resize(DyF0,dstShape[::-1],interpolation=cv2.INTER_AREA)
                SNR_F0 = colfilt(SNR.copy(),(int(Scale+1),int(Scale+1)),2, self.colfiltChunkSize)
                SNR_F0 = cv2.resize(SNR_F0,dstShape[::-1],interpolation=cv2.INTER_AREA)
                
                # Note this DxFM is almost the same as DxFM (same variable) in the light interpolation (only slightly better); however, only small portion of it will be used later at locations specified by M0 and MF that are determined in the light interpolation. So even without the following two lines, the final Dx and Dy result is still the same.
                # to fill out all of the missing values in DxF
                DxFM = colfilt(DxF.copy(), (5,5), 3, self.colfiltChunkSize)
                DyFM = colfilt(DyF.copy(), (5,5), 3, self.colfiltChunkSize)
                SNR_FM = colfilt(SNR_F.copy(), (5,5), 3, self.colfiltChunkSize)
                
                # !!! fill the current-iteration result with previously determined reliable estimates that are not searched in the current iteration
                # idx = np.isnan(DxF) & np.logical_not(np.isnan(DxF0))
                # idx = np.logical_not(np.isnan(DxF0))
                # DxFM[idx] = DxF0[idx]
                # DyFM[idx] = DyF0[idx]
                # SNR_FM[idx] = SNR_F0[idx]
                
                # !!!! Strong interpolation: use filtered estimates wherever the unfiltered estimates do not exist
                # idx = np.isnan(DxF) & np.logical_not(np.isnan(DxFM))
                # DxF[idx] = DxFM[idx]
                # DyF[idx] = DyFM[idx]
                # SNR_F[idx] = SNR_FM[idx]
                
                dstShape = (Dx.shape[0],Dx.shape[1])
                DxF = cv2.resize(DxF,dstShape[::-1],interpolation=cv2.INTER_CUBIC)
                DyF = cv2.resize(DyF,dstShape[::-1],interpolation=cv2.INTER_CUBIC)
                SNR_F = cv2.resize(SNR_FM,dstShape[::-1],interpolation=cv2.INTER_CUBIC)
                
                chpt2 = time.time()
                print(chpt2-chpt1)
                
                snr_gain_idx = (SNR_F >= SNR) 
                idx = np.isnan(Dx) & np.logical_not(DxF)
                print("SNR gain", ChipSizeUniX[i], snr_gain_idx.sum())
                
                MF = cv2.resize(MF.astype(np.uint8),dstShape[::-1],interpolation=cv2.INTER_NEAREST).astype(bool)
                M0 = cv2.resize(M0.astype(np.uint8),dstShape[::-1],interpolation=cv2.INTER_NEAREST).astype(bool)
                
                idxRaw = M0 & (ChipSizeX == 0)
                idxFill = MF & (ChipSizeX == 0)
                valid_idx = (idxRaw | idxFill) & (snr_gain_idx | idx)
                ChipSizeX[valid_idx] = ChipSizeUniX[i]
                InterpMask[idxFill & snr_gain_idx] = True
                Dx[valid_idx] = DxF[valid_idx]
                Dy[valid_idx] = DyF[valid_idx]
                SNR[valid_idx] = SNR_F[valid_idx]
                
        Flag = 1
        ChipSizeY = np.round(ChipSizeX * self.ScaleChipSizeY /2) * 2
        self.Dx = Dx
        self.Dy = Dy
        self.SNR = SNR
        self.InterpMask = InterpMask
        self.Flag = Flag
        self.ChipSizeX = ChipSizeX
        self.ChipSizeY = ChipSizeY
        np.save('chipx.npy', ChipSizeX)

    
    def runAutorift(self):
        '''
        quick processing routine which calls autorift main function (user can define their own way by mimicing the workflow here).
        '''
        # truncate the grid to fit the nested grid
        if np.size(self.ChipSizeMaxX) == 1:
            chopFactor = self.ChipSizeMaxX / self.ChipSize0X
        else:
            chopFactor = np.max(self.ChipSizeMaxX) / self.ChipSize0X
        rlim = int(np.floor(self.xGrid.shape[0] / chopFactor) * chopFactor)
        clim = int(np.floor(self.xGrid.shape[1] / chopFactor) * chopFactor)
        self.origSize = self.xGrid.shape

        self.xGrid = np.round(self.xGrid[0:rlim,0:clim]) + 0.5
        self.yGrid = np.round(self.yGrid[0:rlim,0:clim]) + 0.5
        
        # truncate the search limits as well if they exist
        if np.size(self.SearchLimitX) != 1:
            self.SearchLimitX = self.SearchLimitX[0:rlim,0:clim]
            self.SearchLimitY = self.SearchLimitY[0:rlim,0:clim]
                
        # truncate the chip sizes as well if they exist
        if np.size(self.ChipSizeMaxX) != 1:
            self.ChipSizeMaxX = self.ChipSizeMaxX[0:rlim,0:clim]
            self.ChipSizeMinX = self.ChipSizeMinX[0:rlim,0:clim]
            
        if self.mask is not None:
            self.mask = self.mask[0:rlim,0:clim]
        
        # call autoRIFT main function
        self.glacier_autorift()
        

    def __init__(self):
        
        super(glacier_autoRIFT, self).__init__()
        
        ##Input related parameters
        self.IMG = None
        # self.I2 = None
        self.xGrid = None
        self.yGrid = None
        self.mask = None
        self.pad_img = 400
        self.origSize = None
        self.zeroMask = None

        ##Output file
        self.Dx = None
        self.Dy = None
        self.SNR = None
        self.InterpMask = None
        self.Flag = None
        self.ChipSizeX = None
        self.ChipSizeY = None

        ##Parameter list
        self.StepSize = 12
        self.WallisFilterWidth = 21
        self.ChipSizeMinX = 32
        self.ChipSizeMaxX = 64
        self.ChipSize0X = 32
        self.GridSpacingX = 32
        self.ScaleChipSizeY = 1
        self.SearchLimitX = 12
        self.SearchLimitY = 12
        self.SkipSampleX = 32
        self.SkipSampleY = 32
        self.fillFiltWidth = 3
        self.minSearch = 6
        self.sparseSearchSampleRate = 4
        self.FracValid = 8/25
        self.FracSearch = 0.20
        self.FiltWidth = 5
        self.Iter = 3
        self.MadScalar = 4
        self.colfiltChunkSize = 4
        self.BuffDistanceC = 8
        self.CoarseCorCutoff = 0.01
        self.OverSampleRatio = 16
        self.DataType = 0
        self.MultiThread = 0



var_dict = {}

def initializer(IMG, xGrid, yGrid, SearchLimitX, SearchLimitY, ChipSizeX, ChipSizeY, StepSize):
    var_dict['IMG'] = IMG
    var_dict['xGrid'] = xGrid
    var_dict['yGrid'] = yGrid
    var_dict['SearchLimitX'] = SearchLimitX
    var_dict['SearchLimitY'] = SearchLimitY
    var_dict['ChipSizeX'] = ChipSizeX
    var_dict['ChipSizeY'] = ChipSizeY
    var_dict['StepSize'] = StepSize
    
    
# def compute_dx(ChipShape, ChipI, RefShape, RefI, scale=1, interp='spline'):
#     ChipI = ChipI.reshape(ChipShape)
#     RefI = RefI.reshape(RefShape)
    
#     res = []
#     # start = time.time()
#     for i in range(len(ChipI[:])):
#         res.append(cv2.matchTemplate(RefI[i], ChipI[i], cv2.TM_CCOEFF_NORMED))
#     res2 = np.array(res).mean(0)
    
#     res = res[0]
#     (_, maxVal, _, maxLoc) = cv2.minMaxLoc(res)
#     res[res<0] = 0
#     ratio = maxVal/np.abs(res).mean()

#     (_, maxVal2, _, maxLoc2) = cv2.minMaxLoc(res2)
#     res2[res2<0] = 0
#     ratio2 = maxVal2/np.abs(res2).mean()

#     if ratio2>ratio:
#         res = res2
#         ratio = ratio2
#         maxVal = maxVal2
#         maxLoc = maxLoc2
    
#     dX, dY = maxLoc
#     H, W = res.shape
#     if (maxVal<0.1)|(dX<=W*0.2)|(dY<=H*0.2)|(dX>=W*0.8)|(dY>=H*0.8):
#         return np.nan, np.nan, 0

#     res_up = res.copy()

#     if scale > 1:
#         xstart = maxLoc[0]-3
#         ystart = maxLoc[1]-3
#         res_up = upscale(res[maxLoc[1]-3:maxLoc[1]+4, maxLoc[0]-3:maxLoc[0]+4], scale, type=interp)
    
#         (_, maxVal, _, maxLoc) = cv2.minMaxLoc(res_up)
#         dX, dY = maxLoc
        
#         dX = (dX/scale) + xstart
#         dY = (dY/scale) + ystart
        
#     if np.isnan([dX, dY]).any():
#         return np.nan, np.nan, 0

#     return dX, dY, ratio


def compute_dx(ChipShape, ChipI, RefShape, RefI, scale=1, interp='spline'):
    ChipI = ChipI.reshape(ChipShape)
    RefI = RefI.reshape(RefShape)
    
    res = []
    # start = time.time()
    for i in range(len(ChipI[:])):
        res.append(cv2.matchTemplate(RefI[i], ChipI[i], cv2.TM_CCOEFF_NORMED))
    res = np.array(res).mean(0)
    
    (_, maxVal, _, maxLoc) = cv2.minMaxLoc(res)
    res[res<0] = 0
    ratio = maxVal/np.abs(res).mean()
    
    dX, dY = maxLoc
    H, W = res.shape
    if (maxVal<0.1)|(dX<=W*0.2)|(dY<=H*0.2)|(dX>=W*0.8)|(dY>=H*0.8):
        return np.nan, np.nan, 0

    res_up = res.copy()

    if scale > 1:
        xstart = maxLoc[0]-3
        ystart = maxLoc[1]-3
        res_up = upscale(res[maxLoc[1]-3:maxLoc[1]+4, maxLoc[0]-3:maxLoc[0]+4], scale, type=interp)

        (_, maxVal, _, maxLoc) = cv2.minMaxLoc(res_up)
        dX, dY = maxLoc
        
        dX = (dX/scale) + xstart
        dY = (dY/scale) + ystart
        
    if np.isnan([dX, dY]).any():
        return np.nan, np.nan, 0

    return dX, dY, ratio


def unpacking_column_u(tup):
    
    import numpy as np

    jj, SubPixFlag, oversample, in_shape, I_shape, padX, padY = tup
    
    xGrid = var_dict['xGrid'].copy()
    yGrid = var_dict['yGrid'].copy()
    
    Dx = np.empty(xGrid[:,jj].shape,dtype=np.float32)
    Dx.fill(np.nan)
    Dy = Dx.copy()
    SNR = np.zeros_like(Dx)

    for ii in range(in_shape[0]):
        
        if (var_dict['SearchLimitX'][ii,jj] == 0) & (var_dict['SearchLimitY'][ii,jj] == 0):
            continue
    
        # remember motion terms Dx and Dy correspond to I1 relative to I2 (reference)
        clx = np.floor(var_dict['ChipSizeX'][ii,jj]/2)
        ChipRangeX = slice(int(-clx + padX + xGrid[ii,jj]) , int(clx + padX + xGrid[ii,jj]))
        cly = np.floor(var_dict['ChipSizeY'][ii,jj]/2)
        ChipRangeY = slice(int(-cly + padY + yGrid[ii,jj]) , int(cly + padY + yGrid[ii,jj]))
        ChipI = var_dict['IMG'][:-var_dict['StepSize'],ChipRangeY,ChipRangeX] # Convert Datatype I2
        ChipShape = ChipI.shape
        
        SearchRangeX = slice(int(-clx - var_dict['SearchLimitX'][ii,jj] + padX + xGrid[ii,jj]) , int(clx + padX + var_dict['SearchLimitX'][ii,jj] - 1 + xGrid[ii,jj]))
        SearchRangeY = slice(int(-cly - var_dict['SearchLimitY'][ii,jj] + padY + yGrid[ii,jj]) , int(cly + padY + var_dict['SearchLimitY'][ii,jj] - 1 + yGrid[ii,jj]))
        RefI = var_dict['IMG'][var_dict['StepSize']:,SearchRangeY,SearchRangeX]  # Convert Datatype I1
        RefShape = RefI.shape
        
        minChipI = np.min(ChipI)
        if minChipI < 0:
            ChipI = ChipI - minChipI
        if np.all([ar==ar[0,0] for ar in ChipI]):
            continue
        
        minRefI = np.min(RefI)
        if minRefI < 0:
            RefI = RefI - minRefI
        if np.all([ar==ar[0,0] for ar in RefI]):
            continue

        if SubPixFlag:
            # call C++
            Dx[ii], Dy[ii], SNR[ii] = np.float32(compute_dx(ChipShape, ChipI.astype('uint8'), RefShape, RefI.astype('uint8'), scale=oversample, interp='spline'))
            # SNR[ii] = 10*np.log10(SNR[ii]**2)
        else:
            # call C++
            Dx[ii], Dy[ii], SNR[ii] = np.float32(compute_dx(ChipShape, ChipI.astype('uint8'), RefShape, RefI.astype('uint8'), scale=1, interp='spline'))
            # SNR[ii] = 10*np.log10(SNR[ii]**2)
            
    return Dx, Dy, SNR


def unpacking_column_s(tup):
    
    import numpy as np

    jj, SubPixFlag, oversample, in_shape, I_shape, padX, padY = tup

    xGrid = var_dict['xGrid'].copy()
    yGrid = var_dict['yGrid'].copy()

    Dx = np.empty(xGrid[:,jj].shape,dtype=np.float32)
    Dx.fill(np.nan)
    Dy = Dx.copy()
    SNR = np.zeros_like(Dx)
    
    for ii in range(in_shape[0]):
        
        if (var_dict['SearchLimitX'][ii,jj] == 0) & (var_dict['SearchLimitY'][ii,jj] == 0):
            continue
    
        # remember motion terms Dx and Dy correspond to I1 relative to I2 (reference)
        clx = np.floor(var_dict['ChipSizeX'][ii,jj]/2)
        ChipRangeX = slice(int(-clx + padX + xGrid[ii,jj]) , int(clx + padX + xGrid[ii,jj]))
        cly = np.floor(var_dict['ChipSizeY'][ii,jj]/2)
        ChipRangeY = slice(int(-cly + padY + yGrid[ii,jj]) , int(cly + padY + yGrid[ii,jj]))
        ChipI = var_dict['IMG'][:-var_dict['StepSize'],ChipRangeY,ChipRangeX]
        ChipShape = ChipI.shape
                
        SearchRangeX = slice(int(-clx - var_dict['SearchLimitX'][ii,jj] + padX + xGrid[ii,jj]) , int(clx + padX + var_dict['SearchLimitX'][ii,jj] - 1 + xGrid[ii,jj]))
        SearchRangeY = slice(int(-cly - var_dict['SearchLimitY'][ii,jj] + padY + yGrid[ii,jj]) , int(cly + padY + var_dict['SearchLimitY'][ii,jj] - 1 + yGrid[ii,jj]))
        RefI = var_dict['IMG'][var_dict['StepSize']:,SearchRangeY,SearchRangeX]
        RefShape = RefI.shape
        
        minChipI = np.min(ChipI)
        if minChipI < 0:
            ChipI = ChipI - minChipI
        if np.all([ar==ar[0,0] for ar in ChipI]):
            continue
        
        minRefI = np.min(RefI)
        if minRefI < 0:
            RefI = RefI - minRefI
        if np.all([ar==ar[0,0] for ar in RefI]):
            continue
        
        # print("!! ChipI", ChipI.shape, RefI.shape)
        
        if SubPixFlag:
            # call C++
            Dx[ii], Dy[ii], SNR[ii] = np.float32(compute_dx(ChipShape, ChipI.ravel(), RefShape, RefI.ravel(), scale=oversample, interp='spline'))
            # SNR[ii] = 10*np.log10(SNR[ii]**2)
        else:
            # call C++
            Dx[ii], Dy[ii], SNR[ii] = np.float32(compute_dx(ChipShape, ChipI.ravel(), RefShape, RefI.ravel(), scale=1, interp='spline'))
            # SNR[ii] = 10*np.log10(SNR[ii]**2)
            
    return Dx, Dy, SNR




def arImgDisp_u(IMG, xGrid, yGrid, ChipSizeX, ChipSizeY, SearchLimitX, SearchLimitY, pad_sz, SubPixFlag, oversample, step_size, MultiThread):

    import numpy as np
    import multiprocessing as mp

    if np.size(SearchLimitX) == 1:
        if np.logical_not(isinstance(SearchLimitX,np.float32) & isinstance(SearchLimitY,np.float32)):
            sys.exit('SearchLimit must be float')
    else:
        if np.logical_not((SearchLimitX.dtype == np.float32) & (SearchLimitY.dtype == np.float32)):
            sys.exit('SearchLimit must be float')

    if np.size(ChipSizeX) == 1:
        if np.logical_not(isinstance(ChipSizeX,np.float32) & isinstance(ChipSizeY,np.float32)):
            sys.exit('ChipSize must be float')
    else:
        if np.logical_not((ChipSizeX.dtype == np.float32) & (ChipSizeY.dtype == np.float32)):
            sys.exit('ChipSize must be float')
    
    if np.any(np.mod(ChipSizeX,2) != 0) | np.any(np.mod(ChipSizeY,2) != 0):
        sys.exit('it is better to have ChipSize = even number')
    
    if np.any(np.mod(SearchLimitX,1) != 0) | np.any(np.mod(SearchLimitY,1) != 0):
        sys.exit('SearchLimit must be an integar value')
    
    if np.any(SearchLimitX < 0) | np.any(SearchLimitY < 0):
        sys.exit('SearchLimit cannot be negative')

    if np.any(np.mod(ChipSizeX,4) != 0) | np.any(np.mod(ChipSizeY,4) != 0):
        sys.exit('ChipSize should be evenly divisible by 4')

    if np.size(SearchLimitX) == 1:
        SearchLimitX = np.ones(xGrid.shape, dtype=np.float32) * SearchLimitX
    bwareaopen
    if np.size(SearchLimitY) == 1:
        SearchLimitY = np.ones(xGrid.shape, dtype=np.float32) * SearchLimitY

    if np.size(ChipSizeX) == 1:
        ChipSizeX = np.ones(xGrid.shape, dtype=np.float32) * ChipSizeX
    
    if np.size(ChipSizeY) == 1:
        ChipSizeY = np.ones(xGrid.shape, dtype=np.float32) * ChipSizeY

    # convert from cartesian X-Y to matrix X-Y: X no change, Y from up being positive to down being positive
    SLx_max = np.max(SearchLimitX)
    Px = int(np.max(ChipSizeX)/2 + SLx_max + 2)
    SLy_max = np.max(SearchLimitY)
    Py = int(np.max(ChipSizeY)/2 + SLy_max + 2)
    padx = (pad_sz - Px)
    pady = (pad_sz - Py)

    # adjust center location by the padarray size and 0.5 is added because we need to extract the chip centered at X+1 with -chipsize/2:chipsize/2-1, which equivalently centers at X+0.5 (X is the original grid point location). So for even chipsize, always returns offset estimates at (X+0.5).
    xGrid += (Px + 0.5)
    yGrid += (Py + 0.5)

    Dx = np.empty(xGrid.shape,dtype=np.float32)
    Dx.fill(np.nan)
    Dy = Dx.copy()
    SNR = np.zeros_like(Dx)

    if MultiThread == 0:
        for jj in range(xGrid.shape[1]):
            if np.all(SearchLimitX[:,jj] == 0) & np.all(SearchLimitY[:,jj] == 0):
                continue
            Dx1 = Dx[:,jj]
            Dy1 = Dy[:,jj]
            SNR1 = SNR[:,jj]
            for ii in range(xGrid.shape[0]):
                if (SearchLimitX[ii,jj] == 0) & (SearchLimitY[ii,jj] == 0):
                    continue
                
                # remember motion terms Dx and Dy correspond to I1 relative to I2 (reference)
                clx = np.floor(ChipSizeX[ii,jj]/2)
                ChipRangeX = slice(int(-clx + padx + xGrid[ii,jj]) , int(clx + padx + xGrid[ii,jj]))
                cly = np.floor(ChipSizeY[ii,jj]/2)
                ChipRangeY = slice(int(-cly + pady + yGrid[ii,jj]) , int(cly + pady + yGrid[ii,jj]))
                ChipI = IMG[:-step_size,ChipRangeY,ChipRangeX]
                ChipShape = ChipI.shape

                SearchRangeX = slice(int(-clx - SearchLimitX[ii,jj] + padx + xGrid[ii,jj]) , int(clx + padx + SearchLimitX[ii,jj] - 1 + xGrid[ii,jj]))
                SearchRangeY = slice(int(-cly - SearchLimitY[ii,jj] + pady + yGrid[ii,jj]) , int(cly + pady + SearchLimitY[ii,jj] - 1 + yGrid[ii,jj]))
                RefI = IMG[step_size:,SearchRangeY,SearchRangeX]
                RefShape = RefI.shape
                
                minChipI = np.min(ChipI)
                if minChipI < 0:
                    ChipI = ChipI - minChipI
                if np.all([ar==ar[0,0] for ar in ChipI]):
                    continue
                
                minRefI = np.min(RefI)
                if minRefI < 0:
                    RefI = RefI - minRefI
                if np.all([ar==ar[0,0] for ar in RefI]):
                    continue

                if SubPixFlag:
                    # call C++
                    Dx1[ii], Dy1[ii], SNR1[ii] = np.float32(compute_dx(ChipShape, ChipI.ravel(), RefShape, RefI.ravel(), scale=oversample, interp='spline'))
                    # SNR1[ii] = 10*np.log10(SNR1[ii]**2)
#                   # call Python
                else:
                    # call C++
                    Dx1[ii], Dy1[ii], SNR1[ii] = np.float32(compute_dx(ChipShape, ChipI.ravel(), RefShape, RefI.ravel(), scale=1, interp='spline'))
                    # SNR1[ii] = 10*np.log10(SNR1[ii]**2)
#                   # call Python

    else:
        #   Preparation for parallel
        in_shape = xGrid.shape
        I_shape = IMG[0].shape

        num_cores = min(mp.cpu_count(), 64)

        chunk_inputs = [(jj, SubPixFlag, oversample, in_shape, I_shape, padx, pady)
                        for jj in range(in_shape[1])]

        with mp.Pool(initializer=initializer, initargs=(IMG, xGrid, yGrid, SearchLimitX, SearchLimitY, ChipSizeX, ChipSizeY, step_size), processes=num_cores) as pool:
            Dx, Dy, SNR = zip(*pool.map(unpacking_column_u, chunk_inputs))

        Dx = np.array(Dx).T
        Dy = np.array(Dy).T
        SNR = np.array(SNR).T

    # add back 1) I1 (RefI) relative to I2 (ChipI) initial offset Dx0 and Dy0, and
    #          2) RefI relative to ChipI has a left/top boundary offset of -SearchLimitX and -SearchLimitY
    idx = np.logical_not(np.isnan(Dx))
    Dx[idx] += (-SearchLimitX[idx])
    Dy[idx] += (-SearchLimitY[idx])
    
    # convert from matrix X-Y to cartesian X-Y: X no change, Y from down being positive to up being positive
    Dy = -Dy
    
    return Dx, Dy, SNR






def arImgDisp_s(IMG, xGrid, yGrid, ChipSizeX, ChipSizeY, SearchLimitX, SearchLimitY, pad_sz, SubPixFlag, oversample, step_size, MultiThread):
    
    import numpy as np
    import multiprocessing as mp

    if np.size(SearchLimitX) == 1:
        if np.logical_not(isinstance(SearchLimitX,np.float32) & isinstance(SearchLimitY,np.float32)):
            sys.exit('SearchLimit must be float')
    else:
        if np.logical_not((SearchLimitX.dtype == np.float32) & (SearchLimitY.dtype == np.float32)):
            sys.exit('SearchLimit must be float')

    if np.size(ChipSizeX) == 1:
        if np.logical_not(isinstance(ChipSizeX,np.float32) & isinstance(ChipSizeY,np.float32)):
            sys.exit('ChipSize must be float')
    else:
        if np.logical_not((ChipSizeX.dtype == np.float32) & (ChipSizeY.dtype == np.float32)):
            sys.exit('ChipSize must be float')

    if np.any(np.mod(ChipSizeX,2) != 0) | np.any(np.mod(ChipSizeY,2) != 0):
        sys.exit('it is better to have ChipSize = even number')
    
    if np.any(np.mod(SearchLimitX,1) != 0) | np.any(np.mod(SearchLimitY,1) != 0):
        sys.exit('SearchLimit must be an integar value')

    if np.any(SearchLimitX < 0) | np.any(SearchLimitY < 0):
        sys.exit('SearchLimit cannot be negative')
    
    if np.any(np.mod(ChipSizeX,4) != 0) | np.any(np.mod(ChipSizeY,4) != 0):
        sys.exit('ChipSize should be evenly divisible by 4')
        
    if np.size(SearchLimitX) == 1:
        SearchLimitX = np.ones(xGrid.shape, dtype=np.float32) * SearchLimitX
    
    if np.size(SearchLimitY) == 1:
        SearchLimitY = np.ones(xGrid.shape, dtype=np.float32) * SearchLimitY

    if np.size(ChipSizeX) == 1:
        ChipSizeX = np.ones(xGrid.shape, dtype=np.float32) * ChipSizeX
    
    if np.size(ChipSizeY) == 1:
        ChipSizeY = np.ones(xGrid.shape, dtype=np.float32) * ChipSizeY
    
    SLx_max = np.max(SearchLimitX)
    Px = int(np.max(ChipSizeX)/2 + SLx_max + 2)
    SLy_max = np.max(SearchLimitY)
    Py = int(np.max(ChipSizeY)/2 + SLy_max + 2)
    padx = (pad_sz - Px)
    pady = (pad_sz - Py)
    
    # adjust center location by the padarray size and 0.5 is added because we need to extract the chip centered at X+1 with -chipsize/2:chipsize/2-1, which equivalently centers at X+0.5 (X is the original grid point location). So for even chipsize, always returns offset estimates at (X+0.5).
    xGrid += (Px + 0.5)
    yGrid += (Py + 0.5)
    
    Dx = np.empty(xGrid.shape,dtype=np.float32)
    Dx.fill(np.nan)
    Dy = Dx.copy()
    SNR = np.zeros_like(Dx)
    
    if MultiThread == 0:
        for jj in range(xGrid.shape[1]):
            if np.all(SearchLimitX[:,jj] == 0) & np.all(SearchLimitY[:,jj] == 0):
                continue
            
            Dx1 = Dx[:,jj]
            Dy1 = Dy[:,jj]
            SNR1 = SNR[:,jj]
            for ii in range(xGrid.shape[0]):
                if (SearchLimitX[ii,jj] == 0) & (SearchLimitY[ii,jj] == 0):
                    continue
                
                # remember motion terms Dx and Dy correspond to I1 relative to I2 (reference)
                clx = np.floor(ChipSizeX[ii,jj]/2)
                ChipRangeX = slice(int(-clx + padx + xGrid[ii,jj]) , int(clx + padx + xGrid[ii,jj]))
                cly = np.floor(ChipSizeY[ii,jj]/2)
                ChipRangeY = slice(int(-cly + pady + yGrid[ii,jj]) , int(cly + pady + yGrid[ii,jj]))
                ChipI = IMG[:-step_size,ChipRangeY,ChipRangeX]
                ChipShape = ChipI.shape
                
                SearchRangeX = slice(int(-clx - SearchLimitX[ii,jj] + padx + xGrid[ii,jj]) , int(clx + padx + SearchLimitX[ii,jj] - 1 + xGrid[ii,jj]))
                SearchRangeY = slice(int(-cly - SearchLimitY[ii,jj] + pady + yGrid[ii,jj]) , int(cly + pady + SearchLimitY[ii,jj] - 1 + yGrid[ii,jj]))
                RefI = IMG[step_size:,SearchRangeY,SearchRangeX]
                RefShape = RefI.shape
                
                minChipI = np.min(ChipI)
                if minChipI < 0:
                    ChipI = ChipI - minChipI
                if np.all([ar==ar[0,0] for ar in ChipI]):
                    continue
                
                minRefI = np.min(RefI)
                if minRefI < 0:
                    RefI = RefI - minRefI
                if np.all([ar==ar[0,0] for ar in RefI]):
                    continue
            
                if SubPixFlag:
                    # call C++
                    Dx1[ii], Dy1[ii], SNR1[ii] = np.float32(compute_dx(ChipShape, ChipI.ravel(), RefShape, RefI.ravel(), scale=oversample, interp='spline'))
#                   # call Python
                else:
                    # call C++
                    Dx1[ii], Dy1[ii], SNR1[ii] = np.float32(compute_dx(ChipShape, ChipI.ravel(), RefShape, RefI.ravel(), scale=1, interp='spline'))
#                   # call Python
                    # ssif not np.isnan(Dx1[ii]):
                # SNR1[ii] = 10*np.log10(SNR1[ii]**2)
                
    else:
        #   Preparation for parallel
        in_shape = xGrid.shape
        I_shape = IMG[0].shape
        
        num_cores = min(mp.cpu_count(), 64)
    
        chunk_inputs = [(jj, SubPixFlag, oversample, in_shape, I_shape, padx, pady)
                        for jj in range(in_shape[1])]
            
        with mp.Pool(initializer=initializer, initargs=(IMG, xGrid, yGrid, SearchLimitX, SearchLimitY, ChipSizeX, ChipSizeY, step_size), processes=num_cores) as pool:
            Dx, Dy, SNR = zip(*pool.map(unpacking_column_s, chunk_inputs))
                
        Dx = np.array(Dx).T
        Dy = np.array(Dy).T
        SNR = np.array(SNR).T
        # print('loop_sub', SNR[np.logical_not(np.isnan(Dx))].max(), SNR[np.logical_not(np.isnan(Dx))].mean(), SNR[np.logical_not(np.isnan(Dx))].min())
        # print(Dx.shape, Dy.shape, SNR.shape, xGrid.shape)

    # add back 1) I1 (RefI) relative to I2 (ChipI) initial offset Dx0 and Dy0, and
    #          2) RefI relative to ChipI has a left/top boundary offset of -SearchLimitX and -SearchLimitY
    idx = np.logical_not(np.isnan(Dx))
    Dx[idx] += (-SearchLimitX[idx])
    Dy[idx] += (-SearchLimitY[idx])
    
    # convert from matrix X-Y to cartesian X-Y: X no change, Y from down being positive to up being positive
    Dy = -Dy
    
    return Dx, Dy, SNR





################## Chunked version of column filter
def colfilt(A, kernelSize, option, chunkSize=4):
    
    from skimage.util import view_as_windows as viewW
    import numpy as np
    
    chunkInds = int(A.shape[1]/chunkSize)
    chunkRem = A.shape[1] - chunkSize * chunkInds
    
    O = 0
    
    for ii in range(chunkSize):
        startInds = ii*chunkInds
        if ii == chunkSize-1:
            endInds = (ii+1)*chunkInds + chunkRem
        else:
            endInds = (ii+1)*chunkInds
        
        if (ii == 0)&(ii == chunkSize-1):
            A1 = np.lib.pad(A[:,startInds:endInds],((int((kernelSize[0]-1)/2),int((kernelSize[0]-1)/2)),(int((kernelSize[1]-1)/2),int((kernelSize[1]-1)/2))),mode='constant',constant_values=np.nan)
        else:
            if ii == 0:
                A1 = np.lib.pad(A[:,startInds:np.min((endInds+int((kernelSize[1]-1)/2),A.shape[1]-1))],((int((kernelSize[0]-1)/2),int((kernelSize[0]-1)/2)),(int((kernelSize[1]-1)/2),np.max((0,endInds+int((kernelSize[1]-1)/2)-A.shape[1]+1)))),mode='constant',constant_values=np.nan)
            elif ii == chunkSize-1:
                A1 = np.lib.pad(A[:,np.max((0,startInds-int((kernelSize[1]-1)/2))):endInds],((int((kernelSize[0]-1)/2),int((kernelSize[0]-1)/2)),(np.max((0,0-startInds+int((kernelSize[1]-1)/2))),int((kernelSize[1]-1)/2))),mode='constant',constant_values=np.nan)
            else:
                A1 = np.lib.pad(A[:,np.max((0,startInds-int((kernelSize[1]-1)/2))):np.min((endInds+int((kernelSize[1]-1)/2),A.shape[1]-1))],((int((kernelSize[0]-1)/2),int((kernelSize[0]-1)/2)),(np.max((0,0-startInds+int((kernelSize[1]-1)/2))),np.max((0,endInds+int((kernelSize[1]-1)/2)-A.shape[1]+1)))),mode='constant',constant_values=np.nan)

        B = viewW(A1, kernelSize).reshape(-1,kernelSize[0]*kernelSize[1]).T[:,::1]
    
        Adtype = A1.dtype
        Ashape = A1.shape
        del A1

        output_size = (Ashape[0]-kernelSize[0]+1,Ashape[1]-kernelSize[1]+1)
        C = np.zeros((B.shape[1],),dtype=Adtype)
    
        if option == 0:#    max
            C = np.nanmax(B,axis=0)
            del B
            C = C.reshape(output_size)
        elif option == 1:#  min
            C = np.nanmin(B,axis=0)
            del B
            C = C.reshape(output_size)
        elif option == 2:#  mean
            C = np.nanmean(B,axis=0)
            del B
            C = C.reshape(output_size)
        elif option == 3:#  median
            C = np.nanmedian(B,axis=0, overwrite_input=True)
            del B
            C = C.reshape(output_size)
        elif option == 4:#  range
            C = np.nanmax(B,axis=0) - np.nanmin(B,axis=0)
            del B
            C = C.reshape(output_size)
        elif option == 6:#  MAD (Median Absolute Deviation)
            m = B.shape[0]
            D = np.zeros((B.shape[1],),dtype=Adtype)
            D = np.nanmedian(B,axis=0)
            D = np.abs(B - np.dot(np.ones((m,1),dtype=Adtype), np.array([D])))
            del B
            C = np.nanmedian(D,axis=0, overwrite_input=True)
            del D
            C = C.reshape(output_size)
        elif option[0] == 5:#  displacement distance count with option[1] being the threshold
            m = B.shape[0]
            c = int(np.round((m + 1) / 2)-1)
            #        c = 0
            D = np.abs(B - np.dot(np.ones((m,1),dtype=Adtype), np.array([B[c,:]])))
            del B
            C = np.sum(D<option[1],axis=0)
            del D
            C = C.reshape(output_size)
        else:
            sys.exit('invalid option for columnwise neighborhood filtering')

        C = C.astype(Adtype)

        if np.isscalar(O):
            O = C.copy()
        else:
            O = np.append(O,C,axis=1)

    return O



class DISP_FILT:
    
    def __init__(self):
        ##filter parameters; try different parameters to decide how much fine-resolution estimates we keep, which can make the final images smoother
        
        self.FracValid = 8/25
        self.FracSearch = 0.20
        self.FiltWidth = 5
        self.Iter = 3
        self.MadScalar = 4
        self.colfiltChunkSize = 4
    
    
    def filtDisp(self, Dx, Dy, SearchLimitX, SearchLimitY, M, OverSampleRatio):
        
        import numpy as np
        
        if np.mod(self.FiltWidth,2) == 0:
            sys.exit('NDC filter width must be an odd number')
        
        dToleranceX = self.FracValid * self.FiltWidth**2
        dToleranceY = self.FracValid * self.FiltWidth**2
#        pdb.set_trace()
        Dx = Dx / SearchLimitX
        Dy = Dy / SearchLimitY
        
        DxMadmin = np.ones(Dx.shape) / OverSampleRatio / SearchLimitX * 2;
        DyMadmin = np.ones(Dy.shape) / OverSampleRatio / SearchLimitY * 2;
        
        for i in range(self.Iter):
            Dx[np.logical_not(M)] = np.nan
            Dy[np.logical_not(M)] = np.nan
            M = (colfilt(Dx.copy(), (self.FiltWidth, self.FiltWidth), (5,self.FracSearch), self.colfiltChunkSize) >= dToleranceX) & (colfilt(Dy.copy(), (self.FiltWidth, self.FiltWidth), (5,self.FracSearch), self.colfiltChunkSize) >= dToleranceY)

#        if self.Iter == 3:
#            pdb.set_trace()

        for i in range(np.max([self.Iter-1,1])):
            Dx[np.logical_not(M)] = np.nan
            Dy[np.logical_not(M)] = np.nan
            
            DxMad = colfilt(Dx.copy(), (self.FiltWidth, self.FiltWidth), 6, self.colfiltChunkSize)
            DyMad = colfilt(Dy.copy(), (self.FiltWidth, self.FiltWidth), 6, self.colfiltChunkSize)
            
            DxM = colfilt(Dx.copy(), (self.FiltWidth, self.FiltWidth), 3, self.colfiltChunkSize)
            DyM = colfilt(Dy.copy(), (self.FiltWidth, self.FiltWidth), 3, self.colfiltChunkSize)

            M = (np.abs(Dx - DxM) <= np.maximum(self.MadScalar * DxMad, DxMadmin)) & (np.abs(Dy - DyM) <= np.maximum(self.MadScalar * DyMad, DyMadmin)) & M
        
        return M



def bwareaopen(image,size1):
    
    import numpy as np
    from skimage import measure
    
    # now identify the objects and remove those above a threshold
    labels, N = measure.label(image,connectivity=2,return_num=True)
    label_size = [(labels == label).sum() for label in range(N + 1)]
    
    # now remove the labels
    for label,size in enumerate(label_size):
        if size < size1:
            image[labels == label] = 0

    return image


