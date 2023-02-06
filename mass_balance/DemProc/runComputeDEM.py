#
# Author: Piyush Agram
# Copyright 2016
#


import numpy as np 
import os
import time
import glob
import isceobj
import datetime
import logging
from osgeo import gdal
from isceobj.Planet.Planet import Planet
from isceobj.DemProc.utils import *

logger = logging.getLogger('isce.topsinsar.computedem')

numCoeffs = 10
N_HEIGHTS = 4
N_POINTS_RNG = 10
N_POINTS_AZI = 10

def model(self):
    strt = time.time()
    
    reference_dir = self._insar.referenceSlcProduct
    secondary_dir = self._insar.secondarySlcProduct

    num_swath = len(glob.glob(reference_dir + '/IW*.xml'))
    reference = self._insar.loadProduct(reference_dir + '/IW{0}.xml'.format(1))
    secondary = self._insar.loadProduct(secondary_dir + '/IW{0}.xml'.format(1))
    num_bursts = len(reference.bursts)
    masterWavelength = reference.bursts[0].radarWavelength
    nearRange = reference.bursts[0].startingRange
    nearAzimuth = reference.bursts[0].sensingStart
    pixelsize = reference.bursts[0].rangePixelSize
    deltaT = datetime.timedelta(seconds=reference.bursts[0].azimuthTimeInterval)
    
    refElp = Planet(pname='Earth').ellipsoid
    num_lines = 0
    num_pix = 0
    strt_rng, far_rng, strt_az, far_az = 0, 0, 0, 0
    
    for i in range(num_swath):
        ref = self._insar.loadProduct(reference_dir + '/IW{0}.xml'.format(i+1))
        nbursts = len(ref.bursts)
        if strt_az == 0:
            strt_az = ref.bursts[0].sensingStart
            far_az = ref.bursts[-1].sensingStop
        for j in range(nbursts):
            if strt_rng==0:
                strt_rng = ref.bursts[j].startingRange
                
            strt_az = min(strt_az, ref.bursts[j].sensingStart)
            strt_rng = min(strt_rng, ref.bursts[j].startingRange)
            far_az = max(far_az, ref.bursts[j].sensingStop)
            far_rng = max(far_rng, ref.bursts[j].farRange)
            
    num_lines = np.round((far_az-strt_az)/deltaT) + 1
    num_pix = np.round((far_rng-strt_rng)/pixelsize) + 1

    linMin = 0
    linMax = num_lines
    pixMin = 0
    pixMax = num_pix
    hMin = 0
    hMax = 5000

    cnt = 0
    deltaPixels = reference.bursts[0].numberOfSamples / N_POINTS_RNG
    deltaLines = reference.bursts[0].numberOfLines / N_POINTS_AZI
    deltaHeight = (hMax - hMin) / N_HEIGHTS

    bPerpMat = np.zeros((N_POINTS_AZI * N_POINTS_RNG * N_HEIGHTS, 1))
    bParMat = np.zeros((N_POINTS_AZI * N_POINTS_RNG * N_HEIGHTS, 1))
    thetaMat = np.zeros((N_POINTS_AZI * N_POINTS_RNG * N_HEIGHTS, 1))
    thetaIncMat = np.zeros((N_POINTS_AZI * N_POINTS_RNG * N_HEIGHTS, 1))
    aMat = np.zeros((N_POINTS_AZI * N_POINTS_RNG * N_HEIGHTS, numCoeffs))

    for k in range(N_HEIGHTS):
        height = hMin + k*deltaHeight

        for i in range(N_POINTS_AZI):
            line = linMin + i * deltaLines
            tline = nearAzimuth + (line-1) * deltaT
            referenceSV = reference.orbit.interpolate(tline, method='hermite')
            mxyz = np.array(referenceSV.getPosition())
            mvel = np.array(referenceSV.getVelocity())
            mvelunit = mvel / np.linalg.norm(mvel)

            for j in range(N_POINTS_RNG):
                pixel = pixMin + j * deltaPixels
                rng = nearRange + pixelsize*pixel
                target = reference.orbit.rdr2geo(tline, rng, height=height)
                targxyz = np.array(refElp.LLH(target[0], target[1], target[2]).ecef().tolist())

                slvTime, slvrng = secondary.orbit.geo2rdr(target)
                secondarySV = secondary.orbit.interpolateOrbit(slvTime, method='hermite')
                sxyz = np.array(secondarySV.getPosition())
                svel = np.array(secondarySV.getVelocity())
                svelunit = svel / np.linalg.norm(svel)
                angleOrb = np.arccos(np.clip(np.dot(mvelunit, svelunit), -1, 1))
                orbitConvergence = angleOrb

                bPar = rng - slvrng
                b = np.linalg.norm(mxyz-sxyz)
                r1 = (mxyz-targxyz)
                r2 = (sxyz-targxyz)
                theta = np.arccos(np.dot(mxyz, r1)/(np.linalg.norm(mxyz)*np.linalg.norm(r1))) # Look angle
                theta2 = np.arccos(np.dot(mxyz, r2)/(np.linalg.norm(mxyz)*np.linalg.norm(r2)))
                bPerp = b**2 - bPar**2

                if bPerp<0:
                    bPerp = 0
                elif theta>theta2:
                    bPerp = np.sqrt(bPerp)
                else:
                    bPerp = -np.sqrt(bPerp)

                thetaInc = np.arccos(np.dot(targxyz, r1)/(np.linalg.norm(targxyz)*np.linalg.norm(r1)))

                bPerpMat[cnt,0]    = bPerp
                bParMat[cnt,0]     = bPar
                thetaMat[cnt,0]    = theta
                thetaIncMat[cnt,0] = thetaInc

                aMat[cnt,0] = 1
                aMat[cnt,1] = normalise2(line, linMin, linMax)
                aMat[cnt,2] = normalise2(pixel, pixMin, pixMax)
                aMat[cnt,3] = normalise2(height, hMin, hMax)
                aMat[cnt,4] = normalise2(line, linMin, linMax) * normalise2(pixel, pixMin, pixMax)
                aMat[cnt,5] = normalise2(line, linMin, linMax) * normalise2(height, hMin, hMax)
                aMat[cnt,6] = normalise2(pixel, pixMin, pixMax) * normalise2(height, hMin, hMax)
                aMat[cnt,7] = normalise2(line, linMin, linMax)**2
                aMat[cnt,8] = normalise2(pixel, pixMin, pixMax)**2
                aMat[cnt,9] = normalise2(height, hMin, hMax)**2  
                cnt+=1

                if (bPar==0)&(bPerp==0):
                    alpha = np.nan
                else:
                    alpha = theta - np.arctan2(bPar, bPerp)
                
                # Horizontal / Vertical representation
                # print(360 - (90 - ((theta - alpha)*180/np.pi)))
                bH = b * np.cos(alpha)
                bV = b * np.sin(alpha)

                if bPerp==0:
                    hAmbiguity = np.inf
                else:
                    hAmbiguity = -(masterWavelength*np.linalg.norm(mxyz-targxyz)*np.sin(theta))/(2*bPerp)

    nMat = matTxmat(aMat, aMat)
    rhsBperp = matTxmat(aMat, bPerpMat)
    rhsBpar = matTxmat(aMat, bParMat)
    rhsTheta = matTxmat(aMat, thetaMat)
    rhsThetaInc = matTxmat(aMat, thetaIncMat)

    Qx_hat = np.linalg.cholesky(nMat)
    rhsBperp = np.linalg.solve(Qx_hat.T, np.linalg.solve(Qx_hat, rhsBperp))
    rhsBpar = np.linalg.solve(Qx_hat.T, np.linalg.solve(Qx_hat, rhsBpar))
    rhsTheta = np.linalg.solve(Qx_hat.T, np.linalg.solve(Qx_hat, rhsTheta))
    rhsThetaInc = np.linalg.solve(Qx_hat.T, np.linalg.solve(Qx_hat, rhsThetaInc))
    
    print("Time taken = ", time.time()-strt, "sec")
    
    y_hatBperp  = aMat@rhsBperp
    eHatBperp  = bPerpMat - y_hatBperp
    maxerr = (abs(eHatBperp)).max()
    print("Max. error bperp modeling at 3D datapoints: ", maxerr)
    
    return rhsBperp, rhsBpar, rhsTheta, rhsThetaInc, linMax, pixMax


def computeReference(self, rhsBperp, rhsBpar, rhsTheta, rhsThetaInc, linMax, pixMax):

    seedGridSize = 300
    slopeCalRadius = 4

    inFilename = os.path.join(self._insar.mergedDirname, self._insar.mergedIfgname)
    intImage = isceobj.createIntImage()
    intImage.load(inFilename + '.xml')
    intImage.setAccessMode('read')
    intImage.createImage()
    width = intImage.getWidth()
    length = intImage.getLength()

    # unw = np.memmap(os.path.join(self._insar.mergedDirname, self._insar.unwrappedIntFilename), dtype=np.float32, mode='r+', shape=(length,width))
    unw = (np.fromfile(os.path.join(self._insar.mergedDirname, self._insar.unwrappedIntFilename), dtype=np.float32).reshape(length*2, width))[1:length*2:2, :]

    ds = gdal.Open(os.path.join(self._insar.mergedDirname, 'hgt.rdr.vrt'), gdal.GA_ReadOnly)
    dem = ds.GetRasterBand(1).ReadAsArray()
    ds = None

    img = unw.copy()
    seedGridResY = int((img.shape[0] - 1 - 2 * slopeCalRadius) / (seedGridSize - 1))
    seedGridResX = int((img.shape[1] - 1 - 2 * slopeCalRadius) / (seedGridSize - 1))
    seedList = []
    demNoDataValue = 0
    reference = self._insar.loadProduct(self._insar.referenceSlcProduct + '/IW{0}.xml'.format(1))
    pixelsize = reference.bursts[0].rangePixelSize
    strt_rng = reference.bursts[0].startingRange
    masterWavelength = reference.bursts[0].radarWavelength

    for r in range(seedGridSize):
        y = r * seedGridResY + slopeCalRadius
        for c in range(seedGridSize):
            x = c * seedGridResX + slopeCalRadius;
    #         print(x, y)
            h = dem[y, x]
            if (h!=demNoDataValue)&(h > 0.0) & (unw[y,x]!=0):
                seedList.append([x, y, h, computeSlope(x, y, dem, slopeCalRadius)])

    # sort the seed list in ascending order according to the seed's slope
    seedList = sorted(seedList, key=lambda x: x[3], reverse=False);

    # get the final seed list
    maskSize = 20
    totalFinalSeeds = 300
    mask = np.zeros([maskSize, maskSize])
    finalSeedList = []
    numSeeds = 0
    for seed in seedList:
        maskX = int( seed[0] / img.shape[1] * maskSize)
        maskY = int( seed[1] / img.shape[0] * maskSize)
        if mask[maskY,maskX]==0:
            finalSeedList.append(seed)
            numSeeds+=1
            if numSeeds >= totalFinalSeeds:
                break
            mask[maskY,maskX] = 1

    #####################################

    xc = int(img.shape[1] / 2)
    deg2rad = np.pi/180
    waveNumber = 2*np.pi/masterWavelength
    # print(numSeeds)
    a, b, c, d, e, f = 0, 0, 0, 0, 0, 0
    for i in range(0, numSeeds):
        seed = finalSeedList[i]
        line = getLine(seed[1], self.numberAzimuthLooks)
        pixel_mid = getPixel(xc, self.numberRangeLooks)
        pixel = getPixel(seed[0], self.numberRangeLooks)
        phase = unw[seed[1], seed[0]]

        slantRange = strt_rng + pixel*pixelsize
        incidenceAngle = getThetaInc(line=line, pixel=pixel, linMax=linMax, pixMax=pixMax, rhsThetaInc=rhsThetaInc)
        Bn = getBperp(line=line, pixel=pixel, linMax=linMax, pixMax=pixMax, rhsBperp=rhsBperp) #bn[seed[1], seed[0]]
        Bp = getBpar(line=line, pixel=pixel, linMax=linMax, pixMax=pixMax, rhsBpar=rhsBpar) #bp[seed[1], seed[0]]
        la1 = getTheta(line=line, pixel=pixel, linMax=linMax, pixMax=pixMax, rhsTheta=rhsTheta)
        la2 = getTheta(line=line, pixel=pixel_mid, linMax=linMax, pixMax=pixMax, rhsTheta=rhsTheta)
        flatAngle = (la1 - la2)
        # print(incidenceAngle*180/np.pi, Bn, Bp, flatAngle*180/np.pi)
        # (lookAngles[seed[1],seed[0]] - lookAngles[seed[1],xc])*deg2rad
        alpha = -slantRange * np.sin(incidenceAngle) / (2 * waveNumber * (Bp * np.sin(flatAngle) + Bn * np.cos(flatAngle)))
        # print(alpha)
        # alpha = -slantRange * np.sin(incidenceAngle)/(2*waveNumber * Bn);
        a += -alpha * alpha
        b += alpha
        e += alpha * (seed[2] - alpha * phase)
        f += seed[2] - alpha * phase

    c = -b
    d = numSeeds

    refHeight = (a * f - c * e) / (a * d - c * b)
    refPhase = (e * d - b * f) / (a * d - c * b)

    return refHeight, refPhase


def runComputeDEM(self):

    inFilename = os.path.join(self._insar.mergedDirname, self._insar.mergedIfgname)
    reference = self._insar.loadProduct(self._insar.referenceSlcProduct + '/IW{0}.xml'.format(1))
    strt_rng = reference.bursts[0].startingRange
    pixelsize = reference.bursts[0].rangePixelSize
    intImage = isceobj.createIntImage()
    intImage.load(inFilename + '.xml')
    intImage.setAccessMode('read')
    intImage.createImage()
    width = intImage.getWidth()
    length = intImage.getLength()
    unw = (np.fromfile(os.path.join(self._insar.mergedDirname, self._insar.unwrappedIntFilename), dtype=np.float32).reshape(length*2, width))[1:length*2:2, :]

    rhsBperp, rhsBpar, rhsTheta, rhsThetaInc, linMax, pixMax = model(self)
    refHeight, refPhase = computeReference(self, rhsBperp, rhsBpar, rhsTheta, rhsThetaInc, linMax, pixMax)
    print(refHeight, refPhase)

    xc = int(unw.shape[1] / 2)
    waveNumber = 2*np.pi/0.056
    alpha = unw.copy()

    nrow, ncol = unw.shape
    for i in range(nrow):
        line = getLine(i, self.numberAzimuthLooks)
        
        for j in range(ncol):
            pixel_mid = getPixel(xc, self.numberRangeLooks)
            pixel = getPixel(j, self.numberRangeLooks)
            phase = unw[i, j]

            slantRange = strt_rng + pixel*pixelsize
            incidenceAngle = getThetaInc(line=line, pixel=pixel, linMax=linMax, pixMax=pixMax, rhsThetaInc=rhsThetaInc)
            Bn = getBperp(line=line, pixel=pixel, linMax=linMax, pixMax=pixMax, rhsBperp=rhsBperp) #bn[seed[1], seed[0]]
            Bp = getBpar(line=line, pixel=pixel, linMax=linMax, pixMax=pixMax, rhsBpar=rhsBpar) #bp[seed[1], seed[0]]
            la1 = getTheta(line=line, pixel=pixel, linMax=linMax, pixMax=pixMax, rhsTheta=rhsTheta)
            la2 = getTheta(line=line, pixel=pixel_mid, linMax=linMax, pixMax=pixMax, rhsTheta=rhsTheta)
            flatAngle = (la1 - la2)
            alpha[i,j] = -slantRange * np.sin(incidenceAngle) / (2 * waveNumber * (Bp * np.sin(flatAngle) + Bn * np.cos(flatAngle)))

    height = refHeight + alpha * (unw - refPhase)
    
    elevation_fn = os.path.join(self._insar.mergedDirname, 'elevation.rdr')
    elevation = np.memmap(elevation_fn, dtype=np.float32, mode='w+', shape=(length,width))
    elevation[:,:] = height

    objInt = isceobj.createImage()
    objInt.filename = elevation_fn
    objInt.setWidth(width)
    objInt.setLength(length)
    objInt.scheme = 'BIL'
    objInt.bands = 1
    objInt.dataType = 'FLOAT'
    objInt.setAccessMode('READ')
    objInt.renderHdr()
