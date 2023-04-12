#
# Author: Ayush Gupta
# Copyright 2023
#


import numpy as np 
import os
import time
from glob import glob
import isceobj
import datetime
import logging
from osgeo import gdal
from isceobj.Planet.Planet import Planet
from isceobj.DemProc.utils import *

logger = logging.getLogger('isce.topsinsar.computedem')

numCoeffs1D = 3
nPoints = 200
numPow2D = 5
H_MAX = 5000
N_HEIGHTS = 3
N_POINTS_RNG = 10
N_POINTS_AZI = 10


# Computing reference points
def distribute_points(const, num_points):
    points = []
    wp = np.sqrt(num_points/(const['num_lines']/const['num_pix']))
    wl = (num_points/wp) if (num_points/wp)>wp else wp
    wl = int(wl) # Converting no. of window in line direction to int
    deltal = (const['num_lines']-1)/(wl-1)
    totp = int(const['num_pix']*wl)
    deltap = (totp - 1)/(num_points-1)

    p = -deltap
    l = 0
    lcnt = 0
    for i in range(num_points):
        p += deltap
        
        # Move to next line when p>=num_pix increment lcnt
        while(int(p)>=const['num_pix']):
            p-=const['num_pix']
            lcnt+=1
        
        l = lcnt*deltal
        points.append([int(l), int(p), 0])
        
    points = np.array(points)
    for i in range(nPoints):
        points[i, 0] += const['linMin']
        points[i, 1] += const['pixMin']

    return points


def generateRefPoints(self):

    const = {}

    reference_dir = self._insar.referenceSlcProduct
    secondary_dir = self._insar.secondarySlcProduct

    num_swath = len(glob(reference_dir + '/IW*.xml'))
    reference = self._insar.loadProduct(reference_dir + '/IW{0}.xml'.format(1))
    secondary = self._insar.loadProduct(secondary_dir + '/IW{0}.xml'.format(1))
    num_bursts = len(reference.bursts)

    # Storing imp variables
    masterWavelength = reference.bursts[0].radarWavelength
    slaveWavelength = secondary.bursts[0].radarWavelength
    nearRange = reference.bursts[0].startingRange
    nearAzimuth = reference.bursts[0].sensingStart
    pixelsize = reference.bursts[0].rangePixelSize
    deltaT = datetime.timedelta(seconds=reference.bursts[0].azimuthTimeInterval)

    refElp = Planet(pname='Earth').ellipsoid
    num_lines = 0
    num_pix = 0
    strt_rng, far_rng, strt_az, far_az = 0, 0, 0, 0

    # Iterating the swaths to compute number of lines and azimuth
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


    linMin, linMax = 0, num_lines
    pixMin, pixMax = 0, num_pix
    hMin, hMax = 0, H_MAX
    deltaHeight = (hMax - hMin) / (N_HEIGHTS - 1)
    
    # Storing required constants
    const['num_lines'] = num_lines
    const['num_pix'] = num_pix
    const['linMin'] = linMin
    const['linMax'] = linMax
    const['pixMin'] = pixMin
    const['pixMax'] = pixMax
    const['deltaHeight'] = deltaHeight

    points = distribute_points(const, nPoints)

    # Computing phase for these reference points
    m_picdivlambda = (-4*np.pi)/masterWavelength
    s_picdivlambda = (-4*np.pi)/slaveWavelength
    Phase = np.zeros((nPoints, N_HEIGHTS))

    for hi in range(N_HEIGHTS):
        height = hi*deltaHeight
        for i in range(nPoints):
            line, pixel = points[i, 0], points[i, 1] 

            tline = nearAzimuth + (line-1) * deltaT
            referenceSV = reference.orbit.interpolate(tline, method='hermite')

            # Getting target point position
            rng = nearRange + pixelsize*(pixel-1) # Master range
            target = reference.orbit.rdr2geo(tline, rng, height=height)
            
            # Getting slave position
            slvTime, slvrng = secondary.orbit.geo2rdr(target)
            
            Phase[i, hi] = (m_picdivlambda*rng) - (s_picdivlambda*slvrng)
            print(pixel, line, rng-slvrng, Phase[i,hi])

            if i%20==0:
                print("Height {},  Phase_min: {},  Phase_max: {},  Phase_mean: {}".format(str(height), str(Phase[i].min()), str(Phase[i].max()), str(Phase[i].mean())))

    # Making phase at h=0 as zero 
    for i in range(nPoints):
        offset = Phase[i, 0]
        Phase[i] -= offset

    return Phase, points, const


def compute1DCoeff(Phase, points, const):
    HEI = np.zeros((N_HEIGHTS, 1))
    for i in range(N_HEIGHTS):
        HEI[i,0] = i*const['deltaHeight']

        Phase_norm = Phase.copy()
        Phase_norm = normalise2(Phase_norm, Phase.min(), Phase.max())
        aMat = np.zeros((N_HEIGHTS, numCoeffs1D))

        ALPHA = np.zeros((len(points), numCoeffs1D))

        for i in range(nPoints):
            for hi in range(N_HEIGHTS):
                for j in range(numCoeffs1D):
                    aMat[hi, j] = pow(Phase_norm[i, hi], j)
                    
            nMat = matTxmat(aMat, aMat)
            rhsMat = matTxmat(aMat, HEI)
            Qx_hat = np.linalg.cholesky(nMat)
            rhsMat = np.linalg.solve(Qx_hat.T, np.linalg.solve(Qx_hat, rhsMat))
            
            # Checking for error associated with estimated coeff.
            y_hatBperp  = aMat@rhsMat
            eHatBperp  = HEI - y_hatBperp
            maxerr = (abs(eHatBperp)).max()
            if maxerr>0.01:
                print('High errors in phase', maxerr)
            
            # Storing alpha
            ALPHA[i] = rhsMat[:,0]

    return ALPHA


def compute2DCoeff(points, const, Alpha, numPow2D):
    num_Var = int(((numPow2D + 1)**2 + numPow2D + 1)*0.5)

    if num_Var>nPoints:
        print("!! ERROR: Increase num of reference points or decrease numPow2D")

    aMat_2D = np.zeros((nPoints, num_Var))

    for i in range(nPoints):
        l = normalise2(points[i, 0], const['linMin'], const['linMax'])
        p = normalise2(points[i, 1], const['pixMin'], const['pixMax'])

        aMat_2D[i] = aMat2d(l, p, numPow2D=numPow2D )

    nMat_2D = matTxmat(aMat_2D, aMat_2D)
    rhsMat_2D = matTxmat(aMat_2D, Alpha)
    Qx_hat = np.linalg.cholesky(nMat_2D)

    BETA = np.linalg.solve(Qx_hat.T, np.linalg.solve(Qx_hat, rhsMat_2D))

    return BETA


def runComputeDEM_sch(self):

    # Computing polynomial fit for phase info
    Phase, points, const = generateRefPoints(self)
    ALPHA = compute1DCoeff(Phase, points, const)
    BETA = compute2DCoeff(points, const, ALPHA, numPow2D)

    # Computing 
    inFilename = os.path.join(self._insar.mergedDirname, self._insar.mergedIfgname)
    intImage = isceobj.createIntImage()
    intImage.load(inFilename + '.xml')
    intImage.setAccessMode('read')
    intImage.createImage()
    width = intImage.getWidth()
    length = intImage.getLength()

    unw = (np.fromfile(os.path.join(self._insar.mergedDirname, self._insar.unwrappedIntFilename), dtype=np.float32).reshape(length*2, width))[1:length*2:2, :]
    out_line, out_pix = unw.shape

    firstline = const['linMin'] + (self.numberAzimuthLooks-1)/2
    firstpix = const['pixMin'] + (self.numberRangeLooks-1)/2

    pix = np.array([firstpix + i*self.numberRangeLooks for i in range(out_pix)])
    normalise2(pix, const['pixMin'], const['pixMax'])

    lin = np.array([firstline + i*self.numberAzimuthLooks for i in range(out_line)])
    normalise2(lin, const['linMin'], const['linMax'])

    mat_pix, mat_lin = np.meshgrid(pix, lin)
    pntALPHA = np.zeros((numCoeffs1D, out_line, out_pix))

    num_Var = int(((numPow2D + 1)**2 + numPow2D + 1)*0.5)
    for i in range(numCoeffs1D):
        beta = np.zeros((num_Var, 1))
        beta[:,0] = BETA[:,i]
        pntALPHA[i] = polyVal2d(mat_lin, mat_pix, beta, numPow2D)
        
    height = np.zeros_like(unw)
    unw_norm = unw.copy()
    unw_norm = normalise2(unw_norm, unw.min(), unw.max())
    for i in range(numCoeffs1D):
        height += np.multiply(pntALPHA[i],pow(unw_norm,i))

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


