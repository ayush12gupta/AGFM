import numpy as np


def normalise2(data, mini, maxi):
    data -= (0.5 * (mini + maxi))
    data /= (0.25*(maxi - mini))
    return data


def aMat2d(line, pixel, numPow2D = 5):
    Ai = []
    for i in range(numPow2D+1):
        for j in range(i+1):
            Ai.append( pow(line, (i-j)) * pow(pixel, j) )
    return np.array(Ai)


def matTxmat(mat1, mat2):
    return mat1.T@mat2


def polyVal(C, line, pixel, height):
    line2 = line*line
    pixel2 = pixel*pixel
    height2 = height*height
    return C[0,0] + C[1,0]*line + C[2,0]*pixel + C[3,0]*height + C[4,0]*line*pixel + \
            C[5,0]*line*height + C[6,0]*pixel*height + C[7,0]*line2 + C[8,0]*pixel2 + C[9,0]*height2


def polyVal2d(lin, pix, beta, numcoef2D=5):
    idx = 0
    value = 0
    for i in range(numcoef2D+1):
        for j in range(i+1):
            value += (beta[idx]*pow(lin, i-j)*pow(pix, j))
            idx += 1
    
    return value


def computeSlope(xc, yc, dem, slopeCalRadius):
    slope = 0
    h = 0
    numPoints = 0
    hc = dem[yc, xc]
    halfSlopeCalRadius = slopeCalRadius // 2
    for y in range(yc - slopeCalRadius, yc + slopeCalRadius + 1, slopeCalRadius):
        for x in range(xc - slopeCalRadius, xc + slopeCalRadius + 1, halfSlopeCalRadius):
            h = dem[y, x]
            if h!=0: 
                slope += np.abs(h - hc)
                numPoints+=1

    return slope / numPoints


def getBperp(line, pixel, linMax, pixMax, rhsBperp, height=0, linMin=0, pixMin=0, hMin=0, hMax=5000):
    norm_line = normalise2(line, linMin, linMax)
    norm_pix = normalise2(pixel, pixMin, pixMax)
    norm_height = normalise2(height, hMin, hMax)
    return polyVal(rhsBperp, norm_line, norm_pix, norm_height)


def getBpar(line, pixel, linMax, pixMax, rhsBpar, height=0, linMin=0, pixMin=0, hMin=0, hMax=5000):
    norm_line = normalise2(line, linMin, linMax)
    norm_pix = normalise2(pixel, pixMin, pixMax)
    norm_height = normalise2(height, hMin, hMax)
    return polyVal(rhsBpar, norm_line, norm_pix, norm_height)


def getTheta(line, pixel, linMax, pixMax, rhsTheta, height=0, linMin=0, pixMin=0, hMin=0, hMax=5000):
    norm_line = normalise2(line, linMin, linMax)
    norm_pix = normalise2(pixel, pixMin, pixMax)
    norm_height = normalise2(height, hMin, hMax)
    return polyVal(rhsTheta, norm_line, norm_pix, norm_height)


def getThetaInc(line, pixel, linMax, pixMax, rhsThetaInc, height=0, linMin=0, pixMin=0, hMin=0, hMax=5000):
    norm_line = normalise2(line, linMin, linMax)
    norm_pix = normalise2(pixel, pixMin, pixMax)
    norm_height = normalise2(height, hMin, hMax)
    return polyVal(rhsThetaInc, norm_line, norm_pix, norm_height)


def getB(line, pixel, linMax, pixMax, rhsBperp, rhsBpar, height=0, linMin=0, pixMin=0, hMin=0, hMax=5000):
    Bperp = getBperp(line, pixel, linMax, pixMax, rhsBperp=rhsBperp, height=height, linMin=linMin, pixMin=pixMin, hMin=hMin, hMax=hMax)
    Bpar = getBperp(line, pixel, linMax, pixMax, rhsBpar=rhsBpar, height=height, linMin=linMin, pixMin=pixMin, hMin=hMin, hMax=hMax)
    return np.sqrt(Bperp**2 + Bpar**2)


def getAlpha(line, pixel, linMax, pixMax, rhsBperp, rhsBpar, rhsTheta, height=0, linMin=0, pixMin=0, hMin=0, hMax=5000):
    Bperp = getBperp(line, pixel, linMax, pixMax, rhsBperp=rhsBperp, height=height, linMin=linMin, pixMin=pixMin, hMin=hMin, hMax=hMax)
    Bpar = getBperp(line, pixel, linMax, pixMax, rhsBpar=rhsBpar, height=height, linMin=linMin, pixMin=pixMin, hMin=hMin, hMax=hMax)
    theta = getTheta(line, pixel, linMax, pixMax, rhsTheta=rhsTheta, height=height, linMin=linMin, pixMin=pixMin, hMin=hMin, hMax=hMax)
    if (Bpar==0)&(Bperp==0):
        alpha = np.nan
    else:
        alpha = theta - np.arctan2(Bpar, Bperp)
    return alpha


def getBvert(line, pixel, linMax, pixMax, rhsBperp, rhsBpar, rhsTheta, height=0, linMin=0, pixMin=0, hMin=0, hMax=5000):
    B = getB(line, pixel, linMax, pixMax, rhsBperp=rhsBperp, rhsBpar=rhsBpar, height=height, linMin=linMin, pixMin=pixMin, hMin=hMin, hMax=hMax)
    alpha = getAlpha(line, pixel, linMax, pixMax, rhsBperp=rhsBperp, rhsBpar=rhsBpar, rhsTheta=rhsTheta, height=height, linMin=linMin, pixMin=pixMin, hMin=hMin, hMax=hMax)
    return B*np.sin(alpha)


def getBhor(line, pixel, linMax, pixMax, rhsBperp, rhsBpar, rhsTheta, height=0, linMin=0, pixMin=0, hMin=0, hMax=5000):
    B = getB(line, pixel, linMax, pixMax, rhsBperp=rhsBperp, rhsBpar=rhsBpar, height=height, linMin=linMin, pixMin=pixMin, hMin=hMin, hMax=hMax)
    alpha = getAlpha(line, pixel, linMax, pixMax, rhsBperp=rhsBperp, rhsBpar=rhsBpar, rhsTheta=rhsTheta, height=height, linMin=linMin, pixMin=pixMin, hMin=hMin, hMax=hMax)
    return B*np.cos(alpha)


def getLine(y, azlook):
    line = azlook*y + (azlook//2)
    return line


def getRange(pixel, pixelsize, rng_strt):
    return rng_strt + pixel*pixelsize


def getTline(line, deltaT, nearAzimuth):
    return nearAzimuth + (line-1) * deltaT


def getPixel(x, rnglook):
    return rnglook*x + (rnglook//2)
