import cv2
import numpy as np
import scipy.optimize as opt
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
    # start = time.time()
    res2 = BSpline(res, scale)
    (_, maxVal, _, maxLoc) = cv2.minMaxLoc(res2)
    # chpt1 = time.time()
    # print(chpt1 - start, "1")
    if not fit:
        # print("fit not")
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

    # chpt2 = time.time()
    # print(chpt2 - chpt1, "2")
    # find the optimal Gaussian parameters
    popt, pcov = opt.curve_fit(twoD_Gaussian, (x, y), data, p0=initial_guess)
    # chpt3 = time.time()
    # print(chpt3 - chpt2, "3")
    
    data_fitted = twoD_Gaussian((x, y), *popt).reshape(h,w)
    # chpt4 = time.time()
    # print(chpt4 - chpt3, "4")
    
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

