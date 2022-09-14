import os, shutil
import glob
import time
import pandas as pd
import numpy as np
from osgeo import gdal, osr


def read_vals(fn, nodat=None, band=1):
    ds = gdal.Open(fn)
    disp = ds.GetRasterBand(band).ReadAsArray()
    nodata = np.isnan(disp)
    disp[nodata] = -32767
    if nodat is not None:
        disp[nodat] = -32767
    ds = None
    nodata = (nodat|nodata)
    return disp, nodata


def velocity_correction(csv_fn, tif_dir, deltaT):

    data = pd.read_csv(csv_fn, header=0)
    num = len(data)
    N = ((num + 6)/3) - 1
    A = np.zeros((num, N))
    nodata = None
    L = []
    for i, row in data.iterrows():
        id = row['Id'].split('_')[1]
        st, ed = int(id[0]), int(id[1])
        A[i, st:ed] = 1
        val, nodata = read_vals(os.path.join(tif_dir, id, 'velocity.tif'), nodat=nodata)
        L.append(val)
    
    shp = L[0].shape
    L = np.array(L).reshape(num, -1)
    U, s, VT = np.linalg.svd(A, full_matrices=True)
    S_inv = np.zeros(A.shape).T
    S_inv[:N-1,:N-1] = np.diag(1/s[:N-1])
    # Ab = VT.T.dot(S_inv).dot(U.T)
    x_svd = ((VT.T.dot(S_inv).dot(U.T))@L).reshape(N, shp[0], shp[1])
    x_svd[:,nodata] = -32767
