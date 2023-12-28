from osgeo import gdal, osr       # GDAL support for reading virtual files
import os                         # To create and remove directories
import matplotlib.pyplot as plt   # For plotting
import matplotlib.dates as mdates
import numpy as np
import rasterio
from glob import glob
from rasterio.mask import mask
from tqdm.auto import tqdm
import pandas as pd
import scipy.linalg as la
from datetime import datetime
from scipy.linalg import lstsq as sp_lstsq


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
                return float(line.split(":")[1][1:])


def get_incang(path):
    with open(path, 'r') as fp:
        for l_no, line in enumerate(fp):
            # search string
            if 'Incidence Angle' in line:
                return float(line.split(":")[1][1:])


def read_timeseries(csv_path, base_path, skip_idx=[], xoff=0, yoff=0, xsize=100, ysize=100):
    df = pd.read_csv(csv_path)
    dates = sorted(list(set(df['Master'].tolist() + df['Slave'].tolist())))
    dates = [str(dt) for dt in dates]
    dates = [f'{dt[:4]}/{dt[4:6]}/{dt[6:]}' for dt in dates]
    raster_list = sorted(glob(f'{base_path}/stack*'))
    # print(raster_list[0])
    v_r, v_a = read_velocity(f'{raster_list[0]}/velocity_res.tif', xoff, yoff, xsize, ysize)
    mask = np.zeros_like(v_r)
    # mask = None
    V_a, V_r = {}, {}
    Delta_T = {}
    config = {}
    cnt = 0
    for i, file in enumerate(raster_list[:]):
        idx1, idx2 = file.split('/')[-1].split('_')[1:]
        idx1, idx2 = int(idx1), int(idx2)
        if len(skip_idx)!=0:
            if (idx1 in skip_idx) | (idx2 in skip_idx):
                continue
        
        t1 = datetime.strptime(dates[idx1], '%Y/%m/%d')
        t2 = datetime.strptime(dates[idx2], '%Y/%m/%d')
        delta_t = (t2 - t1).days
        
        v_r, v_a = read_velocity(f'{file}/velocity_res.tif', xoff, yoff, xsize, ysize)
        if delta_t <= 12:
            mask += np.isnan(v_r)
            # continue
        
        if idx1 not in V_a:
            V_a[idx1] = {}
        if idx2 not in V_a[idx1]:
            V_a[idx1][idx2] = v_a
        else:
            print("Repetition on", idx1, idx2)

        if idx1 not in V_r:
            V_r[idx1] = {}
        if idx2 not in V_r[idx1]:
            V_r[idx1][idx2] = v_r
        else:
            print("Repetition on", idx1, idx2)

        # print(delta_t, idx1, idx2)
        if idx1 not in Delta_T:
            Delta_T[idx1] = {}
        if idx2 not in Delta_T[idx1]:
            Delta_T[idx1][idx2] = delta_t
        else:
            print("Delta_T: Repetition on", idx1, idx2)

        cnt += 1

    config['dates'] = dates
    config['V_range'] = V_r
    config['V_azimuth'] = V_a
    config['deltaT'] = Delta_T
    config['Mask'] = np.logical_not(mask!=0)
    config['Num_pairs'] = cnt
    config['flight_angle'] = (get_flightang(f'{base_path}/testGeogrid.txt'))#*(np.pi/180)
    config['incid_angle'] = get_incang(f'{base_path}/testGeogrid.txt')#*(np.pi/180)
    return config


def interpolate_single(config, diff=3):
    indexs = config['V_range'].keys()
    for i in range(max(indexs)+(2-diff)):
        V_r = config['V_range'][i][i+diff]
        V_a = config['V_azimuth'][i][i+diff]
        mask = np.isnan(V_r)
        # print(mask.shape, V_r.shape, i, i+diff)

        V_r_avg = np.zeros_like(V_r)
        V_a_avg = np.zeros_like(V_a)
        for j in range(diff):
            # print(i+j, i+j+1)
            V_r_avg += (config['V_range'][i+j][i+j+1]*config['deltaT'][i+j][i+j+1])
            V_a_avg += (config['V_azimuth'][i+j][i+j+1]*config['deltaT'][i+j][i+j+1])

        V_r_avg /= config['deltaT'][i][i+diff]
        V_a_avg /= config['deltaT'][i][i+diff]
        V_r[mask] = V_r_avg[mask]
        V_a[mask] = V_a_avg[mask]
        config['V_range'][i][i+diff] = V_r
        config['V_azimuth'][i][i+diff] = V_a

    return config


def interpolate_missing(config, min_diff=2):
    indexs = config['V_range']
    max_diff = len(config['V_range'][0].keys())
    for diff in range(min_diff, max_diff+1):
        config = interpolate_single(config, diff)

    return config


def create_Amat(config_asc, config_des):
    indexs = config_asc['V_range'].keys()
    mask = (config_des['Mask']&config_asc['Mask'])
    n_pairs = config_asc['Num_pairs']
    S = get_matS(config_asc, config_des)
    A = np.zeros((4*n_pairs, 3*len(indexs)), dtype=float)
    L = []
    cnt = 0
    pairs = []
    delt = []
    for i in range(max(indexs)+1):
        delt.append(config_asc['deltaT'][i][i+1])
        delt.append(config_asc['deltaT'][i][i+1])
        delt.append(config_asc['deltaT'][i][i+1])
        for j in config_asc['V_range'][i].keys():
            A[(cnt*4):(cnt*4) + 4, i*3:j*3] = np.hstack([S]*(j-i))
            # print((cnt*4),(cnt*4) + 4, i*3,(j)*3, i, j)
            factor = config_asc['deltaT'][i][j] # Converting to displacement
            L.append(config_asc['V_range'][i][j]*factor)
            L.append(config_asc['V_azimuth'][i][j]*factor)
            L.append(config_des['V_range'][i][j]*factor)
            L.append(config_des['V_azimuth'][i][j]*factor)
            pairs.append(f'{i}_{j}')
            cnt += 1
    
    L = np.array(L)
    delt = np.array(delt)
    return A, L, mask, pairs, delt


def create_Amat_trikhonov(config_asc, config_des, lambda_reg = 0.1):
    indexs = config_asc['V_range'].keys()
    mask = (config_des['Mask']&config_asc['Mask'])
    n_pairs = config_asc['Num_pairs']
    S = get_matS(config_asc, config_des)
    A = np.zeros(((4*n_pairs)+3*(len(indexs)-1), 3*len(indexs)), dtype=float)
    L = []
    cnt = 0
    pairs = []
    delt = []
    for i in range(max(indexs)+1):
        delt.append(config_asc['deltaT'][i][i+1])
        delt.append(config_asc['deltaT'][i][i+1])
        delt.append(config_asc['deltaT'][i][i+1])
        for j in config_asc['V_range'][i].keys():
            A[(cnt*4):(cnt+1)*4, i*3:j*3] = np.hstack([S]*(j-i))
            # print((cnt*4),(cnt*4) + 4, i*3,(j)*3, i, j)
            factor = config_asc['deltaT'][i][j] # Converting to displacement
            L.append(config_asc['V_range'][i][j]*factor)
            L.append(config_asc['V_azimuth'][i][j]*factor)
            L.append(config_des['V_range'][i][j]*factor)
            L.append(config_des['V_azimuth'][i][j]*factor)
            pairs.append(f'{i}_{j}')
            cnt += 1

    I3 = np.eye(3)
    v_zero = np.zeros_like(L[0])
    for i in range(len(indexs)-1):
        factor1 = config_asc['deltaT'][i][i+1]
        A[(4*n_pairs)+(3*i):(4*n_pairs)+(3*(i+1)), i*3:(i+1)*3] = (I3/factor1)*lambda_reg
        
        factor2 = config_asc['deltaT'][i+1][i+2]
        A[(4*n_pairs)+(3*i):(4*n_pairs)+(3*(i+1)), (i+1)*3:(i+2)*3] = (-I3/factor2)*lambda_reg
        for n in range(3):
            L.append(v_zero)
        
    L = np.array(L)
    delt = np.array(delt)
    return A, L, mask, pairs, delt


def get_Va_mat(phi):
    return np.array([np.cos(phi), np.sin(phi), 0])


def get_Vr_mat(phi, incid):
    return np.array([(-np.sin(phi)*np.sin(incid)), (np.sin(incid)*np.cos(phi)), np.cos(incid)])


def get_matS(config_asc, config_des):
    incid_asc = config_asc['incid_angle']*(np.pi/180)
    phi_asc = config_asc['flight_angle']*(np.pi/180)
    incid_des = config_des['incid_angle']*(np.pi/180)
    phi_des = config_des['flight_angle']*(np.pi/180)
    
    A = np.zeros((4, 3))
    A[0] = get_Vr_mat(phi_asc, incid_asc)
    A[1] = get_Va_mat(phi_asc)
    A[2] = get_Vr_mat(phi_des, incid_des)
    A[3] = get_Va_mat(phi_des)
    return A


def svd_solve(a, b):
    [U, s, Vt] = la.svd(a, full_matrices=False)
    r = max(np.where(s >= 1e-12)[0])
    temp = np.dot(U[:, :r].T, b) / s[:r]
    return np.dot(Vt[:r, :].T, temp)


def compute_3D_disp(A, L, mask, delT, n):
    h, w = L.shape[1:]
    margin = 0
    V_nev = np.zeros((A.shape[1], h, w))
    # print(V_nev.shape)
    
    for i in tqdm(range(margin, h-margin), desc=f'loop for i = {n}'):
        for j in range(margin, w-margin):
            if mask[i,j]==False:
                continue
            # v_nev = pos_solve(A, L[:,i,j], 1e-25)
            # v_nev = svd_solve(A, L[:,i,j])
            v_nev, _, _, _ = sp_lstsq(A, L[:,i,j], lapack_driver='gelsy', check_finite=False)
            V_nev[:,i,j] = v_nev/delT
            # xhat_lstsq = naive_solve(A, L[:,i,j], 1e-25)
            # print(np.dot(A, v_nev) - L[:,i,j])
            # print(la.norm(np.dot(A, v_nev) - L[:,i,j]))
    return V_nev


if __name__ == '__main__':
    base_path_asc = '/DATA/test_stack_asc3/offset_tracking'
    csv_path_asc = '/DATA/test_stack_asc3/image_pairs.csv'

    base_path_des = '/DATA/test_stack_des3/offset_tracking'
    csv_path_des = '/DATA/test_stack_des3/image_pairs.csv'
    
    asc_files = glob(f'{base_path_asc}/stack*/velocity_res.tif')
    ds = gdal.Open(asc_files[0])
    projs = ds.GetProjection()
    geo = ds.GetGeoTransform()
    v = ds.GetRasterBand(1).ReadAsArray()
    H, W = v.shape
    ds = None
    config_asc = read_timeseries(csv_path_asc, base_path_asc, skip_idx=[], xoff=0, yoff=0, xsize=W, ysize=10)
    config_asc = interpolate_missing(config_asc)
    config_des = read_timeseries(csv_path_des, base_path_des, skip_idx=[], xoff=0, yoff=0, xsize=W, ysize=10)
    config_des = interpolate_missing(config_des)
    A, L, mask, pairs, delt = create_Amat(config_asc, config_des)
    
    H, W = v.shape
    h, w = L.shape[1:]
    V_nev = np.zeros((A.shape[1]//3, 3, H, W))
    bSizeY = 100
    bSizeX = W
    for i in tqdm(range(0, H, bSizeY), desc='Main Loop'):
        # print(i)
        if i + bSizeY < H:
            sizeY = bSizeY
        else:
            sizeY = H - i

        for j in range(0, W, bSizeX):
            if j + bSizeX < W:
                sizeX = bSizeX
            else:
                sizeX = W - j

            config_asc = read_timeseries(csv_path_asc, base_path_asc, skip_idx=[], xoff=j, yoff=i, xsize=sizeX, ysize=sizeY)
            config_asc = interpolate_missing(config_asc)
            config_des = read_timeseries(csv_path_des, base_path_des, skip_idx=[], xoff=j, yoff=i, xsize=sizeX, ysize=sizeY)
            config_des = interpolate_missing(config_des)
            # A, L, mask, pairs, delt = create_Amat(config_asc, config_des)
            A, L, mask, pairs, delt = create_Amat_trikhonov(config_asc, config_des)
            V_nev[:,:,i:(i+sizeY),j:(j+sizeX)] = compute_3D_disp(A, L, mask, delt, i).reshape(-1, 3, sizeY, sizeX)
            
    np.save('/DATA/V_nev_corrected3.npy', V_nev)

