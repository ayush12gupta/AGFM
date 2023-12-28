from shutil import copyfile, move # Utilities for copying and moving files
from osgeo import gdal          # GDAL support for reading virtual files
import os                         # To create and remove directories
import numpy as np
from glob import glob
from rasterio.mask import mask
import pandas as pd
import scipy.linalg as la
from datetime import datetime
from tqdm import tqdm
from scipy.linalg import lstsq as sp_lstsq
# from geogrid_autorift.util import 
from util import get_flightang, get_incang, read_velocity


def read_timeseries(csv_path, base_path, skip_idx=[], xoff=0, yoff=0, xsize=100, ysize=100):
    df = pd.read_csv(csv_path)
    dates = sorted(list(set(df['Master'].tolist() + df['Slave'].tolist())))
    dates = [str(dt) for dt in dates]
    dates = [f'{dt[:4]}/{dt[4:6]}/{dt[6:]}' for dt in dates]
    raster_list = sorted(glob(f'{base_path}/stack*'))
    # print(raster_list[0])
    v_r, v_a = read_velocity(f'{raster_list[0]}/velocity_res_crop.tif', xoff, yoff, xsize, ysize)
    # v_r[v_r==-32767] = np.nan
    # v_a[v_a==-32767] = np.nan
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
        
        v_r, v_a = read_velocity(f'{file}/velocity_res_crop.tif', xoff, yoff, xsize, ysize)
        v_r[v_r==-32767] = np.nan
        v_a[v_a==-32767] = np.nan
        if (idx2-idx1)>4:
            continue

        if ((idx2-idx1)>1)&(delta_t > 36):
            continue
        
        if (delta_t <= 12)|(idx2-idx1)<=1:
            mask += np.isnan(v_r) #|(v_r==-32767)
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
        if i+diff not in config['V_range'][i]:
            continue
            
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


def create_Amat_sbas(config_asc, config_des):
    indexs = config_asc['V_range'].keys()
    mask = (config_des['Mask']&config_asc['Mask'])
    n_pairs = config_asc['Num_pairs']
    # S = get_matS(config_asc, config_des)
    S = np.eye(4)
    A = np.zeros((4*n_pairs, 4*len(indexs)), dtype=float)
    L = []
    cnt = 0
    pairs = []
    delt = []
    DT = []
    
    for i in range(max(indexs)+1):
        delt.append(config_asc['deltaT'][i][i+1])
        delt.append(config_asc['deltaT'][i][i+1])
        delt.append(config_asc['deltaT'][i][i+1])
        delt.append(config_asc['deltaT'][i][i+1])
        for j in config_asc['V_range'][i].keys():
            A[(cnt*4):(cnt*4) + 4, i*4:j*4] = np.hstack([S]*(j-i))
            # print((cnt*4),(cnt*4) + 4, i*3,(j)*3, i, j)
            factor = config_asc['deltaT'][i][j] # Converting to displacement
            L.append(config_asc['V_range'][i][j]*factor)
            L.append(config_asc['V_azimuth'][i][j]*factor)
            L.append(config_des['V_range'][i][j]*factor)
            L.append(config_des['V_azimuth'][i][j]*factor)
            pairs.append(f'{i}_{j}')
            DT.extend([factor]*4)
            cnt += 1
    
    L = np.array(L)
    delt = np.array(delt)
    DT = np.array(DT)
    return A, L, mask, pairs, delt, DT


def create_Amat_sbas2(config_asc, config_des):
    indexs = config_asc['V_range'].keys()
    mask = (config_des['Mask']&config_asc['Mask'])
    n_pairs = config_asc['Num_pairs']
    # S = get_matS(config_asc, config_des)
    S = np.eye(4)
    A = np.zeros((4*n_pairs+4*(len(indexs)-1), 4*len(indexs)), dtype=float)
    L = []
    cnt = 0
    pairs = []
    delt = []
    DT = []
    
    for i in range(max(indexs)+1):
        delt.append(config_asc['deltaT'][i][i+1])
        delt.append(config_asc['deltaT'][i][i+1])
        delt.append(config_asc['deltaT'][i][i+1])
        delt.append(config_asc['deltaT'][i][i+1])
        for j in config_asc['V_range'][i].keys():
            A[(cnt*4):(cnt*4) + 4, i*4:j*4] = np.hstack([S]*(j-i))
            # print((cnt*4),(cnt*4) + 4, i*3,(j)*3, i, j)
            factor = config_asc['deltaT'][i][j] # Converting to displacement
            L.append(config_asc['V_range'][i][j]*factor)
            L.append(config_asc['V_azimuth'][i][j]*factor)
            L.append(config_des['V_range'][i][j]*factor)
            L.append(config_des['V_azimuth'][i][j]*factor)
            pairs.append(f'{i}_{j}')
            DT.extend([factor]*4)
            cnt += 1

    I3 = np.eye(4)
    lambda_reg = 0.1
    v_zero = np.zeros_like(L[0])
    for i in range(len(indexs)-1):
        factor1 = config_asc['deltaT'][i][i+1]
        A[(4*n_pairs)+(4*i):(4*n_pairs)+(4*(i+1)), i*4:(i+1)*4] = (I3/factor1)*lambda_reg
        
        factor2 = config_asc['deltaT'][i+1][i+2]
        A[(4*n_pairs)+(4*i):(4*n_pairs)+(4*(i+1)), (i+1)*4:(i+2)*4] = (-I3/factor2)*lambda_reg
        for n in range(4):
            L.append(v_zero)
            DT.append(config_asc['deltaT'][i][i+1])
    
    L = np.array(L)
    delt = np.array(delt)
    DT = np.array(DT)
    return A, L, mask, pairs, delt, DT


def svd_solve(a, b):
    [U, s, Vt] = la.svd(a, full_matrices=False)
    r = max(np.where(s >= 1e-12)[0])
    temp = np.dot(U[:, :r].T, b) / s[:r]
    return np.dot(Vt[:r, :].T, temp)


def compute_deltaA(config):
    deltaA = []
    consi = []
    idxs = [0]
    cnt = 0
    for i in range(len(config['deltaT'].keys())):
        for j in config['deltaT'][i].keys():
            cnt+=1
            if j==i+1:
                consi.extend([True]*4)
            else:
                consi.extend([False]*4)
            
            deltaA.append(config['deltaT'][i][j])
            deltaA.append(config['deltaT'][i][j])
            deltaA.append(config['deltaT'][i][j])
            deltaA.append(config['deltaT'][i][j])

        idxs.append(cnt)    

    deltaA = np.array(deltaA)
    consi = np.array(consi)
    return deltaA, consi


def compute_corr_vnev(A, L, mask, delt, deltaA):
    h, w = L.shape[1:]
    A_bar = np.linalg.inv(A.T@A)@A.T
    margin = 0
    V_nev = np.zeros((A.shape[1], h, w))
    mdrs = 0
    residues = np.zeros((A.shape[0], h, w))
    wgt = 1/((deltaA)**0.8)
    P = np.diag(wgt)
    L2 = np.zeros((A.shape[0], h, w))

    pbar = tqdm(range(margin, h-margin))
    for i in pbar:
        residue_sc = np.zeros_like(residues[:,i])
        for j in range(margin, w-margin):
            if mask[i,j]==False:
                continue
            # v_nev, res, rank, s = sp_lstsq(A, L[:,i,j], lapack_driver='gelsy', check_finite=False)
            # v_nev, res, rank = sq_lstsq(A_bar, L[:,i,j])
            # v_nev = svd_solve(A, L[:,i,j])
            # print(P.shape, A.shape, L[:,i,j].shape, 'pal')
            v_nev, res, rank, s = sp_lstsq(np.dot(P,A), np.dot(L[:,i,j],P), lapack_driver='gelsy', check_finite=False)
            # print(i, j)
            # print((L[:,i,j]/deltaA)[consi], v_nev/delt)
            residues[:,i,j] = (L[:,i,j] - np.dot(A, v_nev))/deltaA
            residue_sc[:,j] = residues[:,i,j]
            # print(residue_sc[:,j])
            
            L2[:,i,j] = np.dot(A, v_nev)
            V_nev[:,i,j] = v_nev/delt
        
        arrr = np.abs(residue_sc[np.logical_and(np.logical_not(np.isnan(residue_sc)),residue_sc!=0)])
        if arrr.shape[0]>0:
            mdrs = np.median(arrr)
            maxs = np.max(arrr)
            pbar.set_description(f'Median Residual: {mdrs:.02f} {maxs:.02f}')
            pbar.update()

    V_nev[np.isnan(V_nev)] = 0
    residues[np.isnan(residues)] = 0
    L2[np.isnan(L2)] = 0
    
    return V_nev, residues, L2


if __name__=="__main__":

    base_path_asc = '/DATA/test_stack_asc3/offset_tracking'
    csv_path_asc = '/DATA/test_stack_asc3/image_pairs.csv'

    base_path_des = '/DATA/test_stack_des3/offset_tracking'
    csv_path_des = '/DATA/test_stack_des3/image_pairs.csv'
    # os.environ['PROJ_LIB'] = '/home/pc/miniconda3/envs/isce2/share/proj/proj.db'

    ds = gdal.Open(f'{base_path_asc}/stack_1_2/velocity_res.tif')
    projs = ds.GetProjection()
    geo = ds.GetGeoTransform()
    vl = ds.GetRasterBand(1).ReadAsArray()
    H, W = vl.shape
    ds = None
    
    print("Loading offset tracking images")
    config_asc = read_timeseries(csv_path_asc, base_path_asc, skip_idx=[], xoff=0, yoff=0, xsize=W, ysize=H)
    config_asc = interpolate_missing(config_asc)
    config_des = read_timeseries(csv_path_des, base_path_des, skip_idx=[], xoff=0, yoff=0, xsize=W, ysize=H)
    config_des = interpolate_missing(config_des)
    
    print("Creating Design Matrix")
    A, L, mask, pairs, delt, DT = create_Amat_sbas(config_asc, config_des)
    deltaA, _ = compute_deltaA(config_asc)
    
    print("Performing SBAS correction")
    V_nev, residues, L_out = compute_corr_vnev(A, L, mask, delt, deltaA)
    
    print("Saving generated data")
    np.save('/DATA/sbas_out/V_nev_weight.npy', V_nev)
