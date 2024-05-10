from osgeo import gdal            # GDAL support for reading virtual files
import os, shutil, time                      # To create and remove directories
import multiprocessing as mp
import numpy as np
from glob import glob
import scipy.linalg as la
from tqdm import tqdm
from scipy.linalg import lstsq as sp_lstsq
from offset_tracking.util import numpy_array_to_raster, get_incang, get_flightang, get_matA_3D


def read_velocity(path):
    ds = gdal.Open(path)
    V_R_asc = ds.GetRasterBand(1).ReadAsArray()
    V_A_asc = ds.GetRasterBand(2).ReadAsArray()
    projs = ds.GetProjection()
    geo = ds.GetGeoTransform()
    ds = None
    return V_R_asc, V_A_asc, projs, geo


def save_raster_outputs(save_dir, V_nev, projs, geo):
    if os.path.exists(save_dir):
        shutil.rmtree(save_dir)
    os.mkdir(save_dir)
    for i in range(V_nev.shape[0]):
        ds = numpy_array_to_raster(f'{save_dir}/stack_{i}_{i+1}_corrected.tif', V_nev[i], projs, geo, nband=3)
        ds.FlushCache()
        ds = None


def initP(P, resi, sigma=None):
    ab_resi = np.abs(resi)
    sigma = np.median(np.abs(resi)) + 2*np.abs(resi).std()
    med = np.median(np.abs(resi))
    std = max(np.abs(resi).std(), 3)
    # print(resi, sigma, med, std)
    if sigma<5:
        return np.diag(P)
        
    P2 = P.copy()
    for i in range(P.shape[0]):
        if (max(med-std, 0)<ab_resi[i])&((med+std)>ab_resi[i]):
            P2[i] = P[i]
        elif (max(med-3*std, 0)<ab_resi[i])&((med+3*std)>ab_resi[i]):
            # P[i] = 1/((2*sigma/np.abs(resi[i])))*P[i]
            if max(med-std, 0)>=ab_resi[i]:
                cnt = max(med-std, 0)-ab_resi[i]
            elif med+std<=ab_resi[i]:
                cnt = ab_resi[i] - (med+std)
                
            P2[i] = ((med/(med+cnt))**0.3)*P[i]
            # P2[i] = (2*np.abs(resi[i])/std)*P[i]
        else:
            P2[i] = 0

    return np.diag(P2)


def compute_3DVel_single_hermet(V_aoro, A):
    A1, A2 = A[::2], A[1::2]
    b1, b2 = V_aoro[::2], V_aoro[1::2]
    P = np.eye(4)
    # P = initP(np.diag(P), V_aoro)
    P1 = np.diag(np.diag(P)[::2])
    P2 = np.diag(np.diag(P)[1::2])
    for i in range(3):
        N = (A1.T@P1@A1) + (A2.T@P2@A2)
        U = (A1.T@P1@b1) + (A2.T@P2@b2)
        v_nev, res, rank, s = sp_lstsq(N, U, lapack_driver='gelsy', check_finite=False)
        # print('it', it, np.dot(A, v_nev), b)
        # print(np.diag(P1)[0], np.diag(P2)[0], np.diag(P1)[1], np.diag(P2)[1])
        resi1 = (b1 - np.dot(A1, v_nev)) #.reshape(4)    
        sigma01 = (resi1.T@P1@resi1)/2
        
        resi2 = (b2 - np.dot(A2, v_nev)) #.reshape(4)
        sigma02 = (resi2.T@P2@resi2)/2
        if sigma02==0:
            break
        
        P2 = (sigma01/sigma02)*P2
        
        if np.abs(sigma01-sigma02)<1:
            break
    
    res = (V_aoro - np.dot(A, v_nev))
    return v_nev.reshape(3), res 


def unpacking_column(tup):
    jj, in_shape, A, v_aoro = tup
    
    v_nev = np.empty(in_shape[:-1],dtype=np.float32)
    v_nev.fill(np.nan)
    
    for ii in range(in_shape[1]):
        if np.isnan(v_aoro[:,ii]).any():
            continue
            
        v_nev[:,ii], _ = compute_3DVel_single_hermet(v_aoro[:,ii], A)

    return v_nev


def compute_3DVel_mp(V_aoro, config_asc_path, config_des_path):
    incid_asc = (get_incang(config_asc_path))*(np.pi/180)
    fl_ang_asc = (get_flightang(config_asc_path))*(np.pi/180)
    incid_des = (get_incang(config_des_path))*(np.pi/180)
    fl_ang_des = (get_flightang(config_des_path))*(np.pi/180)

    A = get_matA_3D(fl_ang_asc, incid_asc, fl_ang_des, incid_des)
    
    V_nev = np.empty((3, *V_aoro.shape[1:]),dtype=np.float32)
    V_nev.fill(np.nan)

    V_inp = V_aoro.copy()
    # V_inp[V_inp==0] = np.nan

    in_shape = V_nev.shape
    num_cores = min(mp.cpu_count(), 64)
    chunk_inputs = [(jj, in_shape, A, V_inp[:,:,jj]) for jj in range(in_shape[2])]

    with mp.Pool(processes=num_cores) as pool:
        V_nev = pool.map(unpacking_column, chunk_inputs)

    V_nev = np.transpose(np.array(V_nev), (1,2,0))
    
    return V_nev


def batch_inversion(id_master_slave, asc_path, des_path, out_dir):
    
    asc_imgs = glob(f'{asc_path}/offset_tracking/*/velocity.tif')
    des_imgs = glob(f'{des_path}/offset_tracking/*/velocity.tif')
    stacks = list(set([st.split('/')[-2] for st in asc_imgs]) & set([st.split('/')[-2] for st in des_imgs]))
    # assert len(stacks)==len(des_imgs)
    assert len(stacks)!=0
    
    # cwd = os.getcwd()
    os.makedirs(out_dir, exist_ok=True)

    start = time.time()
    for stack in tqdm(stacks):
        try:
            dt_m, dt_s = id_master_slave[np.where(id_master_slave[:,0]==stack)][0,1:]
        except:
            continue

        va_rg, va_az, projs, geo = read_velocity(f'{asc_path}/offset_tracking/{stack}/velocity.tif')
        vd_rg, vd_az, _, _ = read_velocity(f'{des_path}/offset_tracking/{stack}/velocity.tif')
        V_aoro = np.array([va_rg, va_az, vd_rg, vd_az])
        V_nev = compute_3DVel_mp(V_aoro, f'{asc_path}/offset_tracking/testGeogrid.txt', f'{des_path}/offset_tracking/testGeogrid.txt')

        ds = numpy_array_to_raster(f'{out_dir}/{dt_m}_{dt_s}_{stack[6:]}.tif', V_nev, projs, geo, nband=3)
        ds.FlushCache()
        ds = None

    print("Time Taken:", time.time()-start)
