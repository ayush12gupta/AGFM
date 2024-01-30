from shutil import copyfile, move # Utilities for copying and moving files
from osgeo import gdal, osr            # GDAL support for reading virtual files
import os, shutil                      # To create and remove directories
import numpy as np
# from rasterio.mask import mask
# import scipy
from glob import glob
import scipy.linalg as la
from datetime import datetime
from tqdm import tqdm
from scipy.linalg import lstsq as sp_lstsq
# from utils import 
from geogrid_autorift.util import numpy_array_to_raster, get_incang, get_flightang, get_matA_3D


def read_velocity(path):
    ds = gdal.Open(path)
    V_R_asc = ds.GetRasterBand(1).ReadAsArray()
    V_A_asc = ds.GetRasterBand(2).ReadAsArray()
    ds = None
    return V_R_asc, V_A_asc


def save_raster_outputs(save_dir, V_nev):
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
            if max(med-std, 0)>ab_resi[i]:
                cnt = max(med-std, 0)-ab_resi[i]
            elif med+std<=ab_resi[i]:
                cnt = ab_resi[i] - (med+std)
                
            # print(cnt)
            P2[i] = ((med/(med+cnt))**0.3)*P[i]
            # P2[i] = (2*np.abs(resi[i])/std)*P[i]
        else:
            P2[i] = 0

    # print(P2)
        
    return np.diag(P2)


def updateP(P, resi, sigma=None):
    ab_resi = np.abs(resi)
    sigma = np.median(np.abs(resi)) + 2*np.abs(resi).std()
    med = np.median(np.abs(resi))
    std = max(np.abs(resi).std(), 3)
    if sigma<5:
        return np.diag(P)
        
    P2 = P.copy()
    for i in range(P.shape[0]):
        # print(med-std, med+std, ab_resi[i])
        # print(max(med-std, 0)>ab_resi[i], (med+std)<ab_resi[i])
        if ((med+0.5*std)>ab_resi[i]):
            P2[i] = P[i]
        elif ((med+3*std)>ab_resi[i]):
            # P[i] = 1/((2*sigma/np.abs(resi[i])))*P[i]
            if med+0.5*std<ab_resi[i]:
                cnt = ab_resi[i] - (med+0.5*std)
                
            P2[i] = ((med/(med+cnt))**0.3)*P[i]
        else:
            P2[i] = 0
        
    return np.diag(P2)



def compute_3DVel(V_aoro, config_asc_path, config_des_path):
    '''
    V_aoro: [V_asc_rg, V_asc_az, V_des_rg, V_des_az]
    '''
    V_nev_all = []
    MD = []
    v_sbas = []

    N, _, H, W = V_aoro.shape
    for ni in range(N):
        stack = f'stack_{ni}_{ni+1}'
        print('Epoch:',ni, stack)
        # v_asc_path = f'{v_asc_dir}{stack}/velocity_res_crop.tif'
    #     v_des_path = f'{v_des_dir}{stack}/velocity_res_crop.tif'
    #     V_R_asc, V_A_asc = read_velocity(v_asc_path)
        incid_asc = np.zeros_like(V_aoro[ni,0])
        incid = (get_incang(config_asc_path))*(np.pi/180)
        incid_asc[:,:] = incid
        fl_ang_asc = (get_flightang(config_asc_path))*(np.pi/180)
        
    #     V_R_des, V_A_des = read_velocity(v_des_path)
        incid_des = np.zeros_like(V_aoro[ni,0])
        incid = (get_incang(config_des_path))*(np.pi/180)
        incid_des[:,:] = incid
        fl_ang_des = (get_flightang(config_des_path))*(np.pi/180)
        
        # nan = (np.isnan(V_aoro[ni,0]) | np.isnan(V_aoro[ni,0]))|((V_aoro[ni,0]==-32767)|(V_aoro[ni,0]==-32767))
        nan = np.isnan(V_aoro[ni,0])|((V_aoro[ni,0]==-32767)|(V_aoro[ni,0]==0))
        nan2 = np.isnan(V_aoro[ni,2])|((V_aoro[ni,2]==-32767)|(V_aoro[ni,2]==0))
        nan = nan|nan2
        
        margin = 0
        V_nev = np.zeros((3, H, W))
        residue = np.zeros((4, H, W))
        B_new = np.zeros((4, H, W))
        B = np.zeros((4, H, W))
        A = get_matA_3D(fl_ang_asc, incid_asc[0,2], fl_ang_des, incid_des[0,2])
        A_bar = np.linalg.inv(A.T@A)@A.T

        pbar = tqdm(range(margin, H-margin))
        mdrs, mnrs = 0, 0
        nn = 0
        for i in pbar:
            for j in range(margin, W-margin):
                if nan[i,j]==True:
                    continue
                    
                # A = get_matA(fl_ang_asc, incid_asc[i,j], fl_ang_des, incid_des[i,j])
                b = V_aoro[ni,:,i,j]
                P = np.eye(4)
                P = initP(np.diag(P), b)
                # print(i, j, 'IJ')
                for iter in range(10):
                    v_nev, res, rank, s = sp_lstsq(np.dot(P,A), np.dot(b,P), lapack_driver='gelsy', check_finite=False)
                    res = (b - np.dot(A, v_nev)).reshape(4)
                    P = updateP(np.diag(P), res)
                    # print(snr, 'SNR')
                    # print(v_nev, '3D')
                    # print(b, 'orig')
                    # print(np.dot(A, v_nev), 'new')
                    # print(res, 'resi')
                    # print(np.diag(P))
                    # print('++++++++++++++++++++++++++++++++++++++')
                    if np.median(np.abs(res))<1:
                        break
                
                # v_nev, res, rank, s = sp_lstsq(A, b, lapack_driver='gelsy', check_finite=False)
                # print(np.dot(A, v_nev))
                # print('------------------------------------------')
                
                residue[:,i,j] = (b - np.dot(A, v_nev)).reshape(4)
                B[:,i,j] = b.reshape(4)
                B_new[:,i,j] = np.dot(A, v_nev).reshape(4)
                V_nev[:,i,j] = v_nev.reshape(3)

            # mdrs = np.median(np.abs(residue[np.logical_and(np.logical_not(np.isnan(residue)),residue!=0)]))
            resid = residue[:,i,np.logical_not(nan[i])]
            arr = np.abs(resid)
            md = np.median(arr)
            mn = np.mean(arr)
            if not np.isnan(md):
                mdrs+=md
                mnrs+=mn
                nn = nn + 1
                pbar.set_description(f'Median Residual: {(mdrs/nn):.02f}| {(mnrs/nn):.02f}')
                pbar.update()

        V_nev_all.append(V_nev)
        v_sbas.append(V_aoro[ni,:,np.logical_not(nan)])
        # print(residue[:,np.logical_not(nan)].shape)
        MD.append(residue[:,np.logical_not(nan)])


    V_nev_all = np.array(V_nev_all)
    MD = np.array(MD)
    v_sbas = np.array(v_sbas)
    
    return V_nev_all


def compute_3DVel_single(V_aoro, config_asc_path, config_des_path):
    incid_asc = (get_incang(config_asc_path))*(np.pi/180)
    fl_ang_asc = (get_flightang(config_asc_path))*(np.pi/180)
    incid_des = (get_incang(config_des_path))*(np.pi/180)
    fl_ang_des = (get_flightang(config_des_path))*(np.pi/180)
    
    A = get_matA_3D(fl_ang_asc, incid_asc, fl_ang_des, incid_des)
    
    # b = V_aoro[ni,:,i,j]
    # Initialising weights
    P = np.eye(4)
    P = initP(np.diag(P), V_aoro)
    
    for iter in range(10):
        v_nev, res, rank, s = sp_lstsq(np.dot(P,A), np.dot(V_aoro,P), lapack_driver='gelsy', check_finite=False)
        res = (V_aoro - np.dot(A, v_nev)).reshape(4)
        P = updateP(np.diag(P), res)
        
        if np.median(np.abs(res))<1:
            break

    return v_nev.reshape(3), res


# def compute_3DVel_single_hermet(V_aoro, config_asc_path, config_des_path):
#     incid_asc = (get_incang(config_asc_path))*(np.pi/180)
#     fl_ang_asc = (get_flightang(config_asc_path))*(np.pi/180)
#     incid_des = (get_incang(config_des_path))*(np.pi/180)
#     fl_ang_des = (get_flightang(config_des_path))*(np.pi/180)
    
#     A = get_matA_3D(fl_ang_asc, incid_asc, fl_ang_des, incid_des)
#     A1, A2 = A[::2], A[1::2]
#     b1, b2 = V_aoro[::2], V_aoro[1::2]
#     P = np.eye(4)
#     P = initP(np.diag(P), V_aoro)
#     P1 = np.diag(np.diag(P)[::2])
#     P2 = np.diag(np.diag(P)[1::2])
#     for i in range(3):
#         N = (A1.T@P1@A1) + (A2.T@P2@A2)
#         U = (A1.T@P1@b1) + (A2.T@P2@b2)
#         v_nev, res, rank, s = sp_lstsq(N, U, lapack_driver='gelsy', check_finite=False)
#         # print('it', it, np.dot(A, v_nev), b)
#         # print(np.diag(P1)[0], np.diag(P2)[0], np.diag(P1)[1], np.diag(P2)[1])
#         resi1 = (b1 - np.dot(A1, v_nev)) #.reshape(4)    
#         sigma01 = (resi1.T@P1@resi1)/2
#         # sigma01 = (resi1.T@resi1)/2
        
#         resi2 = (b2 - np.dot(A2, v_nev)) #.reshape(4)
#         sigma02 = (resi2.T@P2@resi2)/2
#         P2 = (sigma01/sigma02)*P2
#         if sigma02==0:
#             break
        
#         if np.abs(sigma01-sigma02)<1:
#             break
    
#     res = (V_aoro - np.dot(A, v_nev))
#     return v_nev.reshape(3), res 


def compute_3DVel_single_hermet(V_aoro, A):
    A1, A2 = A[::2], A[1::2]
    b1, b2 = V_aoro[::2], V_aoro[1::2]
    P = np.eye(4)
    P = initP(np.diag(P), V_aoro)
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
        # sigma01 = (resi1.T@resi1)/2
        
        resi2 = (b2 - np.dot(A2, v_nev)) #.reshape(4)
        sigma02 = (resi2.T@P2@resi2)/2
        if sigma02==0:
            break
        
        P2 = (sigma01/sigma02)*P2
        
        if np.abs(sigma01-sigma02)<1:
            break
    
    res = (V_aoro - np.dot(A, v_nev))
    return v_nev.reshape(3), res 
    

if __name__=="__main__":
    ds = gdal.Open('/DATA/test_stack_asc3/offset_tracking/stack_0_3/velocity_res_crop.tif')
    projs = ds.GetProjection()
    geo = ds.GetGeoTransform()
    ds = None
    
    # Reading SBAS corrected velocities
    V_aoro = np.load('/DATA/sbas_out/V_nev_weight.npy')
    V_aoro = np.array([V_aoro[::4], V_aoro[1::4], V_aoro[2::4], V_aoro[3::4]]).transpose(1, 0, 2, 3)
    
    v_asc_dir = '/DATA/test_stack_asc3/offset_tracking/'
    v_des_dir = '/DATA/test_stack_des3/offset_tracking/'
    config_des_path = '/DATA/test_stack_des3/offset_tracking/testGeogrid.txt'
    config_asc_path = '/DATA/test_stack_asc3/offset_tracking/testGeogrid.txt'
    save_dir = '/DATA/3D_disp_sbas2'
    os.makedirs(save_dir, exist_ok=True)
    
    
    stacks_li = sorted(glob('/DATA/test_stack_des3/offset_tracking/stack*'))
    stacks_li = [stack.split('/')[-1] for stack in stacks_li]
    stacks1, stacks = [], []
    for stack in stacks_li:
        i, j = stack.split('_')[1:]
        if (int(j) - int(i))==1:
            stacks1.append(stack)
            
    for i in range(len(stacks1)):
        stacks.append(f'stack_{i}_{i+1}')
    
    # Performing 2D pairs to 3D conversion
    V_nev_all = compute_3DVel(V_aoro, config_asc_path, config_des_path)
    
    # Saving results to a directory
    save_raster_outputs(save_dir, V_nev_all)
    